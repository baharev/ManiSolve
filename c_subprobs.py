# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from collections import namedtuple
from copy import deepcopy
from functools import partial
from itertools import chain, groupby
import os.path

from nx import dfs_preorder_nodes
import nx
nx_reversed = nx.utils.reversed

from dag import create_dag, ntype, get_J_rowwise
from c_codegen import print_fwd_subgraph, print_forward_sweep, \
                      print_backward_sweep
from nl_interpreter import Range
from py3compat import ifilter, imap, irange, izip
from string import Template
from utils import StdoutHijack, clean


TMP_DIR  = '/tmp/stack_machine/' 
# Getting the absolute path of ./va27/
VA27_LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'va27')


# Assumptions: (1) each con / var is in one of the blocks, (2) the block ids are 
# contiguous and 1-based in the .nl, (3) there are equally many con & var blocks
# and (4) the blocks are in block lower Hessenberg form.

def main():
    clean(TMP_DIR)
    probs = ['cse', 'cse2', 'cse3', 'cse4', 'JacobsenTorn', 'mssTornDbg']
    h_max = 0
    def n_fixed_var_func(*_args):
        return 0
    for problem_name in probs:
        g, problem = create_dag(problem_name)
        generate_c_code(g, problem, h_max, n_fixed_var_func)


def generate_c_code(g, problem, h_max, n_fixed_var_func, so_suffix=''):
    # Sparsity pattern of the Jacobian row-wise, codegen is *not* taking problem
    Jrows = get_J_rowwise(problem)
    # All variable bounds
    lbs, ubs = get_var_bnds(problem)
    # node attributes into separate dictionaries, both for clarity and to make
    # codegen *not* to take the g as argument
    names  = {n: g.node[n]['name']  for n in g if 'name'  in g.node[n]}
    values = {n: g.node[n]['value'] for n in g if 'value' in g.node[n]}
    # Invariant: py_point[var_order[vi]] == vi
    # All functions dealing with py_point[] must take var_order as argument!
    var_order = fake_var_order(problem)
    # Invariant: py_residual[con_order[Ci]] == Ci
    con_order = fake_con_order(problem)
    fnames, problem_name = [], problem.name
    for index, con_ids, var_ids in gen_blocks(problem, h_max):
        state = setup(g, con_ids, var_ids, var_order, lbs, ubs)
        with StdoutHijack() as logger:
            codegen(state, con_ids, var_ids, con_order, var_order, index, Jrows, 
                                                names, values, n_fixed_var_func)
        c_code = logger.captured_text()
        # Write the C source file for the subproblem
        fname = problem_name + '_%d.c' % index
        fnames.append(fname)
        with open(TMP_DIR + fname, 'w') as f:
            f.write('// %s\n' % fname)
            f.write(c_code)
    # Generate and add wrapper code
    code = get_wrapper_code(len(fnames))
    #print(code)
    fname = problem_name + '.c'
    fnames.append(fname)
    with open(TMP_DIR + fname, 'w') as f:
        f.write(code)
    # Shell script to compile all C source files into a single .so
    compile_script(problem_name, ' '.join(fnames), h_max, so_suffix)

#-------------------------------------------------------------------------------

def setup(g, con_ids, var_ids, var_order, lbs, ubs):
    # con_ids, var_ids: new constraints and variables, introduced in this block.
    # node_order is needed to maintain the child node order in the subgraphs:
    node_order = {n: i for i, n in enumerate(g)}
    # Expression graph of the block
    g_sub = subgraph_of_deps(g, con_ids, node_order)
    # Hack: finding base and defined vars based on index>=n_vars and NOT ntype
    n_vars = len(var_order)
    # Referenced base and defined variables in the block
    base, defv = get_base_defined_vars(g_sub, n_vars)
    # prev_defv: defined variables NOT depending on any variable in var_ids
    prev_defv = get_prev_defvars_set(g_sub, var_ids, defv)
    # Expression graph of the previously seen defined variables
    g_prev_defv = subgraph_of_deps(g, prev_defv, node_order)
    # Re-write g_sub so that the previously seen defvars are input to the block:
    # Delete the in-edges of the pref_defv, then keep the deps of con_ids only.
    # We will later turn the ntype of these defvar-s to var too (in_defv).
    cut_at_prev_defv(g_sub, prev_defv, con_ids) # in-place reduction of g_sub!
    _base, in_defv = get_base_defined_vars(g_sub, n_vars)
    in_defv = [v for v in in_defv if v in prev_defv]
    # in_base: directly referenced, seen base vars
    in_base = set(base) - set(var_ids)
    in_base = [n for n in in_base if n in g_sub]
    # Sorting is not necessary but guarantees a linear read from py_point
    in_base = sorted(in_base, key=var_order.get)
    # Extract the bounds of the relevant variables
    indices = [int(vi[1:]) for vi in var_ids]
    lo = [lbs[i] for i in indices]
    up = [ubs[i] for i in indices]
    return in_base, in_defv, g_sub, g_prev_defv, lo, up

def codegen(state, con_ids, var_ids, con_order, var_order, index, Jrows, names, 
                                                      values, n_fixed_var_func):
    in_base, in_defv, g_sub, g_prev_defv, lo, up = state
    print_header(len(con_ids),   len(var_ids), n_fixed_var_func, \
                 len(con_order), len(var_order), index)
    print_var_bounds(names, var_ids, lo, up)
    # VA27 x and r (subproblem) to py_point and py_res (full problem)
    print_copy_back(con_ids, var_ids, con_order, var_order, names)
    # ... and the other way around
    print_copy_from(var_ids, var_order, names)
    # All functions dealing with py_point[] must take var_order as argument!
    #--- Function evaluation
    # Store the values of the previously seen base and defined variables.
    # The returned code will read the stored values from the global arrays at 
    # the beginning of the function / Jacobian evaluation.
    fix_input(in_base, var_order, in_defv, g_prev_defv, names, values)
    # Generate global array setup function
    setup_code(in_base, in_defv, index)
    # Start printing the forward evaluation code.
    forward_eval_header(in_base, in_defv, var_ids, index, names)
    # Turn the in_defv into input variables: g_sub.node is hacked from now on!
    turn_prev_defv_to_vars(g_sub, in_defv)
    # The actual function evaluation
    code2 = func_eval_code(g_sub, con_ids)
    print(code2)
    save_residuals(con_ids, names)
    # Forward evaluation done!
    print('}\n') # <- must not be in save_residuals
    #--- Function value checking
    # Only for debugging: checking the residual approx. 0.0 at the solution
    check_code(index, var_ids, var_order)
    #--- Jacobian
    # Start printing Jacobian computation
    jac_eval_header(in_base, in_defv, var_ids, index, names)
    code3 = jac_eval_code(g_sub, con_ids, in_base, in_defv, values)
    print(code3)
    store_results(con_ids, var_ids, Jrows, var_order, names)
    # Jacobian evaluation done!
    #--- Jacobian cross-checking, only for debugging
    jac_check_code(index, var_ids, con_ids, var_order)
    #--- Solve code, calling the VA27 solver
    # n_cons is needed: here start the artificial, variable fixing constraints
    print(VA27_CODE.substitute(index=index, n_cons=len(con_ids)))

#-------------------------------------------------------------------------------
# In the C code, we take (x_0:i-1, x_i) from Python, that is, all the seen base 
# variables and the new base variables x_i. We store the seen stuff into global 
# variables (arrays). Some of the variables of x_i can be fixed (optional), the
# corresponding arrays (indices and values) are also global.
# This group of functions read and write these global arrays.

def print_header(n_cons, n_vars, n_fixed_var_func, n_cons_full, n_vars_full, index):
    n_fixed_vars = n_fixed_var_func(n_cons, n_vars)
    m_actual = n_cons + n_fixed_vars
    row_padding = m_actual % 2 # M_ALLOCATED is even -> 16-byte aligned rows
    print(HEADER.substitute(n_cons=n_cons, n_fixed_vars=n_fixed_vars, 
                            row_padding=row_padding, n_vars=n_vars,
                            n_cons_full=n_cons_full, n_vars_full=n_vars_full,
                            index=index))

HEADER = Template( # file name will be inserted only later, at the top level
'''#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define BIG_M 1.0e20

// 1. Over-allocate rows so that variables can be fixed
//    later by adding additional w*(x_i-c_i)=0 constraints.
// 2. Furthermore, add padding for 16-byte (32-byte?) 
//    alignment of the rows of the J_buffer (J^T * J).

// Full problem sizes for the strides in the x_2D and r_2D arrays
#define N_VARS_FULL ($n_vars_full)
#define M_CONS_FULL ($n_cons_full)

// n_cons == len(con_ids), that is, the constraints in the subproblem
// (n_cons + n_fixed_vars + row_padding)
#define M_ALLOCATED ($n_cons + $n_fixed_vars + $row_padding)

// (n_cons + n_fixed_vars)
#define M_PLUS_FIXED ($n_cons + $n_fixed_vars)

// n_vars == len(var_ids), that is, the variables in the subproblem
#define N ($n_vars)

#define N_FIXED_VARS ($n_fixed_vars)

int asan_bug_workaround_$index[4] = {0}; // <-- ignore this

static   int indices_fixed[N_FIXED_VARS];

static double values_fixed[N_FIXED_VARS];

static double r_buffer[M_ALLOCATED];

static double J_buffer[N][M_ALLOCATED];  // store "transposed"


static int has_insane_values(const double* arr, int n) {
    for (int i=0; i<n; ++i) {
        const double x = arr[i];
        if (!((-BIG_M < x) && (x < BIG_M))) {
            return 1;
        }
    }
    return 0;
}

static void save_fixed_vars(int point_index, 
                            const int idx_fix[][N_FIXED_VARS], 
                            const double val_fix[][N_FIXED_VARS]) 
{
    if (N_FIXED_VARS == 0)
        return;

    const int* idx = idx_fix[point_index];
    const double* val = val_fix[point_index];
    
    for (int i=0; i<N_FIXED_VARS; ++i) {
        indices_fixed[i] = idx[i];
        values_fixed[i] = val[i];
    }
}

static void set_fixed_vars(double* __restrict__ va27_x) {

    assert(!has_insane_values(values_fixed, N_FIXED_VARS));

    for (int i=0; i<N_FIXED_VARS; ++i) {
        int k = indices_fixed[i];
        assert((0 <= k) && (k < N));
        va27_x[k] = values_fixed[i];
    }
}
''')

def fix_input(in_base, var_order, in_defv, g_prev_defv, names, values):
    fix_seen_base(in_base, var_order, names)
    fix_prev_defv(g_prev_defv, var_order, len(var_order), in_defv, values)
    # From the global arrays, copy back the values to v_i, it will be used in 
    # the function and Jacobian evaluation, use: fix_vi_in_func_or_jac_eval()

def fix_seen_base(in_base, var_order, names): # only directly referenced
    if not in_base:
        return
    print()
    print('static double in_base[%d];' % len(in_base))
    print()
    print('static void fix_seen_base(const double* __restrict__ py_point) {\n')
    fmt = '    in_base[{i}] = py_point[{k}];  // {v}, {name}'
    for i, v in enumerate(in_base):
        print(fmt.format(i=i, k=var_order[v], v=v, name=names[v]))
    print('}\n')

def fix_prev_defv(g_prev_defv, var_order, n_vars, in_defv, values):
    if not g_prev_defv:
        return
    print()
    print('static double in_defv[%d];' % len(in_defv))
    print()
    print('static void fix_defv(const double* __restrict__ py_point) {', end='')
    print(' // py_point may be unused due to an AMPL bug \n')
    # Figure out the base vars that we actually need for evaluation prev_defv, 
    # then recompute and store the previously seen defined vars.
    base, _defv = get_base_defined_vars(g_prev_defv, n_vars)
    # Fix the base vars necessary to evaluate the defvars
    fmt = '    const double {v} = py_point[{i}];  // {name}'
    for v in base:
        print(fmt.format(i=var_order[v], v=v, name=g_prev_defv.node[v]['name']))
    print()
    # Compute the defined variables
    print_fwd_subgraph(g_prev_defv, values)
    print()
    # Store the defined variables in a global array
    fmt = '    in_defv[{i}] = {v};  // {name}'
    for i, v in enumerate(in_defv):
        print(fmt.format(i=i, v=v, name=g_prev_defv.node[v]['name']))
    print('}\n')

def fix_vi_in_func_or_jac_eval(in_base, in_defv, names):
    code = [ ]
    fbase = '    const double {v} = in_base[{i}];  // {name}'
    for i, v in enumerate(in_base):
        code.append(fbase.format(i=i, v=v, name=names[v]))
    #
    if in_defv:
        code.append('')
    fdefv = '    const double {v} = in_defv[{i}];  // {name}'
    for i, v in enumerate(in_defv):
        code.append(fdefv.format(i=i, v=v, name=names[v]))
    return '\n'.join(code)

def setup_code(in_base, in_defv, index):
    decl = 'static int setup_{index}({used}const double* __restrict__ py_point)'
    used = '' if in_base or in_defv else '__attribute__((unused)) ' 
    print(decl.format(index=index, used=used) + ' {\n')
    if in_base:
        print('    fix_seen_base(py_point);')
        print('    if (has_insane_values(in_base, %d))' % len(in_base))
        print('        return 1;')
    if in_defv:
        print('    fix_defv(py_point);')
        print('    if (has_insane_values(in_defv, %d))' % len(in_defv))
        print('        return 1;')
    print('    return 0;')
    print('}\n')

def forward_eval_header(in_base, in_defv, var_ids, index, names):
    print()
    print('void func_eval_%d(const double* __restrict__ va27_point,\n'
          '                        double* __restrict__ res        ) {\n' % index)
    print(fix_vi_in_func_or_jac_eval(in_base, in_defv, names))
    assign_actual_base_vars(var_ids, names)
        
def assign_actual_base_vars(var_ids, names):
    print()
    fmt = '    const double {v} = va27_point[{i}];  // {name}'
    for i, v in enumerate(var_ids):
        print(fmt.format(i=i, v=v, name=names[v]))

def save_residuals(con_ids, names):
    print()
    for i, Ci in enumerate(con_ids):
        print('    res[{i}] = {Ci};  // {name}'.format(i=i, Ci=Ci, name=names[Ci]))

def print_var_bounds(names, var_ids, lo, up):
    vnames = ', '.join('%s %s' % (names[v], v) for v in var_ids) 
    print(VAR_BNDS.substitute(vnames=vnames, lbs=', '.join(lo), ubs=', '.join(up)))

VAR_BNDS = Template(
'''extern void randpoint(const double* lb, const double* ub, const int n, double* x);

static void randpoint_wrapper(double x[N]) {

    // $vnames
    static const double lb[N] = { $lbs };  
    static const double ub[N] = { $ubs };    

    randpoint(lb, ub, N, x);
}
''')

def print_copy_back(con_ids, var_ids, con_order, var_order, names):
    # Writing back the VA27 x to the py_point
    fmt = '    py_point[{k}] = x[{i}];  // {name}'
    point_map = '\n'.join(fmt.format(i=i, k=var_order[v], name=names[v])  
                                                 for i, v in enumerate(var_ids))
    # Writing back the VA27 r to the py_residual
    fmt = '    py_residual[{k}] = r[{i}];  // {name}'
    res_map = '\n'.join(fmt.format(i=i, k=con_order[c], name=names[c])  
                                                 for i, c in enumerate(con_ids))
    print(VA27_TO_PY.substitute(va27_x_to_py=point_map, va27_r_to_py=res_map))    

VA27_TO_PY = Template('''
static void copy_back_point(const double* __restrict__ x,
                                  double* __restrict__ py_point) 
{
$va27_x_to_py
}

static void copy_back_residuals(const double* __restrict__ r,
                                      double* __restrict__ py_residual) 
{
$va27_r_to_py
}
''')

def print_copy_from(var_ids, var_order, names):
    # Write the py_point to the VA27 x
    fmt = 'va27_point[{i}] = py_point[{k}];  // {name}'
    py_to_va27 = '\n    '.join(fmt.format(i=i, k=var_order[v], name=names[v]) 
                               for i, v in enumerate(var_ids))
    print(COPY_PY_TO_VA27.substitute(py_to_va27=py_to_va27))

COPY_PY_TO_VA27 = Template('''
static void copy_py_to_x(const double* __restrict__ py_point,
                               double* __restrict__ va27_point) 
{
$py_to_va27
}
''')

#-------------------------------------------------------------------------------

def jac_eval_header(in_base, in_defv, var_ids, index, names):
    print(JAC_EVAL_HEADER.substitute(index=index))
    print(fix_vi_in_func_or_jac_eval(in_base, in_defv, names))
    assign_actual_base_vars(var_ids, names)

JAC_EVAL_HEADER = Template('''
void jac_eval_$index(const double* __restrict__ va27_point,
                            double* __restrict__ res,
                            double jac[][M_ALLOCATED]) // store "transposed"
{
''')

#-------------------------------------------------------------------------------

def func_eval_code(g_sub, con_ids):
    with StdoutHijack() as logger:
        seen_def_vars = set()
        for Ci in con_ids:
            print_forward_sweep(g_sub, seen_def_vars, Ci)
    return logger.captured_text()

#-------------------------------------------------------------------------------

def jac_eval_code(g_sub, con_ids, in_base, in_defvars, values):
    input_vars = set(int(v[1:]) for v in chain(in_base, in_defvars))
    with StdoutHijack() as logger:
        seen_def_vars = set()
        for Ci in con_ids:
            print_forward_sweep(g_sub, seen_def_vars, Ci)
            print_backward_sweep(g_sub, Ci, values, input_vars=input_vars)
    return logger.captured_text()

def store_results(con_ids, var_ids, Jrows, var_order, names):
    save_residuals(con_ids, names)
    print()
    # Store the Jacobian of the block. The indices i_ and j_ are the row and 
    # column indices in this small Jacobian, but we will store it transposed.
    def index_of(n):
        return int(n[1:]) 
    var_j_ = {index_of(vi): j_ for j_, vi in enumerate(var_ids)}
    fmt = '    jac[{j_}][{i_}] = u_{i}_{j};' # Jacobian stored transposed
    # The goal of the loop is to fill out the i, j, i_, j_ in the fmt above.
    for i_, Ci in enumerate(con_ids):
        con_index = index_of(Ci)
        # Get those var indices that are in the block, then get their index j_ 
        # in the small Jacobian.
        for var_index in ifilter(lambda vi: vi in var_j_, Jrows[con_index]):
            j_ = var_j_[var_index]
            print(fmt.format(i=con_index, j=var_index, i_=i_, j_=j_))
    print('}\n')

#-------------------------------------------------------------------------------
# Helper functions, mainly used for computing the g_sub, g_prev_defv subgraphs,
# and the in_base and in_defv variables.

def get_base_defined_vars(g, n_vars):
    # Hack: collect all variables (base and defined) but NOT based on ntype
    var_deps = ifilter(lambda n: n[0]=='v', g)
    # sort by index
    var_deps = sorted(var_deps, key=lambda n: int(n[1:]))
    idx = find_index_of_first_def_var(var_deps, n_vars)
    base = var_deps[:idx]
    defv = var_deps[idx:]
    return base, defv

def find_index_of_first_def_var(var_list, n_vars):
    # Hack: find base and defined variables based on index >= n_vars, 
    # and NOT based on ntype. 
    def defined_var(v):
        assert v[0] == 'v', v 
        return int(v[1:]) >= n_vars
    gen_def_vars = (i for i, v in enumerate(var_list) if defined_var(v))
    return next(gen_def_vars, len(var_list))

def subgraph_of_deps(g, sources, node_order):
    deps = sorted(get_reachable_node_set(g, sources), key=node_order.get)
    return g.subgraph(deps)

def get_prev_defvars_set(g, var_ids, defv):
    node_set = set(var_ids)
    return {v for v in defv if none_reachable(g, v, node_set)}

def none_reachable(g, source, node_set):
    # Returns true if none of the nodes in node_set is reachable from source.
    with nx_reversed(g):
        reachable = dfs_preorder_nodes(g, source)
        return next((0 for n in reachable if n in node_set), 1)

def cut_at_prev_defv(g_sub, prev_defv, con_ids):
    g_sub.remove_edges_from([(n,v) for v in prev_defv for n in g_sub.pred[v]])
    deps = get_reachable_node_set(g_sub, con_ids)
    g_sub.remove_nodes_from([n for n in g_sub if n not in deps])

def turn_prev_defv_to_vars(g_sub, prev_defv):
    # g.subgraph(deps) does NOT copy node attributes and we are about to change
    # them here, therefore we have to make a deep copy of that dictionary.
    g_sub.node = deepcopy(g_sub.node)
    for v in prev_defv:
        g_sub.node[v]['kind'] = ntype.var

def get_reachable_node_set(g, sources):
    # Reverse the edges, and run a DFS from each source, remove duplicate nodes.
    with nx_reversed(g):
        # function, when invoked with a node n as an argument, it returns a 
        # generator of reachable nodes from n 
        reachable_from = partial(dfs_preorder_nodes, g)
        return set(chain.from_iterable(imap(reachable_from, sources)))

#-------------------------------------------------------------------------------
# Iteration logic: the way we extract the blocks from the .nl file, and the way 
# we iterate over the blocks, or sets of consecutive blocks.

def gen_blocks(problem, h_max):
    # yields: last block index, con_ids, var_ids
    assert h_max or _assert_block_lower_hessenberg_form(problem)
    ichain = chain.from_iterable
    #            t = block_index, con_ids, var_ids
    partition = [t for t in _blk_index_cons_vars(problem)]
    for slc in _gen_slices(len(partition), h_max):
        all_con_ids = list(ichain(con_ids for _, con_ids, _ in partition[slc]))
        all_var_ids = list(ichain(var_ids for _, _, var_ids in partition[slc]))
        yield slc.stop-1, all_con_ids, all_var_ids

Slices = namedtuple('Slices', 'seen  subp  new')

def gen_x_r_slices(problem, h_max):
    # Yields: x and r Slices, having seen, subp, and new attributes.
    partition = [ncons_nvars for ncons_nvars in _ncons_nvars_per_block(problem)] 
    ncons_seen, nvars_seen = 0, 0
    for slc in _gen_slices(len(partition), h_max):
        new_ncons, new_nvars = partition[slc.stop-1]
        ncons_seen += new_ncons
        nvars_seen += new_nvars
        ncons_subp = sum(n_cons for n_cons, _ in partition[slc])
        nvars_subp = sum(n_vars for _, n_vars in partition[slc])        
        r_seen = slice(                      0, ncons_seen)
        r_subp = slice(ncons_seen - ncons_subp, ncons_seen)
        r_new  = slice(ncons_seen - new_ncons , ncons_seen)
        x_seen = slice(                      0, nvars_seen)
        x_subp = slice(nvars_seen - nvars_subp, nvars_seen)
        x_new  = slice(nvars_seen - new_nvars , nvars_seen)
        yield Slices(x_seen, x_subp, x_new), Slices(r_seen, r_subp, r_new)

def _gen_slices(stop, lag):
    # gen_slices(4, 2) -> 0:1, 0:2, 0:3, 1:4
    for i in irange(stop):
        yield slice(max(0, i-lag), i+1)

def _blk_index_cons_vars(problem):
    # Yields: (0-based block index, [Ci], [vi]).
    for (cblk, cons), (vblk, vrs) in _gen_blocks(problem):
        assert cblk == vblk
        index = cblk-1
        yield index, ['C%d' % i for i, _ in cons], ['v%i' % i for i, _ in vrs]

def _ncons_nvars_per_block(problem):
    # Yields: (n_cons, n_vars) in the order of the blocks
    for (_cblk, cons), (_vblk, vrs) in _gen_blocks(problem):
        yield sum(1 for _ in cons), sum(1 for _ in vrs)

def _assert_block_lower_hessenberg_form(problem):
    # Either returns True or raises AssertionError
    Jrows = get_J_rowwise(problem)
    seen_vars = set()
    for (cblk, cons), (vblk, vrs) in _gen_blocks(problem):
        assert cblk == vblk
        var_set = {i for i, _ in vrs}
        deps =  set(chain.from_iterable(Jrows[i] for i, _ in cons))
        deps -= var_set
        deps -= seen_vars
        assert not deps, (problem.name, sorted(deps), cblk, cons, vrs)
        seen_vars |= var_set
    unseen_vars = set(irange(problem.nl_header.n_vars))
    unseen_vars -= seen_vars
    assert not unseen_vars, sorted(unseen_vars)
    return True

def _groupby(itr, keyf):
    # eager version of groupby
    return [(k, list(v)) for k, v in groupby(sorted(itr, key=keyf), key=keyf)]

def _gen_blocks(problem):
    # Generates: (cblk, cons), (vblk, vrs), where cons and vrs are generators:
    # (int id, int blk_id), and the block ids (cblk, vblk, blk_id) are 1 based.
    segs = problem.segments    
    # Get [(id, blk_id)] from the S segments
    con_blk_ids, var_blk_ids = segs.con_blocks, segs.var_blocks
    n_cons, n_vars = problem.nl_header.n_cons, problem.nl_header.n_vars 
    # Assumptions: each con / var is in one of the blocks, the block ids are 
    # contiguous and 1-based in the .nl, there are equally many con & var blocks
    # after padding with an empty row / col block at the ends if necessary
    assert len(con_blk_ids) == n_cons, (len(con_blk_ids), n_cons) 
    assert len(var_blk_ids) == n_vars, (len(var_blk_ids), n_vars)
    def blocks(iterable):
        # iterable: [(id, block_id)] 
        def by_block_id(tup):
            return tup[1]
        return _groupby(iterable, by_block_id)
    # constraints and variables grouped by blocks
    cblkid_cid, vblkid_vid = blocks(con_blk_ids), blocks(var_blk_ids)
    # pad with with zero rows or columns at the ends if necessary to have equal 
    # number of blocks
    if cblkid_cid[0][0] == 2: # First var block has no cons, add an empty one 
        cblkid_cid.insert(0, (1, []))
    if vblkid_vid[-1][0] == cblkid_cid[-1][0] - 1: # Last con block has no vars
        vblkid_vid.append((cblkid_cid[-1][0], []))
    cblks = {blk for blk, _ in cblkid_cid}
    assert sorted(cblks) == list(irange(1, len(cblks)+1)), sorted(cblks) 
    vblks = {blk for blk, _ in vblkid_vid}
    assert sorted(vblks) == list(irange(1, len(vblks)+1))
    assert len(cblks) == len(vblks), (len(cblks), len(vblks)) 
    return izip(cblkid_cid, vblkid_vid)

#-------------------------------------------------------------------------------
# TODO We should sort by the real variable order instead of just by block ids?
#      See also the iteration logic above!

def fake_var_order(problem):
    var_blocks = problem.segments.var_blocks
    n_vars = problem.nl_header.n_vars
    return _fake_order('v%d', var_blocks, n_vars)

def fake_con_order(problem):
    con_blocks = problem.segments.con_blocks
    n_cons = problem.nl_header.n_cons
    return _fake_order('C%d', con_blocks, n_cons)

def _fake_order(vi_or_Ci, blocks, expected_length):
    # blocks is a list: [(id, block_id)] 
    def by_block_id(tup):
        return tup[1]
    n_i = [(n,i) for i, (n, _blk) in enumerate(sorted(blocks, key=by_block_id))]
    order =  {vi_or_Ci % n : i for n, i in n_i}
    assert len(order)==expected_length, (len(order),expected_length,'S segments?')
    return order

def var_idx_order(problem):
    var_blocks = problem.segments.var_blocks
    n_vars = problem.nl_header.n_vars
    return _idx_order(var_blocks, n_vars)

def con_idx_order(problem):
    con_blocks = problem.segments.con_blocks
    n_cons = problem.nl_header.n_cons
    return _idx_order(con_blocks, n_cons)

def _idx_order(blocks, length):
    # order: n-th var/con in AMPL is in position i in the permuted problem
    # blocks is a list: [(id, block_id)] 
    def by_block_id(tup):
        return tup[1]
    order = [-1]*length
    for i, (n, _blk) in enumerate(sorted(blocks, key=by_block_id)):
        order[n] = i
    assert sorted(order) == list(range(length)), 'Missing S segments?'
    return order

def get_var_bnds(problem):
    var_bnds = problem.segments.var_bnds
    assert len(var_bnds) == problem.nl_header.n_vars
    assert all(kind == Range.lb_ub for kind, _bnds in var_bnds) # only l<=x<=u implemented
    lbs = [lb for _kind, ( lb, _ub) in var_bnds]
    ubs = [ub for _kind, (_lb,  ub) in var_bnds]
    return lbs, ubs

#-------------------------------------------------------------------------------

def check_code(index, var_ids, var_order):
    n_base = len(var_ids)
    fmt = 'va27_point[{i}] = py_point[{k}];'
    py_to_va27 = '\n    '.join(fmt.format(i=i, k=var_order[v]) 
                               for i, v in enumerate(var_ids))
    code = CHECK.substitute(index=index, n_base=n_base, py_to_va27=py_to_va27)
    print(code)

CHECK = Template('''
void check_$index(const double* py_point, double* py_residual) {

    setup_$index(py_point);

    assert(N == $n_base);

    double va27_point[N];
    
    $py_to_va27
    
    for (int i=0; i<M_ALLOCATED; ++i)
        r_buffer[i] = NAN;
    
    func_eval_$index(va27_point, r_buffer);
                 
    copy_back_residuals(r_buffer, py_residual);
}
''')


def jac_check_code(index, var_ids, con_ids, var_order):
    fmt = 'va27_point[{i}] = py_point[{k}];'
    py_to_va27 = '\n    '.join(fmt.format(i=i, k=var_order[v]) 
                               for i, v in enumerate(var_ids))
    code = JAC_CHECK.substitute(index=index, n_cons=len(con_ids), 
                                py_to_va27=py_to_va27)
    print(code)

JAC_CHECK = Template('''
void jac_check_$index(const double* py_point, 
                      double* py_residual, 
                      double* py_2d_array)
{
    setup_$index(py_point);

    memset(J_buffer, 0, M_ALLOCATED*N*(sizeof J_buffer[0][0]));
    
    double va27_point[N];    
    
    $py_to_va27
    
    for (int i=0; i<M_ALLOCATED; ++i)
        r_buffer[i] = NAN;
    
    jac_eval_$index(va27_point, r_buffer, J_buffer);
    
    copy_back_residuals(r_buffer, py_residual);    
    
    // J_buffer is over-allocated and transposed, hence this copying below.
    
    double (*jac)[N] = (double (*)[N]) py_2d_array;
    
    for (int i=0; i<$n_cons; ++i)
        for (int j=0; j<N; ++j)
            jac[i][j] = J_buffer[j][i];
    
}

//------------------------------------------------------------------------------
''')

#-------------------------------------------------------------------------------
# These functions are concerned with creating the .so file. The top level 
# wrapper functions are necessary to execute the function and jacobian 
# evaluation of the individual subproblems (these functions are accessed by 
# index, instead of name). A shell script is created to compile the code later.

# Name mangling rules, also called by the c_sub_check.py

def compile_sh_path(problem_name, h_max, so_suffix):
    return TMP_DIR + problem_name + '_h_%d%s.sh' % (h_max, so_suffix)

def get_so_name(problem_name, h_max, so_suffix):
    return problem_name + '_h_%d%s.so' % (h_max, so_suffix)

def compile_script(problem_name, fnames, h_max, so_suffix):
    so_name = get_so_name(problem_name, h_max, so_suffix)
    with open(compile_sh_path(problem_name, h_max, so_suffix), 'w') as f:
        f.write(SHELL_SCRIPT.substitute(fnames=fnames, VA27_LIB=VA27_LIB, 
                                        so_name=so_name))

SHELL_SCRIPT = Template('''
set -e
rm -f *.o

# Release build:
gcc -c -Ofast -march=native -std=c99 -fPIC -Wall -Wno-unused-variable $fnames
gfortran -O3 -shared *.o -L$VA27_LIB -lva27 -o ${so_name}

# Debug build:
#clang -c -ggdb3 -fsanitize=address -fno-omit-frame-pointer \
#      -O0 -std=c99 -fPIC -Wall -Wextra -Wno-unused-variable \
#      -Wno-unused-parameter $fnames
#clang -ggdb3 -fsanitize=address -shared -fsanitize=address -shared-libasan \
#      -O0 -shared *.o -L$VA27_LIB -lva27 -lgfortran -o ${so_name}
''')

def get_decls_names(decl_fmt, n_subproblems):
    name_fmt = decl_fmt.split()[2]    # assumes extern void f(...
    name_fmt = name_fmt.split('(')[0] # f(const ... -> f
    decls  = (decl_fmt % i for i in irange(n_subproblems))
    names  = (name_fmt % i for i in irange(n_subproblems))
    return '\n'.join(decls), ',\n    '.join(names)

def get_wrapper_code(n_subproblems):
    check_fmt = 'extern void check_%d(const double* , double* );'
    jac_fmt   = 'extern void jac_check_%d(const double*, double*, double*);'
    setup_fmt = 'extern void setup_%d(const double* );'
    solve_fmt = '''extern void solve_%d(double* x_2d_array, 
                     double* res_2d_array,
                     const int n_points,
                     const int* idx_fixed_2d,
                     const double* val_fixed_2d,
                     const int n_cols_fixed_2d,
                     const double tolerance, // 1.0e-8, >= EPS[i]^2
                     const int n_trials,
                     unsigned int seed,
                     const int iprint);\n'''
    check_decls,     check_names     = get_decls_names(check_fmt, n_subproblems)
    jac_check_decls, jac_check_names = get_decls_names(jac_fmt,   n_subproblems)
    setup_decls,     setup_names     = get_decls_names(setup_fmt, n_subproblems)
    solve_decls,     solve_names     = get_decls_names(solve_fmt, n_subproblems)
    return WRAPPER.substitute(check_decls=check_decls, check_names=check_names,
                              jac_check_decls=jac_check_decls,
                              jac_check_names=jac_check_names,  
                              setup_decls=setup_decls, setup_names=setup_names,
                              solve_decls=solve_decls, solve_names=solve_names,
                              n_subproblems=n_subproblems)

WRAPPER = Template('''#include <assert.h>

//------------------------------------------------------------------------------

$check_decls

typedef void (*check_func_ptr)(const double* , double* );

static check_func_ptr check[$n_subproblems] = { 
    $check_names 
};

//------------------------------------------------------------------------------

/*
$setup_decls

typedef void (*setup_func_ptr)(const double* );

static setup_func_ptr setup[$n_subproblems] = {
    $setup_names
};
*/

//------------------------------------------------------------------------------

$jac_check_decls

typedef void (*jac_check_func_ptr)(const double* , double* , double* );

static jac_check_func_ptr jac_check[$n_subproblems] = { 
    $jac_check_names 
};

//------------------------------------------------------------------------------

$solve_decls

typedef void (*solve_func_ptr)(double* x_2d_array, 
                               double* res_2d_array,
                               const int n_points,
                               const int* idx_fixed_2d,
                               const double* val_fixed_2d,
                               const int n_cols_fixed_2d,
                               const double tolerance, // 1.0e-8, >= EPS[i]^2
                               const int n_trials,
                               unsigned int seed,
                               const int iprint);

static solve_func_ptr solve_subproblem[$n_subproblems] = { 
    $solve_names 
};

//------------------------------------------------------------------------------

void evaluate(int index, const double* py_point, double* residuals) {
    assert(index>=0 && index<$n_subproblems);
    check[index](py_point, residuals);
}

void jacobian_evaluation(int index, 
                         const double* py_point, 
                         double* residuals, 
                         double* py_2d_array) 
{    
    // py_2d_array is the array storing the Jacobian as a dense matrix
    assert(index>=0 && index<$n_subproblems);
    jac_check[index](py_point, residuals, py_2d_array);
}

void solve(int index,
           double* x_2d_array, 
           double* res_2d_array,
           const int n_points,
           const int* idx_fixed_2d,
           const double* val_fixed_2d,
           const int n_cols_fixed_2d,
           const double tolerance, // 1.0e-8, >= EPS[i]^2
           const int n_trials,
           unsigned int seed,
           const int iprint)
{
    assert(index>=0 && index<$n_subproblems);
    solve_subproblem[index](x_2d_array, res_2d_array, n_points, idx_fixed_2d, 
                            val_fixed_2d, n_cols_fixed_2d, tolerance, n_trials, 
                            seed, iprint);
}

''')

#-------------------------------------------------------------------------------

VA27_CODE = Template(
'''// Maps: x -> F(x) (= r), while storing r and J in r_buffer and J_buffer

static void resid(int* m, int* n, double* __restrict__ x, 
                  double* __restrict__ r, int* IFL) 
{
    // Callback function for the FORTRAN solver
    
    assert(*n == N);
    assert(*m == M_PLUS_FIXED);
    
    memset(r_buffer, 0, M_ALLOCATED *     (sizeof r_buffer[0])   );
    memset(J_buffer, 0, M_ALLOCATED * N * (sizeof J_buffer[0][0]));
    
    jac_eval_$index(x, r_buffer, J_buffer);
    
    if (has_insane_values(r_buffer, M_PLUS_FIXED)) {
        *IFL = 1;
        return;
    }
    
    const double* const J_ptr = J_buffer[0];
    
    // Careful: M_ALLOCATED*N elements have to be checked (not M_PLUS_FIXED*N)!
    
    if (has_insane_values(J_ptr, M_ALLOCATED*N)) {
        *IFL = 1;
        return;
    }
    
    // IFL already set to 0 by the caller.
    
    // Now append the fake constraints, essentially fixing
    // the specified variables.
    
    assert(!has_insane_values(values_fixed, N_FIXED_VARS));
    
    const double W = 100.0;
    
    for (int i=0; i<N_FIXED_VARS; ++i) {
        int k = indices_fixed[i];
        assert((0 <= k) && (k < N));
        r_buffer[$n_cons + i] = W*(x[k] - values_fixed[i]);
        J_buffer[k][$n_cons + i] = W; 
    }
    
    memcpy(r, r_buffer, M_PLUS_FIXED*(sizeof r_buffer[0]));
}

static void lsq(int* m, int* n, 
                __attribute__((unused)) double* x, 
                double* __restrict__ r, 
                double A[][N], 
                double* __restrict__ v) 
{   
    assert(*n == N);
    assert(*m == M_PLUS_FIXED);
    
    memcpy(r, r_buffer, M_PLUS_FIXED*(sizeof r_buffer[0]));

    for (int i=0; i<N; ++i) {
        double sum = 0.0;
        double* J_i = J_buffer[i];
        for (int k=0; k<M_PLUS_FIXED; ++k)
            sum += r[k] * J_i[k];
        v[i] = sum;
    }

    for (int i=0; i<N; ++i) {
        double* J_i = J_buffer[i];
        for (int j=0; j<=i; ++j) {  // only need the lower triangular part
            double sum = 0.0;
            double* J_j = J_buffer[j];
            for (int k=0; k<M_PLUS_FIXED; ++k)
                sum += J_i[k] * J_j[k];
            A[j][i] = sum;          // transpose due to Fortran
        }
    }
}

static int is_aligned(const void* pointer, size_t byte_count) { 
    return ((uintptr_t) pointer) % byte_count == 0;
}

typedef void (*RESID)(int* , int* , double* , double* , int* );

typedef void (*LSQ)(int* , int* , double* , double* , double[][N], double* );

//void va27ad_(RESID, LSQ, m, n, x, r, SS, A, D, EPS, IPRINT, MAXFN, MODE, W);
void va27ad_(RESID, LSQ, int*, int*, double*, double*, double*, double[][N], 
             double*, double*, int*, int*, int*, double*);

extern void set_seed(unsigned int seed);

void solve_$index(double* x_2d_array, 
             double* res_2d_array,
             const int n_points,
             const int* idx_fixed_2d,
             const double* val_fixed_2d,
             const int n_cols_fixed_2d,
             const double tolerance, // 1.0e-8, >= EPS[i]^2
             const int n_trials_signed,
             unsigned int seed,
             const int iprint)
{
    int n = N;
    int m = M_PLUS_FIXED;
    double X[N];
    double R[M_PLUS_FIXED];
    double SS;
    double A[N][N];    
    double D[N];    
    double EPS[N];
    int IPRINT = iprint;
    int MAXFN = 500;
    int MODE = 1;
    double W[4*N+M_PLUS_FIXED];
    const int use_first_as_starting_point = (n_trials_signed > 0) ? 0 : 1;
    const int n_trials = (n_trials_signed > 0) ? n_trials_signed : -n_trials_signed; 
    
    assert(n_trials >= 1);
    
    assert(is_aligned(r_buffer, 16));
    assert(is_aligned(J_buffer[0], 16));
#if N > 1 // otherwise we get warnings from clang
    assert(is_aligned(J_buffer[1], 16));
#endif
    assert( ( N_FIXED_VARS &&  idx_fixed_2d &&  val_fixed_2d &&  n_cols_fixed_2d) ||
            (!N_FIXED_VARS && !idx_fixed_2d && !val_fixed_2d && !n_cols_fixed_2d) );
    assert(n_cols_fixed_2d == N_FIXED_VARS);
    assert(N_FIXED_VARS < N && "No degrees of freedom left?");
    assert(N <= M_PLUS_FIXED && "Trying to solve an under-determined problem?");
    assert(M_PLUS_FIXED <= M_ALLOCATED);
    assert(N <= N_VARS_FULL && "Subproblem is at most the original problem");
    assert( (M_PLUS_FIXED - N_FIXED_VARS) <= M_CONS_FULL && 
            "Subproblem is at most the original problem" );  
    
    for (int i=0; i<N; ++i)
        EPS[i] = 1.0e-7;
    
    if (seed)
        set_seed(seed);
    
    //--------------------------------------------------------------------------
    
    double (*py_x)[N_VARS_FULL] = (double (*)[N_VARS_FULL]) x_2d_array;
    
    double (*py_res)[M_CONS_FULL] = (double (*)[M_CONS_FULL]) res_2d_array;

    const int (*idx_fix)[N_FIXED_VARS] = 
                                     (const int (*)[N_FIXED_VARS]) idx_fixed_2d; 
    
    const double (*val_fix)[N_FIXED_VARS] = 
                                  (const double (*)[N_FIXED_VARS]) val_fixed_2d; 

    //--------------------------------------------------------------------------
    
    for (int i=0; i<n_points; ++i) {
    
        double* py_point = py_x[i];
        
        int has_insane = setup_$index(py_point);
        
        if (has_insane) {
            fprintf(stderr, "Insane value(s) among the seen variables!\\n");
            fprintf(stderr, "Skipping point %d\\n", i);
            continue;
        }
        
        save_fixed_vars(i, idx_fix, val_fix);
        
        for (int k=1; k<=n_trials; ++k) {
            
            if ((k == 1) && use_first_as_starting_point) {
            
                copy_py_to_x(py_point, X);
            
                assert(!has_insane_values(X, N));
            }
            else {
                randpoint_wrapper(X);
            }            
            
            set_fixed_vars(X);
            
            /* Turn this into a macro?
            printf("Initial point: [");
            for (int ii=0; ii<N; ++ii)
                printf(" %.4f ", X[ii]);
            printf("]\\n");
            */
            
            // If the solver fails at the initial point, SS remains uninitialized, 
            // and we may not recognize that the solver has failed, and copy back
            // garbage as the "solution". 
            
            SS = BIG_M;
            
            va27ad_(resid, lsq, &m, &n, X, R, &SS, A, D, EPS, &IPRINT, &MAXFN, &MODE, W);
            
            double* py_residual = py_res[i];
            
            //  SS <= N * EPS[i]^2 roughly, so tolerance >= EPS[i]^2  
            if (SS/N < tolerance) {

                //printf("SS/N = %g\\n", SS/N);

                copy_back_point(X, py_point);
                 
                copy_back_residuals(R, py_residual);
                
                break;
            }
            
            // py_point and py_residual is not touched unless the solver succeeds!
        }
    }
}

''')


if __name__ == '__main__':
    main()
