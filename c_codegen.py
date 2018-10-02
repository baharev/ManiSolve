# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division

from nx import dfs_preorder_nodes, topological_sort

from dag_to_py import OPERATOR as FWD_OPERATOR, lin_comb as fwd_lin_comb
from dag import ntype
from py3compat import izip
from utils import StdoutHijack


def print_fwd_subgraph(g_sub, values):
    # The values are necessary only due to a bug in AMPL: Lambda[3,3,1] = 1,
    # so the constant propagation is not performed properly
    _get_value.values = values
    for n in topological_sort(g_sub):
        print_fwd_node(g_sub, n, g_sub.node[n])
    _get_value.values = None

def print_forward_sweep(g, seen_def_vars, Ci):
    print()
    print('    // Forward sweep on constraint {}, {}'.format(Ci, g.node[Ci]['name']))
    print()
    deps = get_deps(g, Ci)
    deps = remove_seen_def_vars(g, seen_def_vars, deps)
    for n in topological_sort(g.subgraph(deps)):
        print_fwd_node(g, n, g.node[n], seen_def_vars=seen_def_vars)

def remove_seen_def_vars(g, seen_def_vars, deps):
    seen = [n for n in deps if n in seen_def_vars]
    remove = set(seen)
    for s in seen:
        remove.update(get_deps(g, s))
    return [n for n in deps if n not in remove]

def print_backward_sweep(g, Ci, values, input_vars=None):
    print('\n    // Backward sweep on constraint {}'.format(Ci))
    _get_value.values = values
    # _assign is hacked so that input variables are not differentiated 
    _assign.input_vars = input_vars if input_vars is not None else set()
    # Jacobian is unused
    seen = set() # Tells whether to use = or += with u_i_j
    # Slight inefficiency: Constraint dependencies were already computed in 
    # print_forward_sweep. We throw it away and recompute deps for clarity.
    deps = get_deps(g, Ci)
    bwd_eval_order = list(topological_sort(g.subgraph(deps)))
    for n in reversed(bwd_eval_order):
        __print_backward_node(g, n, g.node[n], seen)
    _get_value.values  = None
    _assign.input_vars = None

def get_deps(g, source):
    g.reverse(copy=False)
    deps = list(dfs_preorder_nodes(g, source=source))
    g.reverse(copy=False)
    return deps

#-------------------------------------------------------------------------------

def __print_backward_node(g, n, d, seen):
    kind = d['kind']
    op_func = OPERATOR.get(kind, None)
    # operators
    if op_func is not None:
        fwd_in_comment(g, n, d)
        args = tuple(pred for pred in g.predecessors(n))
        op_func(g, n, seen, *args)
    # linear combination
    elif kind == ntype.lin_comb:
        fwd_in_comment(g, n, d)
        lin_comb(g, n, seen)
    # constraint
    elif kind == ntype.con:
        fwd_in_comment(g, n, d)
        start_new_constraint(g, n, d, seen)
    # defined variable
    elif kind == ntype.defvar:
        print()  # It is like fwd_in_comment:
        print('    // defined variable', n, d['name'])
        defined_var(g, n, seen)
    # variable
    elif kind == ntype.var:
        pass
    # number
    elif kind == ntype.num:
        pass # That's OK
    # should never happen
    else:
        raise AssertionError(str((n, d)))

def fwd_in_comment(g, n, d):
    # Only serves for debugging
    with StdoutHijack() as logger:
        print_fwd_node(g, n, d)
        lines = logger.captured_text().splitlines()
    lines = [l for l in lines if l]
    lines[0] = '\n    // ' + lines[0].strip()
    print('\n'.join(lines))

def get_u(node_id):
    assert get_u.con_idx is not None
    first_char = node_id[0]
    assert first_char == 't' or first_char == 'v', node_id 
    return 'u_' + get_u.con_idx + '_' + node_id[1:]

get_u.con_idx = None  # This is set locally, must NOT be set on the top level

def start_new_constraint(g, n, d, seen):
    seen.clear()
    get_u.con_idx = n[1:]
    constraint(g, n, seen)    

def _s(node_id):
    first_char = node_id[0]
    if first_char == 'v' or first_char == 't':
        return node_id
    value = _get_value(node_id)
    return value if value[0] != '-' else '(%s)' % value

def _get_value(node_id):
    assert node_id[0] == 'n', node_id
    return _get_value.values[node_id]

_get_value.values = None # Set by print_backward_sweep or by print_fwd_subgraph  

def _assign(u_k, seen, *args):
    # For input vars v_j we fix u_i_j to 0. We nevertheless put into comment 
    # the code that would be generated if v_j was not fixed.
    # Get column index corresponding to u_k:
    j = int(u_k.rsplit('_', 1)[1])  # 'u_i_j' -> j as int
    prefix = '    // ' if j in _assign.input_vars else '    '
    assign = '%s +=' if u_k in seen else 'double %s ='
    print((prefix + assign) % u_k, *args, end=';\n')
    #
    if j in _assign.input_vars and u_k not in seen: 
            print('    const double %s = 0.0;' % u_k)
    seen.add(u_k)

_assign.input_vars = None  # It will be set by print_backward_sweep

def constraint(g, n, seen):
    # C_k = t_i
    # u_i = 1.0
    (last, ) = tuple(g.predecessors(n))
    _assign(get_u(last), seen, '1.0')

def defined_var(g, n, seen):
    # v_i = t_k
    # u_k (+)= u_i
    (t_k, ) = tuple(g.predecessors(n))
    if t_k[0] != 'n':
        u_k = get_u(t_k)
        _assign(u_k, seen, get_u(n))

def plus(g, n, seen, x, y):
    sumlist(g, n, seen, x, y)

def sumlist(g, n, seen, *args):
    # t_i = t_r + ... + t_z
    # u_r (+)= u_i
    # ...
    # u_z (+)= u_i
    u_i = get_u(n)
    for arg in args:
        if arg[0] != 'n':
            _assign(get_u(arg), seen, u_i)

def minus(g, n, seen, *args):
    t_r, t_s = args
    # t_i = t_r - t_s
    # u_r (+)=  u_i
    # u_s (+)= -u_i
    u_i = get_u(n)
    if t_r[0] != 'n':
        _assign(get_u(t_r), seen, u_i)
    if t_s[0] != 'n':
        _assign(get_u(t_s), seen, '-', u_i)

def log(g, n, seen, x):
    # t_i = log(t_k)
    # u_k (+)= (1/t_k)*u_i
    u_k = get_u(x)
    u_i = get_u(n)
    _assign(u_k, seen, '(1.0/%s) * %s' % (_s(x), u_i))

def power(g, n, seen, x, y):
    # t_i   = pow(base, power)
    # u_k (+)= power*pow(base, power-1)*u_i
    assert y[0] == 'n', (n, x, y)
    u_i = get_u(n)
    u_k = get_u(x)
    pwer = g.node[y]['value'] 
    _assign(u_k, seen,'{pwer}*pow({x}, {pwer}-1)*{u_i}'.format(pwer=pwer,
                                                                  x=x, u_i=u_i))

def divide(g, n, seen, x, y):
    # t_i = t_r/t_s 
    # u_r (+)=  (1/t_s)       * u_i
    # u_s (+)= -(1/t_s) * t_i * u_i
    u_i = get_u(n)
    if x[0] != 'n':
        u_r = get_u(x)
        _assign(u_r, seen, '(1.0/{t_s})*{u_i}'.format(t_s=_s(y), u_i=u_i))
    if y[0] != 'n':
        u_s = get_u(y)
        _assign(u_s, seen,'-(1.0/{t_s})*{t_i}*{u_i}'.format(t_s=_s(y), 
                                                                t_i=n, u_i=u_i))

def negate(g, n, seen, x):
    # t_i = -t_s
    # u_s (+)= -u_i
    u_s = get_u(x)
    u_i = get_u(n)
    _assign(u_s, seen, '-', u_i)

def exp(g, n, seen, x):
    # t_i = exp(t_r)
    # u_r = t_i * u_i
    u_i = get_u(n)
    u_r = get_u(x)
    _assign(u_r, seen, '%s * %s' % (n, u_i))

def multiply(g, n, seen, x, y): 
    # t_i = t_r*t_s 
    # u_r (+)= t_s * u_i
    # u_s (+)= t_r * u_i
    u_i = get_u(n)
    _multiply(x, y, u_i, seen)
    _multiply(y, x, u_i, seen)
    
def _multiply(t_r, t_s, u_i, seen):
    if t_r[0] != 'n':
        # u_r (+)= t_s * u_i
        u_r = get_u(t_r)
        _assign(u_r, seen, _s(t_s), '*', u_i)

def lin_comb(g, n, seen):
    # t_i = sum c_k t_k 
    # for each k: 
    #   u_k (+)= c_k * u_i  
    args   = tuple(g.predecessors(n))
    coeffs = g.node[n]['coeffs']
    assert len(args) == len(coeffs)
    u_i = get_u(n)
    assert u_i in seen, u_i
    for t_k, c_k in izip(args, coeffs):
        if t_k[0] == 'n':
            continue
        u_k = get_u(t_k)
        if c_k == '1.0':
            rhs = u_i
        elif c_k == '-1.0':
            rhs = '-' + u_i
        elif c_k == '0.0':
            rhs = '0.0'
        else:
            rhs = '%s * %s' % ('(%s)' % c_k if c_k[0] == '-' else c_k,  u_i)
        _assign(u_k, seen, rhs)

OPERATOR = {
    'plus':  plus,
    'minus': minus,
    'mult':  multiply,
    'div':   divide,
    'neg':   negate,
    'sum':   sumlist,
    'exp':   exp,
    'log':   log,
    'pow':   power,
}

#===============================================================================

def fwd_assign(n, *args):
    print('    const double', n, '=', *args, end=';\n')

# TODO Keep in sync with dag_to_py.py
def print_fwd_node(g, n, d, seen_def_vars=None):
    seen_def_vars = seen_def_vars if seen_def_vars is not None else set()
    kind = d['kind']
    op_func = FWD_OPERATOR.get(kind, None)
    # operators
    if op_func is not None:
        args = tuple((pred, g.node[pred]) for pred in g.predecessors(n))
        rhs = op_func(*args)
        fwd_assign(n, rhs)
    # linear combination
    elif kind == ntype.lin_comb:
        fwd_assign(n, fwd_lin_comb(g, n))
    # constraint
    elif kind == ntype.con:
        (last, ) = tuple(g.predecessors(n))
        fwd_assign(n, last, ';  // ', d['name'])
    # defined variable
    elif kind == ntype.defvar:
        (last, ) = tuple(g.predecessors(n))
        fwd_assign(n, _s(last))
        print('    // defined variable:', d['name'])
        seen_def_vars.add(n)
    # variable
    elif kind == ntype.var:
        pass # Write variable assignments first?
    # number
    elif kind == ntype.num:
        pass # That's OK
    # should never happen
    else:
        raise AssertionError(str((n, d)))

