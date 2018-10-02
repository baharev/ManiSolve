# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from functools import partial
from itertools import count

from six import iteritems, itervalues
from nx import OrderedDiGraph as nx_OrderedDiGraph, dfs_preorder_nodes

from nl_interpreter import get_problem, interpret, Range, get_J_rowwise
from py3compat import imap, izip

#-------------------------------------------------------------------------------
# Hopefully, we can remove this hack one day when nx.OrderedDigraph is ready.
# Note that the code here won't have to be changed; we will be able to use
# nx.OrderedDiGraph directly. Currently, we have to override the subgraph()
# method. Networkx at:
# ali@X230:~/nx/networkx/networkx$ git log -n 1
# commit e43a7b08bd5ab2640b5a9c3350ed5355bdb82c65
# Date:   Fri Jun 23 17:57:33 2017 +0200

class OrderedDiGraph(nx_OrderedDiGraph):
    
    #@profile
    def subgraph(self, nbunch):
        bunch = self.nbunch_iter(nbunch)
        # create new graph and copy subgraph into it
        H = self.__class__()
        # copy node and attribute dictionaries
        for n in bunch:
            H.node[n]=self.node[n]
        # namespace shortcuts for speed
        H_succ=H.succ
        H_pred=H.pred
        self_succ=self.succ
        self_pred=self.pred
        # add nodes
        for n in H:
            H_succ[n]=H.adjlist_inner_dict_factory()
            H_pred[n]=H.adjlist_inner_dict_factory()
        # add successors
        for u in H_succ:
            Hnbrs=H_succ[u]
            for v,datadict in self_succ[u].items():
                if v in H_succ:
                    Hnbrs[v]=datadict
        # add predecessors
        for u in H_pred:
            Hnbrs=H_pred[u]
            for v,datadict in self_pred[u].items():
                if v in H_pred:
                    Hnbrs[v]=datadict
        H.graph=self.graph
        return H

#-------------------------------------------------------------------------------

class ntype:
    # these have AMPL names and AMPL indices; var and con have bounds
    var    = 'var'
    defvar = 'defvar'
    con    = 'con'
    # ---
    num    = 'num'
    # An operator not in the .nl file, but corresponding to the linear part
    lin_comb = 'lin_comb'
    # operators
    plus     = 'plus'
    minus    = 'minus'
    multiply = 'mult'  
    divide   = 'div'   
    negate   = 'neg'
    sumlist  = 'sum'   
    exp      = 'exp'   
    log      = 'log'   
    power    = 'pow'   


OPERATOR = {
    'plus':  ntype.plus,
    'minus': ntype.minus,
    'mult':  ntype.multiply,
    'div':   ntype.divide,
    'neg':   ntype.negate,
    'sum':   ntype.sumlist,
    'exp':   ntype.exp,
    'log':   ntype.log,
    'pow':   ntype.power,
    # non-AMPL operator
    'lin_comb': ntype.lin_comb, 
}

OPS = set(itervalues(OPERATOR))

#-------------------------------------------------------------------------------

def main():    
    #probname = 'tunnelDiodesSum'
    probname = 'mss10TornScaled'
    g, _problem = create_dag(probname)
    print(g, g.number_of_nodes(), g.number_of_edges())
    for n, d in iteritems(g.node):
        print(n, d)
    #
    gen_ops = (n for n, d in iteritems(g.node) if d['kind'] in OPS)
    for op in gen_ops:
        print(op, '=', g.node[op]['kind'], ' '.join(g.predecessors(op)))

def create_dag(problem_name):
    p = get_problem(problem_name)
    g = OrderedDiGraph()
    g.graph['name'] = problem_name
    add_all_vars(g, p.col_names)
    cons, def_vars = p.segments.cons, p.segments.def_vars
    C_names, V_names = p.row_names, p.def_var_names
    t_start = p.nl_header.n_vars + len(p.def_var_names) 
    t_counter, n_counter = count(start=t_start), count()
    for kind, index in p.segments.eval_order:
        if kind == 'C':
            add_constraint(g, index, cons, t_counter, n_counter, C_names)
        elif kind == 'V':
            add_defined_variable(g, index, def_vars, t_counter, n_counter, V_names)
    #
    Jrows = get_J_rowwise(p)
    check_sparsity_pattern(g, Jrows, p.nl_header.n_vars, p.nl_header.n_nzeros, 
                                                         p.segments.con_jacobian)
    return g, p

def add_constraint(g, index, cons, t_counter, n_counter, C_names):
    nl_part, lin_part, con_range = cons[index]
    nl_last, lin_last = add(g, nl_part, lin_part, t_counter, n_counter)
    rng = add_range(g, con_range, n_counter)
    last = sum_up(g, t_counter, nl_last, lin_last, rng)
    add_last(g, last, 'C%d' % index, ntype.con, C_names[index])

def add_defined_variable(g, index, def_vars, t_counter, n_counter, V_names):
    nl_part, lin_part = def_vars[index]
    nl_last, lin_last = add(g, nl_part, lin_part, t_counter, n_counter)
    last = sum_up(g, t_counter, nl_last, lin_last)
    add_last(g, last, 'v%d' % index, ntype.defvar,  V_names[index])

def add_last(g, last, node_id, kind, name):
    add_node(g, node_id, kind, name=name)
    add_edge(g, last, node_id)    

def sum_up(g, t_counter, nl_last, lin_last, rng=None):
    not_none = tuple(n for n in (nl_last, lin_last, rng) if n)
    if len(not_none) == 1:
        return not_none[0]
    return add_with_children(g, 't%d' % next(t_counter), ntype.sumlist, not_none)

def add_range(g, con_range, n_counter):
    kind = con_range[0]
    assert kind == Range.eq, kind
    (value,) = con_range[1]
    if value == '0':
        return
    # con = rhs - range -> rhs + (-range); flipping the sign:
    value = '-' + value if value[0] != '-' else value[1:]
    return add_num_node(g, n_counter, value)

def add(g, nl_part, lin_part,  t_counter, n_counter):
    nl_last  = add_nonlinear(g, nl_part, t_counter, n_counter)
    lin_last = add_linear(g, lin_part, t_counter)
    assert nl_last or lin_last
    return nl_last, lin_last 

def add_linear(g, lin_part, t_counter):
    indices, coefficients = lin_part
    nzero_coeffs = tuple(c for c in coefficients if c != '0')
    if not nzero_coeffs:
        return
    children = ['v%d' % i for i, c in izip(indices, coefficients) if c != '0']
    return add_with_children(g, 't%d' % next(t_counter), ntype.lin_comb, 
                                                 children, coeffs=nzero_coeffs)

def add_nonlinear(g, nl_part, t_counter, n_counter):
    body = list(interpret(nl_part, t_counter))
    if body == ['n0']:
        return    
    for lhs, op_name, args in body[:-1]:
        args = add_numbers(g, args, n_counter)
        add_with_children(g, lhs, OPERATOR[op_name], args)
    return add_if_num(g, n_counter, body[-1]) # the last node, value of nl_part    

def add_numbers(g, args, n_counter):
    # replace n-1.23 with n42 where 42 is a unique id, else just return arg
    return tuple(imap(partial(add_if_num, g, n_counter), args))

def add_if_num(g, n_counter, arg):
    return add_num_node(g, n_counter, arg[1:]) if arg[0] == 'n' else arg

def add_num_node(g, n_counter, value):
    return add_node(g, 'n%d' % next(n_counter), ntype.num, value=value)

def add_with_children(g, node_id, kind, children, **kwargs):
    # Adds the node_id, the children must already be in g
    add_node(g, node_id, kind, **kwargs)
    for child in children:
        add_edge(g, child, node_id)
    return node_id

def add_all_vars(g, col_names):
    for i, name in enumerate(col_names):
        add_node(g, 'v%d' % i, ntype.var, name=name)

#-------------------------------------------------------------------------------
# Only these functions add nodes / edges to g through g.add_*

def add_node(g, node_id, kind, **kwargs):
    assert node_id not in g, node_id
    g.add_node(node_id, kind=kind, **kwargs)
    return node_id

def add_edge(g, src, dst):
    assert src in g, src
    assert dst in g, dst
    g.add_edge(src, dst)

#-------------------------------------------------------------------------------

def check_sparsity_pattern(g, Jrows, n_vars, n_zeros, k_segment):
    # Cross-check with the k_segment, compute column lengths from rowwise
    col_counts  = get_col_counts(Jrows, n_vars)
    col_lengths = get_intermediate_col_lengths(col_counts)
    assert col_lengths[-1] == n_zeros
    assert k_segment == col_lengths[:-1]
    # Cross-check against DFS from C_i, reaching the base variables
    g.reverse(copy=False)
    for i, cols in enumerate(Jrows):
        check_base_vars(g, n_vars, i, cols)
    g.reverse(copy=False)

def check_base_vars(g, n_vars, i, cols):
    # Used to be  g.node[n]['index'] instead of int(n[1:])
    vrs = [int(n[1:]) for n in dfs_preorder_nodes(g, source='C%d' % i) if n[0] == 'v']
    indices = sorted(index for index in vrs if index < n_vars)
    assert indices == cols, (indices, cols)

def get_col_counts(Jrows, n_vars):
    col_counts = [0]*n_vars
    for cols in Jrows:
        for j in cols:
            col_counts[j] += 1
    return col_counts

def get_intermediate_col_lengths(col_count):
    lengths, accum = [], 0
    for count in col_count:
        accum += count
        lengths.append(accum)
    return lengths


if __name__ == '__main__':
    main()
