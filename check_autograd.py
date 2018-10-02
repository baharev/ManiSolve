# Copyright (C) 2018 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from os.path import isfile
from six import iteritems, exec_
from dag import create_dag, get_J_rowwise
from dag_to_py import print_node as print_fwd_node
from utils import StdoutHijack, serialize, deserialize

def main():
    g, problem = create_dag('eco9')
    numeric_diff(g, problem)

PREAMBLE = \
'''
import autograd.numpy as np
from autograd import jacobian
from autograd.numpy import exp, log

# name: {name}, nodes: {n_nodes}, edges: {n_edges}
'''

POSTAMBLE = \
'''
jac = jacobian(f)
J   = jac(np.array([{x}], dtype=np.double)) 
'''

def numeric_diff(g, problem):
    n_vars = problem.nl_header.n_vars
    n_cons = problem.nl_header.n_cons
    #
    preamble = PREAMBLE.format(name=problem.name, n_nodes=g.number_of_nodes(), 
                               n_edges=g.number_of_edges())
    #
    def_f = 'def f(x):'
    assign = '\n'.join('    v{i} = x[{i}]'.format(i=i) for i in range(n_vars))
    #
    with StdoutHijack() as logger:
        for n, d in iteritems(g.node):
            print_fwd_node(g, n, d)
    code = logger.captured_text()
    code = '\n'.join('    ' + line for line in code.splitlines())
    #
    retval = ', '.join('C%d' % i for i in range(n_cons))
    retval = '    return np.array([%s])' % retval
    #
    sols = problem.solutions
    postamble = POSTAMBLE.format(x=', '.join(sols[0]))
    #
    return '\n'.join((preamble, def_f, assign, '\n', code, retval, postamble))

def get_numeric_jacobian(problem_name):
    # Always return the Jacobian as list of (i, j, J_ij)
    jacfile = 'data/' + problem_name + '.jac.pkl.gz'
    if isfile(jacfile):
        print('Using cached Jacobian for', problem_name)
        return deserialize(jacfile)
    #
    print('Computing then saving the Jacobian for', problem_name)
    g, problem = create_dag(problem_name)
    code = numeric_diff(g, problem)
    globals_ = {}
    try:
        exec_(code, globals_)
    except:
        print('===============================================================')
        print(code)
        print('===============================================================')
        raise
    J = globals_['J']
    #
    Jrows = get_J_rowwise(problem)
    jac = [(i, j, J[i,j]) for i, cols in enumerate(Jrows) for j in cols]
    serialize(jac, 'data/' + problem.name + '.jac.pkl.gz')
    return jac

if __name__ == '__main__':
    main()
