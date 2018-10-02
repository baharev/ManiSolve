# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
import nx #@UnusedImport
from six import iteritems, exec_
from auxfiles import problem_names, mss_ignored
from dag import create_dag, ntype
from py3compat import izip
from utils import StdoutHijack

#-------------------------------------------------------------------------------

def _str(x):
    n, d = x
    kind = d['kind']
    if kind != ntype.num:
        return n
    value = d['value']
    return '(%s)' % value  if value[0] == '-' else value

def plus(x, y):
    return '%s + %s' % (_str(x), _str(y))

def minus(x, y):
    return '%s - %s' % (_str(x), _str(y))

def multiply(x, y):
    return '%s * %s' % (_str(x), _str(y))

def divide(x, y):
    return '%s / %s' % (_str(x), _str(y))

def negate(x):
    return '-%s' % _str(x)

def sumlist(*args):
    return ' + '.join(_str(x) for x in args)

def exp(x):
    return 'exp(%s)' % _str(x)

def log(x):
    return 'log(%s)' % _str(x)

def power(x, y):
    return 'pow(%s, %s)' % (_str(x), _str(y)) 

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

def lin_comb(g, n):
    args   = tuple(g.predecessors(n))
    coeffs = g.node[n]['coeffs']
    assert len(args) == len(coeffs)
    terms = []
    for v_i, c in izip(args, coeffs):
        if c == '1':
            v = ' + %s' % v_i
        elif c == '-1':
            v = ' - %s' % v_i
        elif c == '0':
            v = ''
        else:
            v = ' + %s*%s' % ('(%s)' % c if c[0] == '-' else c,  v_i)
        terms.append(v)
    PLUS_SIGN = ' + '
    if terms[0].startswith(PLUS_SIGN):
        terms[0] = terms[0][len(PLUS_SIGN):]
    return ''.join(terms)

#-------------------------------------------------------------------------------

def main():
    testproblems = problem_names(ignore=['binary_nl_file'] + mss_ignored())
    for probname in testproblems:
        print('---------------------------------')
        check(probname)
    print('Checked', len(testproblems), 'problems')

#-------------------------------------------------------------------------------

def check(probname):
    # the code for evaluation
    with StdoutHijack() as logger:
        _g, problem = print_dag(probname)
    code = logger.captured_text()    
    # generate and execute the python script for each solution vector
    for i, sol in enumerate(problem.solutions):
        print('---------------------------------')
        print(probname, 'solution #%d' % (i+1), '\n')
        # prepend prelude, append prologue to the evaluation code 
        # prelude sets v_i to sol, prologue gets the residuals
        with StdoutHijack() as logger:
            print_script(code, sol, problem.nl_header.n_cons)
        script = logger.captured_text()
        globals_ = {}
        exec_(script, globals_)
        res = globals_['res'] # see get_string_residual
        max_viol = max(map(abs, res))
        assert max_viol < 1.0e-6, max_viol

def print_script(code, sol, res_dim):
    # Prelude
    print('from __future__ import print_function')
    print('from math import exp, log')
    print()
    print('# Solution vector')
    print(get_string_setting_vs(sol))
    print()
    # The actual code for evaluation
    print(code)
    # Prologue
    print()
    print(get_string_residual(res_dim))

def get_string_setting_vs(sol):
    return '\n'.join('v%d = %s' % (i, v) for i, v in enumerate(sol))

def get_string_residual(dim):
    res = 'res = [ %s ]' % ', '.join('C%d' % i for i in range(dim))
    return '\n'.join((res, 'print(max(map(abs, res)))'))

#-------------------------------------------------------------------------------

def print_dag(probname):
    g, problem = create_dag(probname)
    print('#', g, g.number_of_nodes(), g.number_of_edges())
    for n, d in iteritems(g.node):
        print_node(g, n, d)
    return g, problem

# TODO keep in sync with reverse_C.py
def print_node(g, n, d, seen_def_vars=None):
    seen_def_vars = seen_def_vars if seen_def_vars is not None else set()
    kind = d['kind']
    op_func = OPERATOR.get(kind, None)
    # operators
    if op_func is not None:
        args = tuple((pred, g.node[pred]) for pred in g.predecessors(n))
        rhs = op_func(*args)
        print(n, '=', rhs)
    # linear combination
    elif kind == ntype.lin_comb:
        print(n, '=', lin_comb(g, n))
    # constraint
    elif kind == ntype.con:
        (last, ) = tuple(g.predecessors(n))
        print(n, '=', last, ' #', d['name'], '\n')
    # defined variable
    elif kind == ntype.defvar:
        (last, ) = tuple(g.predecessors(n))
        print(n, '=', _str((last, g.node[last])))
        print('# defined variable:', d['name'], '\n')
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

if __name__ == '__main__':
    main()

