# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from collections import OrderedDict
from os import listdir
from toolz import partitionby
from py3compat import imap
from utils import as_pairs, get_lines, warning

# Except nl_interpreter.py (obviously), and check.py (which does not have a 
# Problem object), no other module should use other functions then problem_names

__all__ = ['problem_names', 'mss_ignored']

def problem_names(ignore=[]):
    testproblems = {f[:-4] for f in listdir('data') if f.endswith('.mod')}
    testproblems -= set(ignore)
    return sorted(testproblems, key=lambda s: s.upper())

def mss_ignored():
    ignored = sorted(f[:-4] for f in listdir('data') if f.startswith('mss') \
                     and f.endswith('_B.mod') and int(f[3:5]) > 35)
    warning('Ignoring the MSS family,', len(ignored), 'files')
    return ignored

def nl_fname(probname):
    return 'data/' + probname + '.nl'

def col_fname(probname):
    return 'data/' + probname + '.col'

def row_fname(probname):
    return 'data/' + probname + '.row'

def zero_fname(probname):
    return 'data/' + probname + '.zero'

#-------------------------------------------------------------------------------

def get_solutions(probname, colnames):
    # varnames come from the .zero file
    varnames, old_solutions = _read_old_solutions(probname)
    assert sorted(varnames) == sorted(colnames)
    # permute the old solution vectors into the order as in the colnames
    return _permute_solution_vectors(old_solutions, colnames, varnames)    

def _permute_solution_vectors(old_solutions, colnames, varnames):
    colname_index = {name: i for i, name in enumerate(varnames)}
    solutions = []
    for sol in old_solutions:
        v = [sol[colname_index[name]] for name in colnames]
        solutions.append(v)
    return solutions

def _read_old_solutions(probname):
    # grab the lines
    with open(zero_fname(probname), 'r') as f:
        lines = [l.strip() for l in f]
    # first line gives the dimension
    dim = int(lines[0])
    # chunks are separated by empty lines
    chunks = list(partitionby(lambda l: not l, lines[1:]))
    gen_chunks = as_pairs(chunks)
    # first chunk gives the variable names in order
    blank, varnames = next(gen_chunks)
    assert blank, blank
    assert dim == len(varnames), (dim, varnames)
    # the remaining chunks give the solution vectors
    solutions = list(_gen_solutions(gen_chunks, dim))
    return varnames, solutions

def _gen_solutions(gen_chunks, dim):
    for blank, data in gen_chunks:
        assert blank, blank
        assert dim == len(data), data        
        yield data

#-------------------------------------------------------------------------------

def get_col_names(probname):
    return get_lines(col_fname(probname))

def get_row_names(probname):
    return get_lines(row_fname(probname))

def get_def_var_names(probname):
    lines = get_lines(nl_fname(probname))
    def to_index_name(line):
        if line[0] == 'V':
            # "Vi j k #name" -> i as int
            index = int(line[1:].split(None, 1)[0])
            # "Vi j k #name" -> name
            name  = line.rsplit('#', 1)[1].strip() # Needs: option nl_comments 1;
            return index, name
        return None
    return OrderedDict(i_name for i_name in imap(to_index_name, lines) if i_name)
