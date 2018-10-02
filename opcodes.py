# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from itertools import repeat
from six import itervalues
from py3compat import izip
from utils import as_pairs, duplicates

# ARITY: {opcode: arity}  
# NAME:  {opcode: name}  
# OPERATORS: all the possible operators in human readable format

__all__ = ('ARITY', 'NAME', 'OPERATORS')
 
#-------------------------------------------------------------------------------

__unary = (
    'o13', 'floor',
    'o14', 'ceil',
    'o15', 'abs',
    
    'o16', 'neg',
    'o34', 'not',
    'o37', 'tanh',
    
    'o38', 'tan',
    'o39', 'sqrt',
    'o40', 'sinh',
    
    'o41', 'sin',
    'o42', 'log10',
    'o43', 'log',
    
    'o44', 'exp',
    'o45', 'cosh',
    'o46', 'cos',
    
    'o47', 'atanh',
    'o49', 'atan',
    'o50', 'asinh',
    
    'o51', 'asin',
    'o52', 'acosh',
    'o53', 'acos',
)

__binary = (
    'o0', 'plus',
    'o1', 'minus',
    'o2', 'mult',
    
    'o3', 'div',
    'o4', 'rem',
    'o5', 'pow',
    
    'o6',  'less',
    'o20', 'or',
    'o21', 'and',
    
    'o22', 'lt',
    'o23', 'le',
    'o24', 'eq',
    
    'o28', 'ge',
    'o29', 'gt',
    'o30', 'ne',
    
    'o48', 'atan2',
    'o55', 'intdiv',
    'o56', 'precision',
    
    'o57', 'round',
    'o58', 'trunc',
    'o73', 'iff',
)

__n_ary = (
    'o11', 'min',
    'o12', 'max',
    'o54', 'sum',
    
    'o59', 'count',
    'o60', 'numberof',
    'o61', 'numberofs',
    
    'o70', 'n_ary_and',
    'o71', 'n_ary_or',
    'o74', 'alldiff',
    
    
    'o35', 'if',
    'o65', 'ifs',
    'o72', 'implies',
)

ARITY    = {}

def __add_arities(operators, arity):
    opcodes = operators[::2]    
    dups = duplicates(opcodes)
    assert not dups, dups
    already_added = set(opcodes) & set(ARITY)
    assert not already_added, already_added
    ARITY.update(izip(opcodes, repeat(arity)))

__add_arities(__unary,  1)
__add_arities(__binary, 2)
__add_arities(__n_ary,  None)

NAME = {}

def __add_names(operators):
    seen_names = set(itervalues(NAME))
    for opcode, pretty_name in as_pairs(operators):
        assert opcode not in NAME, (opcode, pretty_name)
        assert pretty_name not in seen_names, (opcode, pretty_name)
        seen_names.add(pretty_name)
        NAME[opcode] = pretty_name

__add_names(__unary)
__add_names(__binary)
__add_names(__n_ary)

assert set(NAME) == set(ARITY)

OPERATORS = set(itervalues(NAME))



