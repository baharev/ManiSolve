# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from collections import namedtuple, OrderedDict
from itertools import chain
from toolz import partitionby
from auxfiles import get_col_names, get_row_names, get_def_var_names, nl_fname,\
                     get_solutions
from opcodes import ARITY, NAME
from py3compat import imap, irange, izip
from utils import warning, as_pairs, get_lines

Problem = namedtuple('Problem', '''name  nl_header  segments  col_names  
                                   def_var_names  row_names  solutions''')

Header = namedtuple('Header', 'n_vars  n_cons  n_nzeros  n_objs')

#-------------------------------------------------------------------------------

def get_nl_file_content(full_nl_file_name):
    lines = get_lines(full_nl_file_name)
    lines = remove_comments(lines)
    header, content = split_to_header_and_content(lines)
    nl_header = parse_nl_header(header) 
    return nl_header, get_segments(nl_header, content)    

def parse_nl_file(probname):
    return get_nl_file_content(nl_fname(probname))

def remove_comments(lines):
    # 'v1  #x[6]'  ->  'v1'
    return list(imap(lambda l: l.rsplit('#', 1)[0].rstrip(), lines))

#-------------------------------------------------------------------------------
# Stuff that mainly concerns the header
#
# 0:  g3 1 1 0 # problem silly
# 1:   9 4 3 1 1 2 # vars, constraints, objectives, ranges, eqns, lcons
# 2:   3 3 # nonlinear constraints, objectives
# 3:   0 0 # network constraints: nonlinear, linear
# 4:   7 8 5 # nonlinear vars in constraints, objectives, both
# 5:   0 1 0 1 # linear network variables; functions; arith, flags
# 6:   0 0 0 0 0 # discrete variables: binary, integer, nonlinear (b,c,o)
# 7:   17 5 # nonzeros in Jacobian, gradients
# 8:   5 4 # max name lengths: constraints, variables
# 9:   2 0 0 3 1 # common exprs: b,c,o,c1,o1
# The lines below the header have a non-blank first character

def split_to_header_and_content(lines):
    if lines[0][0] != 'g':
        raise IOError('ASCII .nl file is expected.\nPlease generate the'
                      ' .nl file with "write gname;" in the AMPL .mod file\n'
                      'or as "ampl -ogname name.mod" in the command line.')
    parts = list(partitionby(lambda l: l[0]==' ', lines))
    header, content = list(chain(parts[0], parts[1])), parts[2] 
    return header, content

def parse_nl_header(header):
    assert len(header) == 10, 'See expected format above in comment'
    header = list(imap(to_ints, header[1:])) # <- discarding the first line!
    n_vars, n_cons, n_objs = header[0][:3]
    nzeros = header[6][0]
    if any(header[5]):
        warning('The model has discrete variables but we treat them as continuous!')
    return Header(n_vars=n_vars, n_cons=n_cons, n_nzeros=nzeros, n_objs=n_objs)

def to_ints(line):
    # '1 2' -> [1, 2]
    entries = line.split()
    return list(imap(int, entries))

#-------------------------------------------------------------------------------
# F  imported function description
# S  suffix values
# V  defined variable definition (must precede V,C,L,O segments where used)
# C  algebraic constraint body
# L  logical constraint expression
# O  objective function
# d  dual initial guess
# x  primal initial guess
# r  bounds on algebraic constraint bodies ("ranges")
# b  bounds on variable
# k  Jacobian column counts (must precede all J segments)
# J  Jacobian sparsity, linear terms
# G  Gradient sparsity, linear terms

SEGMENTS = {'F', 'S', 'V', 'C', 'L', 'O', 'd', 'x', 'r', 'b', 'k', 'J', 'G'}

# FIXME Start documenting the shape of the members!
Segments = namedtuple('Segments', '''def_vars  cons  objs  var_bnds 
                                     initial_guess  con_jacobian  eval_order
                                     con_blocks  var_blocks''')

def get_segments(nl_header, content):
    # Returns Segments.
    n_cons = nl_header.n_cons
    n_objs = nl_header.n_objs
    allsegs = Segments(def_vars=OrderedDict(), cons=[None]*n_cons, 
                       objs=[None]*n_objs, var_bnds=[], initial_guess={}, 
                       con_jacobian=[], eval_order=[], con_blocks=[],  
                       var_blocks=[])
    # Creates a sequence of tuples. Two adjacent tuples are the header and the
    # body of the segment:
    #     V68 0 0
    #     o3         -> ('V68 0 0',), ('o3', 'v2', 'v53')
    #     v2
    #     v53
    segments_w_header = list(partitionby(lambda l: l[0] in SEGMENTS, content))
    segments_w_header = fix_empty_segments(segments_w_header)
    # Assumption: Constraints precede r segments. Enforce it:
    segments_w_header = move_r_beyond_last_C(segments_w_header)
    # Populate allsegs by dispatching to the corresponding segment handler. 
    for header, segment in segments_w_header:
        ADD_SEGMENT.get(header[0], ignored)(header, segment, allsegs)
    return allsegs

def move_r_beyond_last_C(segments_w_header):
    r_index = first_index_of(segments_w_header, 'r')
    if r_index is None:
        warning('Apparently the problem does not have r segments.')
        return segments_w_header
    last_C_index = last_index_of(segments_w_header, 'C')
    assert last_C_index is not None
    if r_index > last_C_index:
        return segments_w_header
    pos = last_C_index + 1
    new_lst = segments_w_header[:r_index] # elems before r segment
    new_lst.extend(segments_w_header[r_index+1:pos]) # elems after r till last C
    new_lst.append(segments_w_header[r_index]) # move r segment here
    new_lst.extend(segments_w_header[pos:])    # copy the rest after last C
    assert len(new_lst) == len(segments_w_header)
    return new_lst    

def first_index_of(segments_w_header, header):
    # returns None if not found
    return next((i for i, e in enumerate(segments_w_header) if e[0][0]==header), None)

def last_index_of(segments_w_header, header):
    # returns None if not found
    idx = first_index_of(segments_w_header[::-1], header)
    return len(segments_w_header) - idx - 1 if idx is not None else None

# TODO Write an appropriate algorithm to replace partitionby, which then 
#      eliminates the need for this post-processing.
def fix_empty_segments(segments_w_header):
    # The assumption is that each header is followed by a segment. However,
    # if the segment is empty, then the header is followed by an another header.
    new_segments_with_header = [ ]
    for headers, segment in as_pairs(segments_w_header):
        n_headers = len(headers)
        assert n_headers
        if n_headers == 1:
            # Everything is fine, just one header with a segment
            new_segments_with_header.append((headers[0], segment))
        else:
            # Insert fake empty segments after each header
            segments = [tuple() for _ in irange(n_headers-1)]
            segments.append(segment)
            for h, s in izip(headers, segments):
                new_segments_with_header.append((h, s))
    return new_segments_with_header

def V_segment(header, segment, allsegs):
    index, length, _k = to_ints(header[1:])
    linear_part, nl_part = segment[:length], segment[length:] 
    linear_part = parse_linear_part(linear_part)
    V_segs = allsegs.def_vars
    assert index not in V_segs, index
    V_segs[index] = (nl_part, linear_part)
    allsegs.eval_order.append(('V', index))

def C_segment(header, segment, allsegs):
    i = int(header[1:])
    cons = allsegs.cons
    assert cons[i] is None
    cons[i] = (segment, None, None)
    allsegs.eval_order.append(('C', i))

def J_segment(header, segment, allsegs):
    i, j = to_ints(header[1:])
    indices, coefficients = parse_linear_part(segment)
    assert len(indices) == j, i
    nl_part, lin_part, rng = allsegs.cons[i]
    assert lin_part is None, (i, lin_part)
    allsegs.cons[i] = (nl_part, (indices, coefficients), rng)

def O_segment(header, segment, allsegs):
    # similar to the C segment (nl_part of the constraints)
    i, j = to_ints(header[1:])
    O_segs = allsegs.objs
    assert O_segs[i] is None
    direction = 'minimize' if j == 0 else 'maximize'
    O_segs[i] = (direction, segment, None)
    allsegs.eval_order.append(('O', i))

def G_segment(header, segment, allsegs):
    # similar to the J segment (lin_part of the constraints)
    i, j = to_ints(header[1:])
    indices, coefficients = parse_linear_part(segment)
    assert len(indices) == j, i
    direction, nl_part, lin_part = allsegs.objs[i]
    assert lin_part is None, (i, lin_part)
    allsegs.objs[i] = (direction, nl_part, (indices, coefficients))

def r_segment(_header, segment, allsegs):
    cons = allsegs.cons
    for i, rng in enumerate(imap(to_range, segment)):
        nl_part, lin_part, _None = cons[i]
        cons[i] = (nl_part, lin_part, rng)

def b_segment(_header, segment, allsegs):
    allsegs.var_bnds.extend(imap(to_range, segment))

def x_segment(_header, segment, allsegs):
    allsegs.initial_guess.update(imap(to_index_value, segment))

def k_segment(_header, segment, allsegs):
    allsegs.con_jacobian.extend(imap(int, segment))

def S_segment(header, segment, allsegs):
    kind, _length, name = header.split()
    if name != 'blockid':
        warning('S segment "%s" ignored' % header)
        return
    kind = int(kind[1:])
    assert kind == 0 or kind == 1, (kind, header)
    block_ids = allsegs.con_blocks if kind == 1 else allsegs.var_blocks
    block_ids.extend(parse_blockids(segment))

def ignored(header, *_args, **_kwargs):
    warning('%s segment ignored' % header[0])

def to_index_value(line):
    parts = line.split()
    assert len(parts) == 2
    return int(parts[0]), parts[1]

def parse_linear_part(segment):
    indices, coeffs = [ ], [ ]
    for i, c in imap(lambda l: l.split(), segment):
        indices.append(int(i))  # make the indices integers 
        coeffs.append(c)        # but keep the coefficients as strings
    return indices, coeffs

def parse_blockids(segment):
    return [(i, j) for i, j in imap(to_ints, segment)]

ADD_SEGMENT = {
    'C': C_segment,
    'V': V_segment,
    'J': J_segment,
    'O': O_segment,
    'G': G_segment,
    'r': r_segment,
    'b': b_segment,
    'x': x_segment,
    'k': k_segment,
    'S': S_segment,
}

class Range:
    lb_ub  = 0  # lb <= body <= ub
    ub     = 1  # body <= ub
    lb     = 2  # lb <= body
    no_bnd = 3  # body is unbounded
    eq     = 4  # body = c
    compl  = 5  # "5 k i" in the .nl file, complementary constraint, ignoring it

TO_RANGE_KIND = {
    '0': Range.lb_ub,
    '1': Range.ub,
    '2': Range.lb,
    '3': Range.no_bnd,
    '4': Range.eq,
    '5': Range.compl
}

def to_range(line):
    parts = line.split()    
    return TO_RANGE_KIND[parts[0]], tuple(parts[1:])

#-------------------------------------------------------------------------------
# Interpreting nonlinear constraint segments
# https://en.wikipedia.org/wiki/Reverse_Polish_notation#Postfix_algorithm
#
# While there are input tokens left
#     Read the next token from input.
#     If the token is a value
#         Push it onto the stack.
#     Otherwise, the token is an operator (operator here includes both operators and functions).
#         It is already known that the operator takes n arguments.
#         If there are fewer than n values on the stack
#             (Error) The user has not input sufficient values in the expression.
#         Else, Pop the top n values from the stack.
#         Evaluate the operator, with the values as arguments.
#         Push the returned results, if any, back onto the stack.
# If there is only one value in the stack
#     That value is the result of the calculation.
# Otherwise, there are more values in the stack
#     (Error) The user input has too many values.

def interpret(body, counter):
    segment = list(body)
    stack = []
    evaluated = ('n', 'v', 't')
    while segment:
        token = segment.pop()
        kind = token[0]
        if kind in evaluated or token.isdigit():
            stack.append(token)
        else:
            assert kind == 'o', kind
            assert token in NAME, token
            arity = ARITY[token]
            if not arity: # n-ary operators have their arity as the first argument
                arity = int(stack.pop())
            assert len(stack) >= arity
            args = [stack.pop() for _ in irange(arity)]
            tmp = 't%d' % next(counter)
            #print(tmp, '=', NAME[token], ' '.join(args))
            yield tmp, NAME[token], args
            stack.append(tmp)
    assert len(stack) == 1, stack
    assert stack[0][0] in evaluated, stack
    yield stack[0]

#-------------------------------------------------------------------------------

def get_problem(name):
    nl_header, segments = parse_nl_file(name)   
    col_names     = get_col_names(name)
    def_var_names = get_def_var_names(name)
    row_names     = get_row_names(name)
    solutions     = get_solutions(name, col_names)
    return Problem(name=name, nl_header=nl_header, segments=segments, 
                   col_names=col_names, def_var_names=def_var_names, 
                   row_names=row_names, solutions=solutions)

def get_J_rowwise(problem):
    cons = problem.segments.cons
    # Get the sparsity pattern as a sparse matrix rowwise
    return [indices for _nl_part, (indices, _coefficients), _rng in cons]
