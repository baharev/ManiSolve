# Copyright (C) 2017, 2018 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from collections import namedtuple
from functools import partial
from glob import glob
from os import mkdir, remove, makedirs
from os.path import isdir, isfile
import numpy as np
from cffi_solve import CProblem
from c_sub_check import create_dag, generate_c_code, compile_c_code, \
                        TMP_DIR, clean_up_intermediate_files,  get_so_name, \
                        gen_x_r_slices
from c_subprobs import var_idx_order, con_idx_order, get_var_bnds
from perturb import perturb_C
from subsampling2 import subsample2
from utils import print_timing, warning


TestCase = namedtuple('TestCase', 'interesting  h_max  so_suffix')

PROBLEMS = {
    'blockEx': TestCase(interesting=['x[%d]' % i for i in range(1, 21)],
                         h_max=5, so_suffix='_plot'),
            
    'spider2D': TestCase(interesting=['x[%d]' % i for i in range(1, 41)],
                         h_max=5, so_suffix='_plot'),
    'mss10_4': 
        TestCase(interesting=['x[3,%d]' % i for i in range(10, 0, -1)],
                 #interesting=['T[%d]' % i for i in range(10, 0, -1) ],
                 h_max=4, so_suffix='_plot'),
    
    'mss20_4': 
        TestCase(interesting=['x[3,%d]' % i for i in range(20, 0, -1)],
                 #interesting=['x[3,%d]' % i for i in range(20, 0, -1)],
                 #interesting=['T[%d]' % i for i in range(20, 0, -1) ],
                 h_max=4, so_suffix='_plot'),
    
    'mss40_4': 
        TestCase(interesting=['x[3,%d]' % i for i in range(40, 0, -1)],
                 #interesting=['T[%d]' % i for i in range(40, 0, -1) ],
                 h_max=4, so_suffix='_plot'),
            
    'mss60_4': 
        TestCase(interesting=['x[3,%d]' % i for i in range(60, 0, -1)],
                 #interesting=['T[%d]' % i for i in range(60, 0, -1) ],
                 h_max=8, so_suffix='_plot'),
            
    'mss60_A': 
        TestCase(interesting=['x[3,%d]' % i for i in range(60, 0, -1)],
                 #interesting=['T[%d]' % i for i in range(60, 0, -1) ],
                 h_max=8, so_suffix='_plot'),

    'mss60_B': 
        TestCase(interesting=['x[3,%d]' % i for i in range(60, 0, -1)],
                 #interesting=['T[%d]' % i for i in range(60, 0, -1) ],
                 h_max=6, so_suffix='_plot'),

    'mss75_B': 
        TestCase(interesting=['x[3,%d]' % i for i in range(75, 0, -1)],
                 #interesting=['T[%d]' % i for i in range(75, 0, -1) ],
                 h_max=8, so_suffix='_plot'),

    'extrSeq': 
        TestCase(interesting=['x[3,%d]' % i for i in range(30, 0, -1)],
                 #interesting=['T[%d]' % i for i in range(60, 0, -1) ],
                 h_max=8, so_suffix='_plot'),
}

#-------------------------------------------------------------------------------
# WARNING: We assume an upper envelope with square blocks except the first and 
# the last block which are rectangular (under- and over-determined, resp.)
#-------------------------------------------------------------------------------

_N_PTS = 100

N_INITIAL_POINTS = 400
N_POINTS_BACKSLV = int(10*_N_PTS)
N_POINTS_PERTURB = int(10*_N_PTS)
NEW_PTS_PER_SUBP = 200

K_OVERSAMPLE = 5

K_MATCHES = 20

# Note: Fixing x_new and a *strict* tol -> poor performance ~ rejection sampling
TOL_BACK =   1.0E-6 # When we are inserting new points out of the blue
TOL_PERT =   1.0E-6 # When we apply some artifical constraint violations too
TOL_LAST =   1.0E-6 # At the last h_max subproblems

TOL_TRY    = 1.0e-3 # Try to fix bound violations less then this threshold
TOL_ACCEPT = 1.0e-6 # Acceptable constraint violation after fixing bound violations

# And we are subsampling the point clump!

#-------------------------------------------------------------------------------

_fwd = None
_bwd = None

def main():
    #name = 'spider2D'
    #name = 'blockEx'
    name = 'mss60_B'
    #name = 'extrSeq'
    interesting, h_max, so_suffix = PROBLEMS[name]
    #---
    if not isdir(TMP_DIR):
        print('Creating folder "{}"'.format(TMP_DIR))
        makedirs(TMP_DIR)
    #---
    g, problem = create_dag(name)
    manifold_dim = get_manifold_dim(problem)
    print('Manifold dimension:', manifold_dim)
    global bwd_n_fixed_vars
    bwd_n_fixed_vars = partial(__bwd_n_fixed_vars, manifold_dim)
    #---
    fwd_so_path = TMP_DIR + get_so_name(problem.name, 0, so_suffix)
    bwd_so_path = TMP_DIR + get_so_name(problem.name, h_max, so_suffix)
    if not isfile(fwd_so_path) or not isfile(bwd_so_path):
        print('Generating native code!')
        # Forwardsolve stuff
        generate_c_code(g, problem, 0, fwd_n_fixed_vars, so_suffix=so_suffix)
        compile_c_code(problem.name, 0, so_suffix)
        clean_up_intermediate_files()
        # Backsolve stuff
        generate_c_code(g, problem, h_max, bwd_n_fixed_vars, so_suffix=so_suffix)
        compile_c_code(problem.name, h_max, so_suffix)
        clean_up_intermediate_files()
    else:
        print('Using cached native code')
    #---
    global _fwd, _bwd
    _fwd = CProblem(fwd_so_path)
    _bwd = CProblem(bwd_so_path)
    #---
    cascading_solve(problem, h_max, so_suffix, interesting)

def fwd_n_fixed_vars(n_cons, n_vars):
    # Used both by the C code generation and by the iteration logic
    return max(0, n_vars - n_cons)

bwd_n_fixed_vars = None  # will be set when the manifold dim is available

def __bwd_n_fixed_vars(manifold_dim, n_cons, n_vars):
    dof = n_vars - n_cons
    if dof > 0:
        return dof
    elif dof == 0:
        return manifold_dim
    else:
        return 0

def get_manifold_dim(problem):
    x_slc, r_slc = next(gen_x_r_slices(problem, 0))
    n_cons, n_vars = length(r_slc.subp), length(x_slc.subp)
    assert n_cons == 0 and n_vars > 0, (n_cons, n_vars)
    return n_vars

#-------------------------------------------------------------------------------

@print_timing
def cascading_solve(problem, h_max, so_suffix, interesting):
    # More or less copies solve_check_fixed from c_sub_check.py
    meta = get_ordered_info(problem)
    #---------------------------------------------------------------
    # var names, bounds, solutions, and indices in x are dumped
    setup_data_for_plotter(interesting, meta)
    #---------------------------------------------------------------
    np.random.seed(1)
    x_2D, r_2D = solve_setup(problem, N_INITIAL_POINTS)
    fwd_slices = [slcs for slcs in gen_x_r_slices(problem, 0)]
    bwd_slices = [slcs for slcs in gen_x_r_slices(problem, h_max)]
    n_slices = len(fwd_slices)
    assert n_slices == len(bwd_slices), (n_slices, len(bwd_slices))
    assert h_max >= 1, (h_max, 'backsolves require this')
    last = n_slices - 1
    lb, ub = meta.lb, meta.ub
    # Initialize the torn variables
    x_slc, r_slc = fwd_slices[0]
    n_cons, n_vars = length(r_slc.subp), length(x_slc.subp)
    assert n_cons == 0 and n_vars > 0, (n_cons, n_vars)
    x_2D[:,x_slc.subp] = random_sample(lb[x_slc.subp], ub[x_slc.subp], len(x_2D))
    dump_data('Initialization', 0, x_2D[:,x_slc.seen])
    #
    for index in range(1, last):
        old_points = true_mask(len(x_2D)) # <-- The index of the last old point  
        n_pts = len(x_2D)                 #     would have been sufficient
        # Keep the next line here: crash due to ASAN -> we know where it happened
        print('Forward')
        detailed_info(problem, index, len(x_2D), meta, fwd_slices, h_max, so_suffix)
        print('Backward')
        detailed_info(problem, index, len(x_2D), meta, bwd_slices, h_max, so_suffix)
        #dbg_x_2D = x_2D.copy() # for plotting those points where VA27 failed
        #
        # Do the forward solve
        x_slc, r_slc = fwd_slices[index]
        n_cons, n_vars = length(r_slc.subp), length(x_slc.subp)
        assert n_cons > 0 and n_vars > 0, (index, n_cons, n_vars)
        _fwd.solve(index, x_2D, r_2D, iprint=0)
        #
        # Discard the failed ones, but keep the ones with bound violations
        x_2D, r_2D, mask = discard_failed_ones(x_2D, r_2D, x_slc, r_slc)
        old_points = old_points[mask]
        log_losses('solve failed', n_pts, len(x_2D))
        #dump_data('Fwd solve failed', index, dbg_x_2D[~mask,x_slc.seen])
        dump_data('After fwd solve with bound violations', index, x_2D[:,x_slc.seen])
        #
        # Do the backsolve which inserts new points
        x_slc, r_slc = bwd_slices[index]
        x_2D, r_2D, old_points = backsolve(index, x_2D, r_2D, old_points, h_max, 
                                                               bwd_slices, meta)
        # It is important that we discard bound violations only after backsolve
        repair_bnds(index, x_2D, r_2D, lb, ub, x_slc.subp, r_slc.subp)
        n_pts = len(x_2D)
        x_2D, r_2D, mask = discard_failed_ones(x_2D, r_2D, x_slc, r_slc)
        old_points = old_points[mask]
        log_losses('solve failed and bound infeas', n_pts, len(x_2D))
        dump_data('No bound violations', index, x_2D[:,x_slc.seen])
        #
        # Discard those newly inserted points that are too close in x_slc.new,
        # but keep all old_points
        selected = subsample_new_points(x_2D, x_slc.new, old_points)
        x_2D, r_2D = x_2D[selected], r_2D[selected]
        old_points = old_points[selected]
        dump_data('After subsampling backsolve', index, x_2D[:,x_slc.seen])
        #
        # Diagnostics
        assert x_2D.shape == r_2D.shape
        assert len(old_points) == len(x_2D), 'old_points not updated properly'
        # FIXME Enable dbg_violations when done with parameter tuning
        #dbg_violations(index, x_2D, r_2D, x_slc, r_slc, lb, ub)
        assert np.isfinite(x_2D[:,x_slc.seen]).all(), index
        short_info(problem, index, n_cons, n_vars, h_max, so_suffix)
    #
    # Solve the final, overdetermined block
    index = last
    x_slc, r_slc = bwd_slices[index]
    print('Backward')
    detailed_info(problem, index, len(x_2D), meta, bwd_slices, h_max, so_suffix)
    x_2D, r_2D = final_block(index, x_2D, r_2D, lb, ub, bwd_slices)
    # Also discarded the failed and bound infeasible ones.
    #
    # Diagnostics
    assert x_2D.shape == r_2D.shape
    dbg_violations(index, x_2D, r_2D, x_slc, r_slc, lb, ub)
    assert np.isfinite(x_2D[:,x_slc.seen]).all(), index
    short_info(problem, index, n_cons, n_vars, h_max, so_suffix)

#===============================================================================

def backsolve(index, x_2D, r_2D, old_points, h_max, bwd_slices, meta):
    x_slc, r_slc = bwd_slices[index]
    if index <= h_max:
        # Insert new points out of the blue        
        new_x_2D, new_r_2D = backsolve_initial(index, h_max, x_slc, r_slc, meta)
        msg = 'After backsolve'
    else:
        new_x_2D, new_r_2D = perturb(index, x_2D, r_2D, h_max, bwd_slices, meta)
        msg = 'After perturbation'
    # Both branches discarded failed and bound infeasible points.
    # WARNING: new_r_2D has NaN for those values that are not in subp!
    x_2D, r_2D = np.vstack((x_2D, new_x_2D)), np.vstack((r_2D, new_r_2D))
    old_points = np.hstack((old_points, false_mask(len(new_x_2D))))
    dump_data(msg, index, x_2D[:,x_slc.seen])
    return x_2D, r_2D, old_points 

def final_block(index, x_2D, r_2D, lb, ub, bwd_slices):
    # We skip the forward solve, and we don't subsample.
    # Why only just clip? No good reasons: I would have to change the C codegen.
    assert index == len(bwd_slices)-1, (index, len(bwd_slices))
    x_slc, r_slc = bwd_slices[index]
    # Assumed to be overdetermined, with no new variables:
    assert length(r_slc.subp) == length(x_slc.subp) + length(r_slc.new)
    assert length(x_slc.new) == 0
    #dbg_x_2D = x_2D.copy() # for plotting those points where VA27 failed
    
    print('Solving the last h_max subproblems simultaneously')
    # Solve would work too but it is slower and fails more often.
    #_bwd.solve(index, x_2D, r_2D, tol=TOL_LAST, iprint=0)
    _bwd.solve_from(index, x_2D, r_2D, tol=TOL_LAST, iprint=0)

    n_pts = len(x_2D)
    x_2D, r_2D, _mask = discard_failed_ones(x_2D, r_2D, x_slc, r_slc)
    log_losses('solve failed', n_pts, len(x_2D))
    #dump_data('Last bwd solve failed', index, dbg_x_2D[~mask,x_slc.seen])
    dump_data('After bwd solve with bound violations', index, x_2D[:,x_slc.seen])

    repair_bnds(index, x_2D, r_2D, lb, ub, x_slc.subp, r_slc.subp, just_clip=True)
    n_pts = len(x_2D)
    x_2D, r_2D, _mask = discard_failed_ones(x_2D, r_2D, x_slc, r_slc)
    log_losses('solve failed and bound infeas', n_pts, len(x_2D))
    dump_data('No bound violations', index, x_2D[:,x_slc.seen])
    return x_2D, r_2D

def subsample_new_points(x_2D, x_slc_new, old_points):
    # Handle the edge cases
    if old_points.all() and old_points.any():
        warning('*** Failed to insert any new point ***')
        selected = true_mask(len(old_points))
    elif old_points.any():
        spidx = find_splitpoint(old_points)
        new_selected = subsample2(x_2D[spidx:,x_slc_new], NEW_PTS_PER_SUBP)
        assert len(x_2D) == len(old_points), (len(x_2D), len(old_points))
        selected = true_mask(len(old_points))
        selected[spidx:] = new_selected
    else:
        warning('All old points were lost!')
        selected = subsample2(x_2D[:,x_slc_new], NEW_PTS_PER_SUBP)
    return selected

#===============================================================================

def find_splitpoint(old_points):
    # FIXME We would only need to track the index of the last old point. I did 
    # not do this refactoring. This function assumes that there are new points, 
    # and that the old_points mask is [True, ..., True, False, ... False].
    arr = np.select([old_points], [1,])
    diffs = np.ediff1d(arr, to_begin=0)   
    run_ends, = np.where(diffs != 0)
    old, new = np.split(arr, run_ends)
    assert len(old) > 0
    assert old.all()
    assert len(new) > 0
    assert not new.any()
    return len(old)

#-------------------------------------------------------------------------------

def length(slc):
    return slc.stop - slc.start

def true_mask(shape):
    return np.full(shape, np.True_, np.bool_)

def false_mask(shape):
    return np.full(shape, np.False_, np.bool_)

def log_losses(msg, n_pts, curr_pts):
    if n_pts == 0:
        assert curr_pts == 0
        print('%s: (all points have been lost already)' % msg)
        return
    lost = n_pts - curr_pts
    print('%s: %d/%d' % (msg, lost, n_pts), '(%.1f%%)' % (100.0*lost/n_pts))

def detailed_info(problem, index, n_pts, meta, x_r_slc, h_max, so_suffix):
    x_slc, r_slc = x_r_slc[index]
    n_cons, n_vars = length(r_slc.subp), length(x_slc.subp)
    fmt = '{}, index: {}, size: {}x{}, h_max: {}, so_suffix: "{}"'
    print(fmt.format(problem.name, index, n_cons, n_vars, h_max, so_suffix))
    print('Number of points:', n_pts)
    print('Vars:', ', '.join(meta.var_names[x_slc.subp]))
    print('New: ', ', '.join(meta.var_names[x_slc.new]))
    print('Cons:', ', '.join(meta.con_names[r_slc.subp]))
    print('New: ', ', '.join(meta.con_names[r_slc.new]))
    
def short_info(problem, index, n_cons, n_vars, h_max, so_suffix):
    fmt = '{}, index: {}, fwd size: {}x{}, h_max: {}, so_suffix: "{}"'
    print(fmt.format(problem.name, index, n_cons, n_vars, h_max, so_suffix))
    print()

def dbg_violations(index, x_2D, r_2D, x_slc, r_slc, lb, ub):
    assert len(x_2D), 'Lost all points...'
    for x, r in zip(x_2D, r_2D):
        _bwd.evaluate(index, x, r)
    con_viol = np.linalg.norm(r_2D[:,r_slc.subp].reshape(-1), ord=np.inf)
    lb_, ub_ = lb[x_slc.subp], ub[x_slc.subp]
    x = x_2D[:,x_slc.subp]
    bnd_viol = max(0.0, (lb_-x).max(), (x-ub_).max())
    assert np.isfinite(con_viol) and np.isfinite(bnd_viol), (con_viol, bnd_viol)
    print('Number of points:', len(x_2D))
    print('Max con viol inf:', con_viol)
    print('Max con viol L2: ', np.linalg.norm(r_2D[:,r_slc.subp], axis=1).max())
    print('Max bnd viol:', bnd_viol)

def discard_failed_ones(x_2D, r_2D, x_slc, r_slc):
    x_new = x_2D[:,x_slc.new]
    r_new = r_2D[:,r_slc.new]
    x_mask = np.isfinite(x_new).all(axis=1)
    r_mask = np.isfinite(r_new).all(axis=1)
    mask = x_mask & r_mask 
    return x_2D[mask], r_2D[mask], mask

def solve_setup(problem, n_points):
    n_cons, n_vars = problem.nl_header.n_cons, problem.nl_header.n_vars
    x_2D = np.full((n_points, n_vars), np.nan)  
    r_2D = np.full((n_points, n_cons), np.nan)
    return x_2D, r_2D    

#-------------------------------------------------------------------------------

def nothing(*args, **kwargs):
    pass

log = nothing

def perturb(index, x_2D, r_2D, h_max, bwd_slices, meta):
    assert index > h_max and index < len(bwd_slices) - 1, (index, len(bwd_slices))
    x_slc, r_slc = bwd_slices[index]
    n_cons_subp, n_vars_subp = length(r_slc.subp), length(x_slc.subp)
    assert n_cons_subp == n_vars_subp, (n_cons_subp, n_vars_subp)
    lb, ub = meta.lb, meta.ub
    # Downsample the point clump:
    if len(x_2D) > NEW_PTS_PER_SUBP:
        #x_next, _r_next = bwd_slices[index + 1]
        #slc = slice(x_slc.subp.start, x_next.subp.start) # seen for the last time 
        selected = subsample2(x_2D[:,x_slc.subp], NEW_PTS_PER_SUBP)
        x_2D, r_2D = x_2D[selected], r_2D[selected]  
    #
    A_Ainv_J33 = get_pinv(index, x_2D, r_2D, x_slc, r_slc)
    x_2D_pert, index_pert = get_linear_perturbed_pts(index, x_2D, A_Ainv_J33, 
                                                     bwd_slices, lb, ub)
    #
    # Re-evaluate at the perturbed points, and throw away the failed ones
    #---
    # Make bound feasible first
    #x_slc_subp = x_slc.subp
    for i in range(len(x_2D_pert)):
        # Small random perturbations can improve the resolution below the feedstage
        # if the tolerances are too permissive; otherwise it seems to make matters worse
        #x_2D_pert[i,x_slc_subp] += 0.01*np.random.randn(*(x_2D_pert[i,x_slc_subp].shape))
        x_2D_pert[i] = np.clip(x_2D_pert[i], meta.lb, meta.ub)
    #---
    assert x_2D_pert.shape[1] == x_2D.shape[1]
    r_2D_pert = np.full((len(x_2D_pert), r_2D.shape[1]), np.nan)
    for x, r in zip(x_2D_pert, r_2D_pert):
        _bwd.evaluate(index, x, r)
    #print('Inf norm of pertubed point:', np.linalg.norm(r[r_slc.subp], np.inf))
    r_subp = r_2D_pert[:,r_slc.subp]
    mask = np.isfinite(r_subp).all(axis=1)
    x_2D_pert, r_2D_pert, index_pert = x_2D_pert[mask], r_2D_pert[mask], index_pert[mask]
    #
    #print('<<<')
    #print('before backsolve, perturbed points')
    #dbg_violations(index, x_2D_pert, r_2D_pert, x_slc, r_slc, lb, ub)
    x_2D_pert, r_2D_pert = backsolve_middle(index, x_2D_pert, r_2D_pert, index_pert,
                                            h_max, x_slc, r_slc, meta)
    #print('after backsolve, perturbed points')
    #if len(x_2D_pert) > 0:
    #    dbg_violations(index, x_2D_pert, r_2D_pert, x_slc, r_slc, lb, ub)
    #print('>>>')
    dump_data('Perturbed', index, x_2D_pert[:,x_slc.seen])
    return x_2D_pert, r_2D_pert

def get_pinv(index, x_2D, r_2D, x_slc, r_slc):
    n_cons_subp, n_vars_subp = length(r_slc.subp), length(x_slc.subp)
    x3 = length(x_slc.new)
    assert length(x_slc.new)==length(r_slc.new), (x_slc, r_slc)
    A_Ainv_J33 = []
    for x, r in zip(x_2D, r_2D):
        jac = np.full((n_cons_subp, n_vars_subp), np.nan)
        _bwd.jacobian_evaluation(index, x, r, jac)
        # Assumes that J33 is square
        A_Ainv_J33.append((jac[:,:-x3], np.linalg.pinv(jac[:,:-x3], rcond=1.0E-4), jac[-x3:,-x3:]))
    return A_Ainv_J33

def random_idx_mask(n_points, n_fixed, n_vars):
    # idx in current subproblem slice, starting from 0
    # Admittedly inefficient implementation
    assert n_fixed <= n_vars, (n_fixed, n_vars)
    assert n_fixed > 0, n_fixed
    indices = np.arange(n_vars)
    idx  = np.full((n_points, n_fixed), -1, dtype=np.intc)
    # With the mask, we want to zero out the NOT-selected elements in delta x
    mask = np.full((n_points, n_vars), np.True_, dtype=bool)
    for i in range(n_points):
        select = np.sort(np.random.choice(indices, size=n_fixed, replace=False))
        idx[i] = select
        mask[i, select] = np.False_ 
    return idx, mask

def get_linear_perturbed_pts(index, x_2D, A_Ainv_J33, bwd_slices, lb, ub):
    print('Starting linear perturbations')
    #
    x_slc, r_slc = bwd_slices[index]
    #
    n_fixed = bwd_n_fixed_vars(length(r_slc.subp), length(x_slc.subp))
    # indices in x_new, and not in x_subp, we will have to fix it later!
    n_new = length(x_slc.new)
    indices, idx_mask = random_idx_mask(N_POINTS_PERTURB, n_fixed, n_new)  
    values = random_sample(lb[x_slc.new], ub[x_slc.new], N_POINTS_PERTURB)
    #b = np.full(length(x_slc.subp), 0.0)
    
    #np.set_printoptions(formatter={'float': lambda x: '%.4f' % x}, linewidth=1000)
    
    dr_norm_2D = np.full((len(x_2D), len(values)), np.nan)
    x_pert_2D  = np.full((len(x_2D), len(values), length(x_slc.subp)), np.nan)
    for x, x_pert, dr_norm, (A, Ainv, J33) in zip(x_2D, x_pert_2D, dr_norm_2D, A_Ainv_J33):
        dx_new = values - x[x_slc.new] # Sign here: A*x=b; later: setting dx_full
        dx_new[idx_mask] = 0.0
        perturb_C(dx_new, J33, Ainv, A, x[x_slc.subp], dr_norm, x_pert)
#         # perturb_C does this loop (up until the if statement) but in C:
#         for k, dx3 in enumerate(dx_new):
#             b[-n_new:] = J33 @ dx3
#             dx1_dx2 = Ainv @ b
#             dr = A @ dx1_dx2 - b
#             dr_norm[k] = np.dot(dr, dr)
#             x_pert[k] = x[x_slc.subp] + np.concatenate((-dx1_dx2, dx3))
#             #--------------------------------------------
#             # Only for debugging:
#             if dr_norm[k] < TOL_PERT*length(r_slc.subp):
#                 cnt += 1
#                 y = x.copy()
#                 y[x_slc.subp] = x_pert[k]
#                 r = np.full(x.shape, np.nan, dtype=np.double)
#                 _bwd.evaluate(index, y, r)
#                 print('y:', y[x_slc.seen])
#                 print('r:', r[r_slc.subp])

#    print('Candidates:', cnt)

    print('Candidates computed')
    
    assert np.isfinite(x_pert_2D).all()
    
    #print('dr_norm', dr_norm_2D.T.shape)
    #print('x_pert', np.swapaxes(x_pert_2D, 0, 1).shape)
    pert_x, pert_idx = [], []
    
    offset = x_slc.new.start - x_slc.subp.start
    for k, (dr_norm, x_pert) in enumerate(zip(dr_norm_2D.T, np.swapaxes(x_pert_2D, 0, 1))):
        sorter = np.argsort(dr_norm, kind='mergesort')        
        cutoff = np.searchsorted(dr_norm, TOL_PERT*length(r_slc.subp), sorter=sorter)
        cutoff = max(1, cutoff) # Always add the best point, even if not promising
        small_r = sorter[:cutoff]
        x_perturbed = x_2D[small_r]
        x_perturbed[:,x_slc.subp] = x_pert[small_r]
        idx_perturbed = np.tile(indices[k]+offset, (len(small_r), 1)) # index_pert in .subp, hence the offset
        if cutoff > K_MATCHES:
            # Do the downsampling here
            mask = subsample2(x_perturbed[:,x_slc.new], K_MATCHES)
            x_perturbed   = x_perturbed[mask]
            idx_perturbed = idx_perturbed[mask]
        pert_x.extend(x_perturbed)
        pert_idx.extend(idx_perturbed)
    
    return np.array(pert_x), np.array(pert_idx)

def random_indices(n_points, n_fixed, indices):
    # idx in current subproblem slice, starting from 0
    # Admittedly inefficient implementation
    assert indices.ndim == 1
    assert n_fixed <= len(indices), (n_fixed, len(indices))
    assert n_fixed > 0, n_fixed
    idx = np.full((n_points, n_fixed), -1, dtype=np.intc)
    for i in range(n_points):
        idx[i] = np.sort(np.random.choice(indices, size=n_fixed, replace=False))
    return idx

def idx_val_uniform(n_points, n_fixed, indices, lb_subp, ub_subp):
    # idx in current subproblem slice, starting from 0
    # idx must be valid indices in lb and ub (the current subprolem slice)
    idx = random_indices(n_points, n_fixed, indices)
    lb, ub = np.take(lb_subp, idx), np.take(ub_subp, idx)
    val = np.random.uniform(lb, ub)
    assert idx.shape == val.shape, (idx.shape, val.shape)
    return idx, val

#-------------------------------------------------------------------------------
# In the next 3 functions the first assert tells when the function is called,
# and the second assert checks our assumptions regarding the sparsity pattern

def backsolve_initial(index, h_max, x_slc, r_slc, meta):
    assert index >= 1 and index <= h_max, index
    # Assumed to be underdetermined
    assert length(r_slc.subp) < length(x_slc.subp)
    log_on_enter(index, x_slc, r_slc, meta)
    # We invent new points out of the blue
    x_2D = np.full((N_POINTS_BACKSLV, len(meta.var_names)), np.nan)  
    r_2D = np.full((N_POINTS_BACKSLV, len(meta.con_names)), np.nan) 
    # The subsampling in random_sample scales poorly if n_points > 5000:
    lb, ub = meta.lb, meta.ub
    idx, val =  fix_x_new_uniformly(lb, ub, len(x_2D), x_slc, r_slc)
    #---
#     # FIXME Hack!
#     name = meta.problem_name
#     if (name == 'mss60_4') or (name.startswith('mss60_') and index != 1):
#         x1_i, x3_i = None, None
#         for i in range(x_slc.new.start, x_slc.new.stop):
#             name = meta.var_names[i]
#             if name.startswith('x[1,'):
#                 print('x1:', name)
#                 print('idx:', i - x_slc.subp.start)
#                 x1_i = i - x_slc.subp.start
#             if name.startswith('x[3,'):
#                 print('x3:', name)
#                 print('idx:', i - x_slc.subp.start)
#                 x3_i = i - x_slc.subp.start
#         # idx = i - x_slc.subp.start
#         n_hacked = 50 + 6
#         idx[-n_hacked:,0] = x1_i - x_slc.subp.start
#         val[-n_hacked:,0] = 0.0
#         idx[-n_hacked:,1] = x3_i - x_slc.subp.start
#         val[-n_hacked:-6,1] = np.linspace(lb[x3_i], ub[x3_i], num=n_hacked-6)
#         eps = 0.005
#         val[-6, 0] = ub[x1_i] 
#         val[-6, 1] = lb[x3_i]
#         val[-5, 0] = ub[x1_i] - eps
#         val[-5, 1] = lb[x3_i]
#         val[-4, 0] = ub[x1_i] - 2*eps
#         val[-4, 1] = lb[x3_i] 
#         val[-3, 0] = ub[x1_i] - eps
#         val[-3, 1] = lb[x3_i] + eps
#         val[-2, 0] = ub[x1_i] - 2*eps
#         val[-2, 1] = lb[x3_i] + 2*eps
#         val[-1, 0] = ub[x1_i] - 2*eps
#         val[-1, 1] = lb[x3_i] + eps
    #---
    _bwd.solve_fixed(index, x_2D, r_2D, idx, val, tol=TOL_BACK)
    #---
#     for x, r in zip(x_2D, r_2D):
#         if not np.isfinite(r[r_slc.seen]).all() or np.absolute(r[r_slc.seen]).max() > TOL_BACK:
#             continue
#         show_hline = True
#         for i in range(x_slc.seen.start, x_slc.seen.stop):
#             name, lo, up = meta.var_names[i], lb[i], ub[i]
#             if x[i] < lo - 1.0e-4 or x[i] > up + 1.0e-4:
#                 if show_hline:
#                     show_hline = False 
#                     print('---------------------------------------------------')
#                 if x[i] < lo - 1.0e-4:
#                     print('%.3f' % x[i], '<', '%.3f' % lo, ' %s' % name)
#                 if x[i] > up + 1.0e-4:
#                     print('%.3f' % x[i], '>', '%.3f' % up, ' %s' % name)
    #---
    return get_good_points(index, x_2D, r_2D, x_slc, r_slc, lb, ub)

def backsolve_middle(index, x_2D, r_2D, index_pert, h_max, x_slc, r_slc, meta):
    assert index > h_max, index
    # Assumed to be square before fixing vars:
    assert length(r_slc.subp) == length(x_slc.subp) # and we fix in addition vars
    log_on_enter(index, x_slc, r_slc, meta)
    idx, val =  fix_x_new_to_their_current_value(x_2D, x_slc.subp, index_pert)
    r_2D[:] = np.nan # We do not recognize failed ones otherwise (it was evaluated)
    _bwd.solve_fixed_from(index, x_2D, r_2D, idx, val, tol=TOL_PERT, iprint=0)
    return get_good_points(index, x_2D, r_2D, x_slc, r_slc, meta.lb, meta.ub)

#-------------------------------------------------------------------------------

def repair_bnds(index, x_2D, r_2D, lb, ub, x_slc_subp, r_slc_subp, just_clip=False):
    # x_slc and r_slc *must* be a backward slice since we call backsolve 
    lo, up = lb[x_slc_subp], ub[x_slc_subp]
    lb_viol = np.maximum(lo - x_2D[:,x_slc_subp], 0.0)
    ub_viol = np.maximum(x_2D[:,x_slc_subp] - up, 0.0)
    lb_norm = np.linalg.norm(lb_viol, axis=1)
    ub_norm = np.linalg.norm(ub_viol, axis=1)
    err_norm = np.maximum(lb_norm, ub_norm) 
    error = err_norm**2 / length(x_slc_subp)
    # Candidate: finite (not NaN), violated, and it is less then tol
    r_finite = np.isfinite(r_2D[:,r_slc_subp]).all(axis=1)
    x_finite = np.isfinite(x_2D[:,x_slc_subp]).all(axis=1)
    solver_ok = x_finite & r_finite 
    violated  = error > 0.0
    small_viol= error < TOL_TRY
    candidate = solver_ok & violated &  small_viol
    too_bad   = solver_ok & violated & ~small_viol
    x_2D[too_bad, x_slc_subp] = np.nan
    r_2D[too_bad, r_slc_subp] = np.nan
    print('Bound violation too large:', too_bad.sum())
    print('Trying to repair', candidate.sum(), 'points')
    # Try the dumb clipping first
    (indices,) = np.where(candidate)
    project_back_to_box(index, x_2D, r_2D, x_slc_subp, r_slc_subp, indices, lo, up)
    clipping_fixed = np.isfinite(r_2D[candidate, r_slc_subp]).all(axis=1)
    print('Clipping fixed:', clipping_fixed.sum())
    if just_clip:
        return
    
    # We try again, but now with the local solver. We fix the |J| most violated
    # variables; if we have less, we pick the remaining ones at random.
    try_again = np.compress(~np.isfinite(r_2D[indices,r_slc_subp]).all(axis=1), indices)
    card_J = bwd_n_fixed_vars(length(r_slc_subp), length(x_slc_subp))
    assert card_J > 0, (card_J, index) # Must call with a backward slice!
    bnd_viol = np.maximum(lb_viol, ub_viol)
    for i in try_again:
        x, r, viol = x_2D[i], r_2D[i], bnd_viol[i]
        idx = np.argsort(viol, kind='quicksort')[-card_J:] # quicksort makes a random choice for us
        idx = np.sort(idx)
        val = x[x_slc_subp][idx]
        x.shape = (1, x.shape[0])
        r.shape = (1, r.shape[0])
        idx = np.array(idx, dtype=np.intc)
        idx.shape = (1, idx.shape[0])
        val.shape = (1, val.shape[0])
        _bwd.solve_fixed_from(index, x, r, idx, val, tol=TOL_ACCEPT, iprint=0)
    # We still have to clip them again!
    project_back_to_box(index, x_2D, r_2D, x_slc_subp, r_slc_subp, try_again, lo, up)
    solver_repaired = np.isfinite(r_2D[try_again, r_slc_subp]).all(axis=1)
    print('Solver repaired:', solver_repaired.sum())
    # assert that all succeeded points are bound feasible 
    succeeded = np.isfinite(r_2D[:,r_slc_subp]).all(axis=1)
    x_2D_good = x_2D[succeeded,x_slc_subp]
    bound_feas = (lo <= x_2D_good).all(axis=1) & (x_2D_good <= up).all(axis=1)
    assert bound_feas.all()

def project_back_to_box(index, x_2D, r_2D, x_slc_subp, r_slc_subp, indices, lo, up):
    x_2D[indices,x_slc_subp] = np.clip(x_2D[indices,x_slc_subp], lo, up)
    for i in indices:
        x, r = x_2D[i], r_2D[i]
        _bwd.evaluate(index, x, r)
        resid = r[r_slc_subp]
        error = np.dot(resid, resid) / len(resid)
        if (not np.isfinite(resid).all()) or (error >= TOL_ACCEPT):
            r[r_slc_subp] = np.nan

#-------------------------------------------------------------------------------

# FIXME Grep for x_new and fix: we are fixing only a random subset of it

def _x_new_indices(x_slc):
    offset = x_slc.subp.start
    start, stop = x_slc.new.start - offset, x_slc.new.stop - offset
    assert start > 0, start
    assert stop <= length(x_slc.subp), (stop, length(x_slc.subp)) 
    indices = np.arange(start, stop)
    return indices

def fix_x_new_uniformly(lb, ub, n_points, x_slc, r_slc):
    # backsolve initial
    n_fixed = bwd_n_fixed_vars(length(r_slc.subp), length(x_slc.subp))
    indices = _x_new_indices(x_slc)
    return idx_val_uniform(n_points, n_fixed, indices, lb[x_slc.subp], ub[x_slc.subp])

def random_sample(lb, ub, n_points):
    x = np.random.uniform(lb, ub, (K_OVERSAMPLE*n_points, len(lb)))
    mask = subsample2(x, n_points)
    return x[mask]

def fix_x_new_to_their_current_value(x_2D, x_slc_subp, idx):
    # backsolve middle
    val = np.full(idx.shape, np.nan)
    for i, (x_val, indices) in enumerate(zip(x_2D[:, x_slc_subp], idx)):
        val[i] = x_val[indices]
    return idx, val

def log_on_enter(index, x_slc, r_slc, meta):
    global log  # don't forget: log = nothing in get_good_points()
    log = print
    log()
    log('### In backsolve ###')
    log('index:', index)
    log('size: {}x{}'.format(length(r_slc.subp), length(x_slc.subp)))
    log('fixed:', bwd_n_fixed_vars(length(r_slc.subp), length(x_slc.subp)))
    log('bwd vars:', ', '.join(meta.var_names[x_slc.subp]))
    log('bwd cons:', ', '.join(meta.con_names[r_slc.subp]))
    log('x_new:   ', ', '.join(meta.var_names[x_slc.new]))
    log()

def get_good_points(index, x_2D_orig, r_2D_orig, x_slc, r_slc, lb, ub):
    global log
    n_points = len(x_2D_orig)
    repair_bnds(index, x_2D_orig, r_2D_orig, lb, ub, x_slc.subp, r_slc.subp)
    x_2D, r_2D, _mask = discard_failed_ones(x_2D_orig, r_2D_orig, x_slc, r_slc)
    log_losses('solve failed and bound infeas', n_points, len(x_2D))
    log('### End of backsolve ###')
    log()
    log = nothing
    return x_2D, r_2D

#-------------------------------------------------------------------------------

PLOTTER_DIR = '/tmp/plotter/'

def clean_plotter_dir():
    if isdir(PLOTTER_DIR):
        for f in glob(PLOTTER_DIR + '*'):
            remove(f)
    else:
        mkdir(PLOTTER_DIR)

def setup_data_for_plotter(interesting, meta, delete_dir=True):
    if delete_dir:
        clean_plotter_dir()
    lb, ub, sol_2D = meta.lb, meta.ub, meta.sol_2D
    name_to_xindex = meta.name_to_xindex
    indices = np.fromiter((name_to_xindex[n] for n in interesting), np.int)
    v_min = np.min(lb[indices])
    v_max = np.max(ub[indices])
    print('x bounds for plotting:', v_min, v_max)
    print('indices to watch in x:', indices)
    dump('name.txt', meta.problem_name)
    dump('varnames.txt',     '\n'.join(interesting))
    dump('bounds.txt',       '%f  %f\n' % (v_min, v_max))
    np.save(PLOTTER_DIR + 'indices_in_x.npy', indices)
    write_name_to_index_map('name_to_index_map.txt', name_to_xindex)
    np.save(PLOTTER_DIR + 'solutions.npy', sol_2D)
    print()

def dump(fname, string):
    with open(PLOTTER_DIR + fname, 'w') as f:
        f.write(string)

def write_name_to_index_map(fname, name_to_xindex):
    lst = list('%s  %d\n' % (k, v)  for k, v in sorted(name_to_xindex.items()))
    with open(PLOTTER_DIR + fname, 'w') as f:
        f.writelines(lst)

def dump_data(msg, index, x_2D_slc):
    file_index = 1000*index + dump_data.counter
    dump_data.counter += 1
    msg += ' (index=%d, n_pts=%d)' % (index, x_2D_slc.shape[0])
    with open(PLOTTER_DIR + 't_%d.txt' % file_index, 'w') as f:
        f.write(msg)
    np.save(PLOTTER_DIR + 'x_%d.npy' % file_index, x_2D_slc)

dump_data.counter = 0

#-------------------------------------------------------------------------------

OrderedInfo = namedtuple('OrderedInfo', '''problem_name  con_names  var_names  
                                           sol_2D  lb  ub  name_to_xindex''')

def get_ordered_info(problem):
    var_order = var_idx_order(problem)
    lbs, ubs = get_var_bnds(problem)
    colnames = problem.col_names
    name_to_xindex = dict(zip(colnames, var_order))   
    n_vars = problem.nl_header.n_vars
    perm = np.array(var_order, np.int)
    perm = invert_permutation(perm)
    # perm: from AMPL order to permuted order
    var_names = names_ordered(colnames, perm)
    lb, ub = bounds_ordered(lbs, ubs, perm)
    sol_2D = solutions_ordered(problem.solutions, n_vars, perm)
    # constraints:
    con_perm = np.array(con_idx_order(problem), np.int)
    con_perm = invert_permutation(con_perm)
    con_names = names_ordered(problem.row_names, con_perm)
    return OrderedInfo(problem_name=problem.name, con_names=con_names, 
                       var_names=var_names, sol_2D=sol_2D, lb=lb, ub=ub, 
                       name_to_xindex=name_to_xindex)

def bounds_ordered(lbs, ubs, perm):
    return to_ndarray(lbs)[perm], to_ndarray(ubs)[perm]

def solutions_ordered(solutions, n_vars, perm):
    sol_2D = np.full((len(solutions), n_vars), np.nan) 
    for i, sol in enumerate(solutions):
        sol_2D[i,:] = to_ndarray(sol)[perm]
    return sol_2D

def names_ordered(list_of_strings, perm):
    return to_str_ndarray(list_of_strings)[perm]

def invert_permutation(p):
    '''The argument p is assumed to be some permutation of 0, 1, ..., len(p)-1. 
    Returns an array s, where s[i] gives the index of i in p.'''
    s = np.empty(p.size, p.dtype)
    s[p] = np.arange(p.size)
    return s

def to_ndarray(list_of_strings):
    return np.fromiter(map(float, list_of_strings), np.double)

def to_str_ndarray(list_of_strings):
    max_len = max(map(len, list_of_strings))
    return np.fromiter(list_of_strings, 'S%d' % max_len).astype('U')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
