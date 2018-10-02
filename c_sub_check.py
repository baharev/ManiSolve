# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from random import Random
from subprocess import check_call
import numpy as np
from c_subprobs import TMP_DIR, fake_var_order, gen_blocks, generate_c_code, \
                       compile_sh_path, get_so_name, gen_x_r_slices
from c_wrapper import *
from dag import create_dag
from py3compat import imap, izip
from check_autograd import get_numeric_jacobian
from utils import clean, delete_files, warning


def main():
    np.set_printoptions(formatter={'float': lambda x: '%.4f' % x}, linewidth=1000)
    clean(TMP_DIR)
    for problem_name, h_max, fix_func, so_suffix in gen_testinput():
        g, problem = create_dag(problem_name)
        J = get_J(problem_name)
        generate_c_code(g, problem, h_max, fix_func, so_suffix)
        compile_c_code(problem_name, h_max, so_suffix)
        check_residuals(problem, h_max, so_suffix)
        cross_check_jacobians(problem, J, h_max, so_suffix)
        slice_check(problem, h_max, so_suffix)
        solve_va27_check(problem, h_max, fix_func, so_suffix)
        clean_up_intermediate_files()

def gen_testinput():
    # Only those problems can be tested that have proper S segments
    #probs = ['mss10TornScaled', ]
    probs = ['tunnelDiodes', 'tunnelDiodesSum', 'eco9', 'blockEx', 'bratu20', 
             'cse', 'cse2', 'cse3', 'cse4', 'JacobsenTorn', 'mssTornDbg',
             'pad_nothing', 'pad_con', 'pad_var', 'pad_both']
    start_fmt = 'Started: {}, h_max: {}, fixed: {}'
    done_fmt  = 'Done: {}, h_max: {}, fixed: {}\n'
    for fixed in False, True:
        fix_func = get_n_fixed_vars if fixed else nothing_fixed
        so_suffix = '_f' if fixed else ''
        for h_max in 0, 1, 2, 100:
            for problem_name in probs:
                print(start_fmt.format(problem_name, h_max, fixed))
                #
                yield problem_name, h_max, fix_func, so_suffix
                #
                print(done_fmt.format(problem_name, h_max, fixed))
                print('-----------------------------------------------------\n')

#-------------------------------------------------------------------------------

def get_J(problem_name):
    jac_good = get_numeric_jacobian(problem_name)
    return {(i, j): ds for i, j, ds in jac_good}

#-------------------------------------------------------------------------------

def setup(problem, h_max, so_suffix):
    #
    open_so(TMP_DIR + get_so_name(problem.name, h_max, so_suffix))
    #
    sol = {'v%d' % i: v for i, v in enumerate(problem.solutions[0])}
    var_order = fake_var_order(problem)  # py_point[var_order[vi]] == vi
    n_cons, n_vars = problem.nl_header.n_cons, problem.nl_header.n_vars
    x = np.array([float('NaN')]*n_vars)    
    r = np.array([float('NaN')]*n_cons)
    return x, r, sol, var_order

def check_residuals(problem, h_max, so_suffix):
    x, r, sol, var_order = setup(problem, h_max, so_suffix)
    for index, _con_ids, var_ids in gen_blocks(problem, h_max):
        copy_sol_slice_to_x(sol, x, var_ids, var_order)
        print('x =', x)
        evaluate(index, x, r)
        print('r =', r, '\n')
    r_max = np.max(np.absolute(r)) # Gives NaN if r has NaN
    assert r_max < 1.e-5, (problem.name, r_max)

def cross_check_jacobians(problem, J, h_max, so_suffix):
    x, r, sol, var_order = setup(problem, h_max, so_suffix)
    for index, con_ids, var_ids in gen_blocks(problem, h_max):
        for v in var_ids:
            x[var_order[v]] = sol[v]
        print('x =', x)
        jac = np.full((len(con_ids), len(var_ids)), np.nan)
        #print(jac)
        jacobian_evaluation(index, x, r, jac)
        # No point in checking the Jacobian if the r is wrong
        print('r =', r)
        if con_ids:
            r_max = np.nanmax(np.absolute(r))
            assert r_max < 1.e-5, (problem.name, r_max)
        print('jac =', jac)
        # Now we can cross-check the Jacobian
        compare_elementwise(jac, J, con_ids, var_ids)
        print()
    assert np.allclose(r, 0.0, atol=1.e-5), (problem.name, r)
    print('{}: Jacobian OK!\n'.format(problem.name))

def compare_elementwise(jac, J, con_ids, var_ids):
    def index_of(n):
        return int(n[1:])
    print(jac.shape)
    for i, cidx in enumerate(imap(index_of, con_ids)):
        for j, vidx in enumerate(imap(index_of, var_ids)):
            assert_close(i, j, cidx, vidx, jac, J)

def assert_close(i, j, cidx, vidx, jac, J):
    dr = jac[i,j]
    if (cidx, vidx) not in J:
        assert dr == 0.0, (i, j, cidx, vidx, dr)
    else:
        ds = J[(cidx, vidx)] 
        delta = abs(ds-dr)
        denom = max(abs(ds), abs(dr), 1)
        #print(i, j, ds, dr)
        # Compare with cross_check in reverse_ad.py
        assert delta / denom < 1.0e-12, (i, j, cidx, vidx, ds, dr, delta) 

#-------------------------------------------------------------------------------

def solve_va27_check(problem, h_max, fix_func, so_suffix):
    solve_func = solve_check if not so_suffix else solve_check_fixed
    solve_func(problem, h_max, fix_func, so_suffix)

def solve_setup(problem, h_max, so_suffix):
    n_cons, n_vars = problem.nl_header.n_cons, problem.nl_header.n_vars
    x, r, sol, var_order = setup(problem, h_max, so_suffix)
    x_2D = x.reshape((1, n_vars))
    r_2D = r.reshape((1, n_cons))
    return x, r, sol, var_order, x_2D, r_2D    

def print_info(problem, index, con_ids, var_ids, h_max, so_suffix):
    fmt = '{}, index: {}, size: {}x{}, h_max: {}, so_suffix: {}'
    print(fmt.format(problem.name, index, len(con_ids), len(var_ids), h_max, so_suffix))

def fix_accumulating_error(index, x, r, sol, var_ids, var_order):
    # Sanity check first: The subproblem's solution x_sol must be sane!
    indices = [var_order[v] for v in var_ids]
    x_sol = x[indices] 
    assert np.isfinite(x_sol).all(), (x_sol, x) # VA27 failed?
    # Do the actual work: fix the accumulating error by setting the values in
    # x to their true values, retrieved from sol. It should not change too much.
    x_before = np.copy(x)
    for i, v in izip(indices, var_ids):
        x[i] = sol[v]
    if not np.isclose(x_before, x, equal_nan=True).all():
        warning('Significant deviation from the solution! Multiplicity or singularity?')
        warning(x_before)
        warning(x)
        evaluate(index, x, r)

#-------------------------------------------------------------------------------

def nothing_fixed(*_args):
    return 0

def solve_check(problem, h_max, _fix_func, so_suffix):
    assert so_suffix == ''
    x, r, sol, var_order, x_2D, r_2D = solve_setup(problem, h_max, so_suffix)
    #
    for index, con_ids, var_ids in gen_blocks(problem, h_max):
        # Keep the next line here: crash due to ASAN -> we know where it happened
        print_info(problem, index, con_ids, var_ids, h_max, so_suffix)
        if len(var_ids) > len(con_ids):
            print('Under-determined subproblem, just setting x and r')
            copy_sol_slice_to_x(sol, x, var_ids, var_order)
            evaluate(index, x, r)
        elif var_ids:
            solve(index, x_2D, r_2D, iprint=100)
        else:
            print('Subproblem has no variables, just evaluating r')
            evaluate(index, x, r)
        print('x =', x)
        print('r =', r)
        print_info(problem, index, con_ids, var_ids, h_max, so_suffix)
        fix_accumulating_error(index, x, r, sol, var_ids, var_order)
        print()
    assert np.allclose(r, 0.0, atol=1.e-5), (problem.name, r)

def copy_sol_slice_to_x(sol, x, var_ids, var_order):
    for v in var_ids:
        x[var_order[v]] = sol[v]

#-------------------------------------------------------------------------------

def get_n_fixed_vars(n_cons, n_vars):
    # over-determined -> do nothing
    if n_vars < n_cons:
        return 0
    # only 1 variable to be fixed -> do nothing
    if n_vars == 1:      
        return 0
    # under-determined -> make it square
    if n_vars > n_cons:
        return n_vars - n_cons
    # square system -> make it over-determined by 1 DoF
    return 1

def solve_check_fixed(problem, h_max, fix_func, so_suffix):
    assert so_suffix == '_f'
    x, r, sol, var_order, x_2D, r_2D = solve_setup(problem, h_max, so_suffix)
    rng = Random(3)
    #
    for index, con_ids, var_ids in gen_blocks(problem, h_max):
        # Keep the next line here: crash due to ASAN -> we know where it happened
        print_info(problem, index, con_ids, var_ids, h_max, so_suffix)
        n_fixed_vars = fix_func(len(con_ids), len(var_ids))
        if con_ids and var_ids:  # The normal case: non-degenerate subproblems
            if n_fixed_vars:
                idx, val = get_fixed_vars(n_fixed_vars, rng, var_ids, sol)
                solve_fixed(index, x_2D, r_2D, idx, val, iprint=100)
            else:
                solve(index, x_2D, r_2D, iprint=100)
        elif var_ids or con_ids:
            print('Subproblem has no variables or constraints')
            copy_sol_slice_to_x(sol, x, var_ids, var_order)
            evaluate(index, x, r)
        else:
            raise AssertionError(str((index, n_fixed_vars, con_ids, var_ids)))
        print('x =', x)
        print('r =', r)
        print_info(problem, index, con_ids, var_ids, h_max, so_suffix)
        fix_accumulating_error(index, x, r, sol, var_ids, var_order)
        print()
    assert np.allclose(r, 0.0, atol=1.e-5), (problem.name, r)

def get_fixed_vars(n_fixed_vars, rng, var_ids, sol):
    idx = np.array([-1]*n_fixed_vars, dtype=np.intc)
    val = np.array([float('NaN')]*n_fixed_vars)
    pos_vi = [(pos, vi) for pos, vi in enumerate(var_ids)]
    v_fixed = rng.sample(pos_vi, n_fixed_vars)
    for i, (pos, vi) in enumerate(v_fixed):
        idx[i] = pos
        val[i] = sol[vi]
    print('fixed idx:', idx)
    print('fixed val:', val)
    idx = idx.reshape((1, n_fixed_vars))
    val = val.reshape((1, n_fixed_vars))
    return idx, val

#-------------------------------------------------------------------------------

def slice_check(problem, h_max, so_suffix):
    # Slices = namedtuple('Slices', 'seen  subp  new')
    x_r_slices = [slcs for slcs in gen_x_r_slices(problem, h_max)]
    x, r, sol, var_order = setup(problem, h_max, so_suffix)
    for index, con_ids, var_ids in gen_blocks(problem, h_max):
        print_info(problem, index, con_ids, var_ids, h_max, so_suffix)
        print()
        x_slc, r_slc = x_r_slices[index]
        check_new_part_is_uninitialized(x, r, x_slc, r_slc)
        copy_sol_slice_to_x(sol, x, var_ids, var_order)
        print('x =', x)
        evaluate(index, x, r)
        print('r =', r)
        check_slices(x, r, x_slc, r_slc)
    assert np.allclose(r, 0.0, atol=1.e-5), (problem.name, r)

def check_slices(x, r, x_slc, r_slc):
    print('Seen:')
    x_seen = x[x_slc.seen]
    r_seen = r[r_slc.seen]
    print('x =', x_seen)
    print('r =', r_seen)
    assert (x_seen == x[np.isfinite(x)]).all()
    assert (r_seen == r[np.isfinite(r)]).all()
    assert (x_seen == x[:len(x_seen)]).all()
    assert (r_seen == r[:len(r_seen)]).all()
    #---------------------------------------------
    print('Subproblem:')
    x_subp = x[x_slc.subp]
    r_subp = r[r_slc.subp]
    print('x =', x_subp)
    print('r =', r_subp)
    assert np.isfinite(x_subp).all()
    assert np.isfinite(r_subp).all()
    assert x_subp.size == 0 or (x_subp == x_seen[-len(x_subp):]).all()
    assert r_subp.size == 0 or (r_subp == r_seen[-len(r_subp):]).all()
    #---------------------------------------------
    print('New:')
    x_new = x[x_slc.new]
    r_new = r[r_slc.new]
    print('x =', x_new)
    print('r =', r_new)
    assert np.isfinite(x_new).all()
    assert np.isfinite(r_new).all()
    assert x_new.size == 0 or (x_new == x_seen[-len(x_new):]).all()
    assert r_new.size == 0 or (r_new == r_seen[-len(r_new):]).all()
    #---------------------------------------------
    print()

def check_new_part_is_uninitialized(x, r, x_slc, r_slc):
    x_new = x[x_slc.new]
    r_new = r[r_slc.new]
    assert (~np.isfinite(x_new)).all()
    assert (~np.isfinite(r_new)).all()

#-------------------------------------------------------------------------------

def compile_c_code(problem_name, h_max, so_suffix):
    sh_name = compile_sh_path(problem_name, h_max, so_suffix)
    check_call(['/bin/bash', sh_name], cwd=TMP_DIR)    

def clean_up_intermediate_files():
    delete_files(TMP_DIR + '*.sh')
    delete_files(TMP_DIR + '*.o')
    delete_files(TMP_DIR + '*.c')

if __name__ == '__main__':
    main()
