import pytest
import sys
import os
import importlib
from argparse import Namespace

current_dir = os.path.dirname(os.path.abspath(__file__))
top_level_dir = os.path.abspath(os.path.join(current_dir, '..'))
if current_dir not in sys.path:
    sys.path.append(current_dir)
if top_level_dir not in sys.path:
    sys.path.append(top_level_dir)
tests_dir = os.path.join(top_level_dir, 'tests/problems')
if tests_dir not in sys.path:
    sys.path.append(tests_dir)

from src.dfndfl.run import run_problem, Problem

PROBLEM_FOLDER = 'problems'
PROBLEM_NAMES = [
    "banex",
    "cb2",
    "cb3var20",
    "cb3var30",
    "cb3var40",
    "cb3var50",
    "colville1",
    "crescent",
    "davidon2",
    "elattar",
    "evd61",
    "gill",
    "goffin",
    "hs78",
    "kowalik",
    "l1hilb",
    "l1hilb20",
    "l1hilb30",
    "l1hilb40",
    "lukexp",
    "lukfilter",
    "lukgamma",
    "maxl",
    "maxq",
    "maxq30",
    "maxq40",
    "maxq50",
    "maxquad",
    "mxhilb",
    "oet5",
    "oet6",
    "osborne2",
    "pbc1",
    "polak2",
    "polak3",
    "polak6",
    "prob10",
    "prob102",
    "prob107",
    "prob109",
    "prob110",
    "prob113",
    "prob115",
    "prob116",
    "prob206",
    "prob208",
    "prob210",
    "rosen",
    "shelldual",
    "shor",
    "steiner2",
    "tr48",
    "transformer",
    "watson",
    "wong1",
    "wong2",
    "wong3"
]

@pytest.mark.parametrize('problem_name', PROBLEM_NAMES)
def test_dfn_dfl(problem_name):
    problem_module = importlib.import_module(f'{PROBLEM_FOLDER}.{problem_name}')
    
    problem_args = {
        'name': problem_module.name,
        'startp': problem_module.startp,
        'lb': problem_module.lb,
        'ub': problem_module.ub,
        'nint': problem_module.nint,
        'ncont': problem_module.ncont,
        'lbmix': problem_module.lbmix,
        'ubmix': problem_module.ubmix,
        'x_initial': problem_module.x_initial,
        'feval': problem_module.feval,
    }

    problem = Problem(**problem_args)
    
    args = Namespace(
        alg='DFN_DFL',
        NM_memory=3,
        max_fun=5000,
        outlev=0,
        constrained=False
    )

    if problem.nint>=2 and problem.ncont>=2:
        try:
            results = run_problem(args, problem)
        except Exception as e:
            pytest.fail(f"Running problem {problem_name} failed with exception: {e}")
        
        assert results is not None
        assert len(results) > 0
        for result in results:
            assert 'algorithm' in result
            assert 'problem' in result
            assert 'constraint_type' in result
            assert 'n' in result
            assert 'm' in result
            assert 'f' in result
            assert 'stopfl' in result
            assert 'CACHE.F_length' in result
            assert 'nf' in result
            assert 'x' in result
    else:
        with pytest.raises(Exception) as exc_info:   
            results = run_problem(args, problem)
        assert str(exc_info.value) == "The problem should contain at least 2 integer variables and 2 continuous variables."