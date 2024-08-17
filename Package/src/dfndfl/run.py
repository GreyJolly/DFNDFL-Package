import sys
import os
import argparse
import importlib.util
import logging

current_dir = os.path.dirname(os.path.abspath(__file__))
top_level_dir = os.path.abspath(os.path.join(current_dir, '..'))
if top_level_dir not in sys.path:
    sys.path.append(top_level_dir)

from dfndfl.problem_runner import Cache, Problem
from dfndfl import mixed_DFN_DFLINT as dfn_dflint

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run a DFO')
    parser.add_argument('--problem-path', type=str, required=True, help='Path to the Python file containing the problem definition')
    parser.add_argument('--constrained', action='store_true', help='Solve constrained problem')
    parser.add_argument('--alg', type=str, default='DFN_DFL', choices=["DFN_DFL", "DFL"], help='Algorithm to be used')
    parser.add_argument('--max_fun', type=int, default=5000, help='Maximum number of function evaluations')
    parser.add_argument('--outlev', type=int, default=1, help='Output level')
    parser.add_argument('--NM_memory', type=int, default=3, help='History size for nonmonotone linesearch')
    return parser.parse_args()

def load_problem_from_path(path):
    spec = importlib.util.spec_from_file_location("problem_module", path)
    problem_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(problem_module)

    required_vars = ['name', 'startp', 'lb', 'ub', 'nint', 'ncont', 'lbmix', 'ubmix', 'x_initial', 'feval', 'xmix']
    problem_args = {}

    for var in required_vars:
        if hasattr(problem_module, var):
            problem_args[var] = getattr(problem_module, var)
        else:
            raise ValueError(f"Required variable '{var}' not found in the problem file.")

    problem = Problem(**problem_args)
    
    return problem

def run_problem(args, problem):
    alg = args.alg
    M = args.NM_memory
    max_fun = args.max_fun
    outlev = args.outlev
    results = []
    
    if problem.nint >= 2 and problem.ncont >= 2:
        for J, m in problem.m.items():
            if (args.constrained and m > 0) or (not args.constrained and m == 0):
                CACHE = Cache(max_fun, problem.n + m + 1)
                x, f, stopfl, Dused, nf = dfn_dflint.box_DFN_DFL(alg, M, J, problem, CACHE, max_fun, outlev)
                result = {
                    "algorithm": alg,
                    "problem": problem.name,
                    "constraint_type": J,
                    "n": problem.n,
                    "m": problem.m[J],
                    "f": f,
                    "stopfl": stopfl,
                    "CACHE.F_length": len(CACHE.F),
                    "nf": nf,
                    "x": x
                }
                logger.info(result)
                results.append(result)
        return results
    else:
        raise ValueError("The problem should contain at least 2 integer variables and 2 continuous variables.")

def main():
    args = parse_arguments()
    problem = load_problem_from_path(args.problem_path)
    run_problem(args, problem)

if __name__ == "__main__":
    main()