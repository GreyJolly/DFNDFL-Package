# DFNDFL Python Package
The DFNDFL Python package provides tools for solving Direct Function Optimization (DFO) problems. The package allows users to define custom optimization problems and solve them using the algorithms DFNDFL or DFL.

## Installation
The package can be installed with pip:
```
python3 -m pip install dfndfl
```
## Usage
### Problem definition
To use the DFNDFL package, first define an optimization problem in a Python file. The problem file should define the following variables:

- `name`: string - name of the problem
- `startp`: numpy array - the starting point for the continuous problem
- `lb`: numpy array - the lower bounds of the continuous problem
- `ub`: numpy array - the upper bounds of the continuous problem
- `nint`: int - the number of discrete variables (>= 2)
- `ncont`: int - the number of continuous variables (>= 2), (the variables are so intended: x[0] ... x[ncont-1] are the continuous variables and x[ncont] ... x[n-1]     are the discrete variables)
- `lbmix`: numpy array - the actual lower bounds of the mixed integer problem
- `ubmix`: numpy array - the actual upper bounds of the mixed integer problem
- `x_initial`: numpy array - the actual initial point of the mixed integer problem
- `xmix`: numpy array - a temporary array used for computationv
- `feval`: function handle - function to compute the objective function value
(N.B. the point must be reconstructed through the use of reconstruct_xmix before calling feval)

Here's a sample problem definition:
```
import numpy as np

name      = 'MISO prob. 10'
nint   = 30
ncont  = n-nint
lb     =-15.0*np.ones(n)
ub     = 30.0*np.ones(n)
lbmix  =-15.0*np.ones(n)
ubmix  = 30.0*np.ones(n)
startp =   7.0*np.ones(n) 
x_initial =7.0*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):  
    f = - 20*np.exp(-0.2*np.sqrt(np.sum(x**2)/15)) - np.exp(np.sum(np.cos(2*np.pi*x))/15)
    return f
```

### Running a problem
To run a problem from the command line, use:
```
python3 -m dfndfl.run --problem-path path/to/your/problem.py [options]
```

### Command-Line Arguments

The `run.py` script accepts several command-line arguments to control the optimization process:

 ```
usage: run --problem-path PROBLEM_PATH [-h] [--constrained] [--alg {DFN_DFL,DFL}] [--max_fun MAX_FUN]
[--outlev OUTLEV] [-M NM_MEMORY] 

--problem-path PROBLEM_PATH
						Path to the Python file containing the problem definition

optional arguments:
-h, --help			show this help message and exit
--constrained		Solve constrained problems
--alg {DFN_DFL,DFL}	Name of algorithm to be used
--max_fun MAX_FUN	Maximum number of function evaluations
--outlev OUTLEV		Output level
-M NM_MEMORY, --NM_memory NM_MEMORY
						History size for nonmonotone linesearch
```

## Logging

The DFNDFL package uses Python's `logging` module to log results and progress. The log level is set to `INFO` by default. You can adjust the log level by modifying the `logging.basicConfig(level=logging.INFO)` line in `run.py`.

## License
This project is licensed under the GNU General Public License - see the [LICENSE](./LICENSE) for details.

## Useful links
- https://github.com/GreyJolly/DFNDFL-Package
- https://pypi.org/project/dfndfl
- https://github.com/DerivativeFreeLibrary/DFNDFL
