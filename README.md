# NewtonKrylov

This function implements a Newton-Krylov solver to find the root of a function.

The basic idea is to compute the inverse of the Jacobian with an iterative Krylov method. These methods require only evaluating the Jacobian-vector products, which are conveniently approximated by a finite difference:

Jv≈(f(x+ω∗v/|v|)−f(x))/ω
Jv≈(f(x+ω∗v/|v|)−f(x))/ω
Due to the use of iterative matrix inverses, these methods can deal with large nonlinear problems.

LGMRES (a variant of restarted GMRES iteration) reuses some of the information obtained in the previous Newton steps to invert Jacobians in subsequent steps.

This repository contains a C++ implementation of the Newton-Krylov algorithm, with Python bindings. It functionally equivalent to SciPy's implementation, but is around 20 times faster.

## Installation (Linux)

If you have Python 3 installed, you can simply add `NewtonKrylov/bindings` to your path, and import the precompiled module with `import NewtonKrylov`.
Alternatively, use the following commands to compile the Python bindings yourself (making sure you start in the top NewtonKrylov directory):

```sh
$ sudo apt-get install python3-dev
$ cd bindings
$ bash compile_bindings.sh
```

## Python Usage

The NewtonKrylov module has one function, solve:

`NewtonKrylov.solve(F, xin, f_tol, f_rtol, x_tol, x_rtol)`

### Parameters:	
* F : (function(x) -> f)
Function whose root to find; should take and return an array
* x0 : (array)
Initial guess for the solution
* f_tol : (float, optional)
Absolute tolerance (in max-norm) for the residual. If omitted, default is 6e-6.
* f_rtol : (float, optional)
Relative tolerance for the residual. If omitted, not used.
* x_tol : (float, optional)
Absolute minimum step size, as determined from the Jacobian approximation. If the step size is smaller than this, optimization is terminated as successful. If omitted, not used.
* x_rtol : (float, optional)
Relative minimum step size. If omitted, not used.

### Returns:	
* x : array
The final solution.

### NOTE
Do not call `solve` with `arg=None` explicitly for any arguments `arg`. Simply leave out these arguments.

## C++ Usage

The function to which you want to find the roots must take and return an Eigen VectorXd. x0 must be a VectorXd also.
Your code must be compiled against the eigen3 and pybind11 libraries, found in the `include` folder.
