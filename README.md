# nabla-c ∇

> Optimization and machine learning built from scratch in C, from numerical derivatives to a neural network, one stage at a time.

No external libraries beyond the C standard library and `libm`. Every stage reuses the previous one, so the whole project is a single path from calculus to a working neural network.

## Why

Most ML tutorials start with a framework. This project starts with a partial derivative and asks: what is the minimum I need to build next? Each stage ends with something that runs on its own.

## Roadmap

- [ ] **Stage 1: Optimizers.** Numerical gradient and Hessian, Newton's method, gradient descent (fixed step, momentum, backtracking line search), compared on the same test functions.
- [ ] **Stage 2: Linear regression.** Cost function over a dataset, gradient descent vs. the closed-form solution (normal equations).
- [ ] **Stage 3: Logistic regression / softmax.** Classification, cross-entropy, analytic gradients checked against the numerical ones.
- [ ] **Stage 4: Autodiff.** Small reverse-mode engine; the numerical gradient from stage 1 becomes the gradient checker.
- [ ] **Stage 5: MLP.** Dense layers, backpropagation, SGD and Adam. First goal: XOR.
- [ ] **Stage 6: Digit classifier.** MNIST in pure C, with train / save / load from the CLI.

## Stage 1 layout

```
nabla-c/
├── Makefile
├── README.md
├── .gitignore
├── include/
│   └── nabla/
│       ├── matrix.h          # matrix_t, basic ops, LU solve
│       ├── calculus.h        # numerical gradient and Hessian
│       ├── objective.h       # objective_fn signature + test functions
│       └── optim.h           # Newton, gradient descent, line search
├── src/
│   ├── matrix.c
│   ├── calculus.c
│   ├── objective.c           # f, g, h, Rosenbrock, ...
│   ├── newton.c
│   ├── gradient_descent.c
│   └── line_search.c
├── apps/
│   └── optimize.c            # CLI: pick function, method, start point
├── tests/
│   ├── test_matrix.c
│   ├── test_calculus.c       # numerical vs. known analytic derivatives
│   └── test_optim.c          # convergence on known minima
├── scripts/
│   └── plot_paths.py         # reads docs/*.csv, plots trajectories
└── docs/
    └── .gitkeep              # generated CSVs go here (ignored by git)
```

## Design notes

- **Objective signature:** `double f(const double *x, size_t n, void *ctx)`. The `ctx` pointer carries data or parameters, so the same optimizers work later for regression and networks.
- **Gradient:** central differences, `h ≈ 1e-6`. **Hessian:** `h ≈ 1e-4`.
- **Newton step:** solves `H · dx = -∇f` with LU and partial pivoting instead of building the inverse.
- **Stopping:** `‖∇f‖ < tol` or maximum iterations, whichever comes first.
- **Output:** each run writes the trajectory (`iter, x1..xn, f`) and a sampled grid of `f` to `docs/` so it can be plotted.

## Build

```sh
make
make test
make clean
```

Requires a C11 compiler (`gcc` or `clang`) and `make`.

## Usage

Target interface for stage 1 (adjust as the CLI takes shape):

```sh
./build/bin/optimize --func rosenbrock --method gd --lr 0.001 --iters 5000 --x0 -1.2,1.0
./build/bin/optimize --func rosenbrock --method newton --iters 50 --x0 -1.2,1.0
python3 scripts/plot_paths.py docs/data.csv docs/malla.csv
```

## Status

Work in progress. Stage 1 is built on an existing n-dimensional Newton's method implementation with hand-written partial derivatives.

## License

MIT