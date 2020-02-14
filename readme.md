# Mixed Precision Path Tracking for Polynomial Homotopy Continuation

This repository contains the reference implementation accompanying the article
[arXiv:1902.02968](https://arxiv.org/abs/1902.02968) as well as the
scripts and data to perform the experiments.

This project requires Julia 1.3. In the terminal switch to this folder
and execute
```sh
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
```
to install all necessary dependencies.
To run the experiments see the `scripts` folder.
