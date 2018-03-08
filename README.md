# FLSA

[![Build Status](https://travis-ci.org/EQt/FLSA.jl.svg?branch=master)](https://travis-ci.org/EQt/FLSA.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/e28l9al5h3r0hmcu/branch/master?svg=true)](https://ci.appveyor.com/project/EQt/flsa-jl/branch/master)

**At the moment, the code is being refactored!!!**

Computing a graph induced *Fused LASSO Signal Approximator*.
You can use it to denoise data.
The package includes some utility methods to assess the algorithms and
apply it to *images* as well es *ion-mobilty spectrometry* (IMS) data sets.

## Mathematical Formulation

![flsa formula](resources/flsa-formula.svg?sanitize=true "objective function")

For the one dimensional version of the Johnson's dynamic programming algorithm, have a look into
[Lasso.jl](https://github.com/simonster/Lasso.jl)

## Denoising
The fused LASSO signal approximator can be used to denoise e.g. images:
#### Noisy Input
![demo noise](resources/demo_noise.png?raw=true "noisy input data")

#### Cleaned by FLSA
![demo flsa](resources/demo_flsa.png?raw=true "after cleaning with FLSA")


## Algorithms

### Fast Gradient Projection (FGP)
Also known as Fast Iterative Shrinkage Algorithm (FISTA).

### Maximum Gap Tree (MGT)
Own algorithm based on a iterative approximation by dynamic programming algorithm minimizing a sub-tree-graph.


## Copyright

[Elias Kuthe](mailto:elias.kuthe@tu-dortmund.de) 2015, 2016, 2017, 2018.
This work is provided under the simplified BSD license (see [LICENSE.md](/LICENSE.md)).
