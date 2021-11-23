# NLEigenJacobiDavidsonModified
Routine to solve the nonlinear eigenvalue problem for real symmetric eigenproblem using the generalzied Jacobi-Davidson algorithm as proposed by [Dumont(2007)](https://doi.org/10.1002/nme.1997). 

It is used the C++ linear algebra and logging libraries [Eigen](https://gitlab.com/libeigen/eigen) and [spdlog](https://github.com/gabime/spdlog), respectively.


# Usage
The nonlinear eingesolver based on the modified Jacobi-Davidson method can be added in your repository by the following command:

```sh
git clone --recursive https://github.com/capnrenan/NLEigenJacobiDavidsonModified
```

# References
Dumont, N.A. (2007), On the solution of generalized non-linear complex-symmetric eigenvalue problems. Int. J. Numer. Meth. Engng., 71: 1534-1568. https://doi.org/10.1002/nme.1997
