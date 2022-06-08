# NLEigenJacobiDavidsonModified
Routine to solve the nonlinear eigenvalue problem for real symmetric eigenproblem using the generalzied Jacobi-Davidson algorithm as proposed by [Dumont(2007)](https://doi.org/10.1002/nme.1997). 

It is used the C++ linear algebra library, [Eigen](https://gitlab.com/libeigen/eigen), and the C++ logging library, [spdlog](https://github.com/gabime/spdlog). The eigensolver code can also support the [Blaze-lib](https://bitbucket.org/blaze-lib/blaze/src/master/) C++ math library. To enable [Blaze-lib](https://bitbucket.org/blaze-lib/blaze/src/master/), you have to define the macro ```#define ENABLE_BLAZE``` in the pre-compiled header file  ```nlpch.h```. Check the [requirements for Blaze-lib](https://bitbucket.org/blaze-lib/blaze/wiki/Configuration%20and%20Installation).

# Usage
The nonlinear eingesolver based on the modified Jacobi-Davidson method can be added in your repository by the following command:

```sh
git clone --recursive https://github.com/renancsales/NLEigenJacobiDavidsonModified
```

# References
Dumont, N.A. (2007), On the solution of generalized non-linear complex-symmetric eigenvalue problems. Int. J. Numer. Meth. Engng., 71: 1534-1568. https://doi.org/10.1002/nme.1997
