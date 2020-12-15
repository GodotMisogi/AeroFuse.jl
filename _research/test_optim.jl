##
using Optim

##
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

##
x0 = [0.0, 0.0]
optimize(f, x0)

##
optimize(f, x0, LBFGS(); autodiff = :forward)