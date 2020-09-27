# %%
import numpy as np
from xfoil import model
from aeromdao import xfoil_noise as xfn
from aeromdao import airfoil_parametrization as afp
from smt import sampling_methods

# %% LHS sampling
limits = np.array([
                   [0.15, 0.4], # Upper surface LE
                   [-0.025, 0.3], # Upper surface TE
                   [-0.4, -0.1], # Lower surface LE pt 2
                  ])

num = 50
sampling = sampling_methods.LHS(xlimits=limits)
xs = [ np.round(x, 6) for x in sampling(num) ]

# %% Computations
polar = xfn.XFoil()
data, airfoils = [], []

# Flow conditions
alpha = 8
mach = 0.4
reynolds = 2e6

for (n, (au0, au4, al1)) in enumerate(xs):

    # Compute airfoil
    foil = afp.airfoil_CST([au0, 0.3, 0.2, 0.15, au4], 
                           [-0.2, al1, -0.1, -0.001, -0.02],
                           0.001, 0.001 - 0.001, -0.2)
    airfoil = model.Airfoil(foil[:,0], foil[:,1])
    analysis = xfn.xfoil_obj(polar, airfoil, mach, reynolds)
    result = analysis.a(alpha)

    airfoils.append(foil)
    data.append(result)

# %% 
results = filter(lambda x: np.nan not in x[0], zip(data, airfoils))

# %%
