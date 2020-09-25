#%%
import numpy as np
import xfoil_noise
import airfoil_parametrization as afp
from smt import sampling_methods
import os, importlib
import itertools
import random
from matplotlib import pyplot as plt
from matplotlib import rc, cm
rc('font',**{'family':'serif'})
rc('text', usetex=True)

importlib.reload(afp)
importlib.reload(xfoil_noise)
#%%

# Instantiate airfoil method using Bernstein class function and basis by default
airfoil = afp.CSTAirfoilLEM(afp.bernstein_class, afp.bernstein_basis, -0.001)

# CST coefficient definitions
num_points = 20
num_desvars = 5
alpha_u = np.linspace(0.001, 0.005, num_desvars)
alpha_l = np.linspace(-0.005, -0.01, num_desvars)
dz_u = 0.0002
dz_l = 0.0002

# Foil surface generation test
xs = np.linspace(0, 1, num_points)
upper_surf = [ airfoil.z(x, alpha_u, dz_u) for x in xs ]
lower_surf = [ airfoil.z(x, alpha_l, dz_l) for x in xs ]

# Zipping upper and lower surfaces into coordinates
foil_upper = np.column_stack((xs, upper_surf))
foil_lower = np.column_stack((xs[::-1], lower_surf[::-1]))

# Cosine spacing test
cos_coords = afp.cosine_airfoil(foil_upper, foil_lower, n=80)
afp.write_foil(cos_coords, 'Cosine spacing CST', '../coordinates/cosine_test.dat')

# Plotting
fig1 = plt.figure(1, figsize=(6,1), dpi=300)
plt.plot(xs, upper_surf, '.-', color='orange') 
plt.plot(xs, lower_surf, '.-', color='blue')
plt.plot(cos_coords[:,0], cos_coords[:,1], '.', color='green', zorder=3)
plt.tight_layout()
plt.savefig('airfoil.png', dpi=300)

#%%
# Hicks-Henne
airfoil_path = '../coordinates/E432.dat'
foil_coords = afp.read_foil(airfoil_path)

n = 30
locations = [3, 15, 50 ,60]
alphas = random.sample(list(np.linspace(-0.02, 0.02, num=20)), len(locations))

old_coords, coords = afp.hicks_henne2D(foil_coords, locations, alphas, bump_width=8)
    
fig1 = plt.figure(1, figsize=(6,1), dpi=300)
plt.plot(old_coords[0], old_coords[1], '.', coords[0], coords[1], '.', markersize=0.7)
plt.legend(['Baseline', 'Perturbed'])

# out = np.column_stack((coords[0], coords[1]))
# export_path = airfoil_path.replace(".dat","-perturbed.dat")
# np.savetxt(export_path, out, header='perturbed', comments='', fmt='%1.6f %10.6f')

#%%

# Output directory to write dummy temp files and final results
out_path = 'results'

# Design space definition for NACA foils
digits = '24'
angles = np.linspace(0., 9., 5)
machs = np.linspace(0.1, 0.6, 5)
reynolds = [ xfoil_noise.reynoldsNumber(1., mach, 'air') for mach in machs ]

dataset = []

# Computes aerodynamic coefficients and noise in Cartesian product of sets
for (ang, mach, re) in itertools.product(angles, machs, reynolds):
    print(f"Angle: {ang}, Mach: {mach}, Reynolds: {re}")
    results = xfoil_noise.polar(f'{digits}16',
                                np.round(re, 6),
                                np.round(mach, 3),
                                np.round(ang, 4),
                                ref_pressure=1e5, rho=1.225, sound_speed=330.)
    dataset.append(list(results.values()))

# Filename to save dataset
file=f'{out_path}/NACA-{digits}16-dataset.txt'
# Saves dataset into a file
np.savetxt(file, np.array(dataset), fmt='%1.6f', delimiter=', ', header='Angle, Mach, Cl, Cd, Cm, Noise:')

#%%

# Filename to load dataset
file = f'{out_path}/NACA-{digits}16-dataset.txt'
# Loads previous data into NumPy array
loader = np.loadtxt(file, delimiter=', ')

# %%
# Sampling with CST method

# LHS sampling
limits = np.array([
                   [0.15, 0.4], # Upper surface LE
                   [-0.025, 0.3], # Upper surface TE
                   [-0.4, -0.1], # Lower surface LE pt 2
                  ])

num = 50
sampling = sampling_methods.LHS(xlimits=limits)
xs = [ np.round(x, 6) for x in sampling(num) ]

# Create file
header = f"a_u_0, a_l_1, a_l_4"
filename = "../testing/CSTAirfoils.txt"
np.savetxt(filename, xs, delimiter=', ', fmt='%.5f', header=header)


# Flow conditions
mach = 0.4
reynolds = 2e6
alpha = 8

storage, airfoils = [], []

for (n, (au0, au4, al1)) in enumerate(xs):

    # Compute airfoils
    airfoil = afp.airfoil_CST([au0, 0.3, 0.2, 0.15, au4], 
                              [-0.2, al1, -0.1, -0.001, -0.02], 
                              0.001, 0.001 - 0.001, -0.2)
    
    # File-naming and writing
    foilname = f'../testing/CST-{n}'
    filename = f'{foilname}.dat'
    afp.write_foil(airfoil, foilname, filename=f'{foilname}.dat')

    # XFOIL computations
    data = xfoil_noise.polar(filename, 
                             reynolds, mach, alpha, 
                             refine=True, graphics=True, max_iter=500)
    
    # Convergence check and append
    if data['cd'] != 0: 
        storage.append(data)
        airfoils.append(airfoil) 
    else: 
        print(f"Airfoil CST-{n} not converged.")
        del xs[n] # Delete unconverged sample point
    
#%%
# Plotting
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')

qoi = 'cd'
qoi_vals = np.array([ data[qoi] for data in storage if data[qoi] ])
qoi_max, qoi_min = max(qoi_vals), min(qoi_vals)
normColorIndex = (qoi_vals - qoi_min)/(qoi_max - qoi_min)
bowcolors = cm.viridis(normColorIndex)

fig = plt.figure(1, figsize=(12,5), dpi=300)
ax = fig.add_subplot(121,
                     projection='3d',
                    )
ax.axis(option='equal')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'Airfoil ID')
# ax.set_title('Foil Profiles')
ax.view_init(60, 240)

# fig2 = plt.figure(2, figsize=(12,5), dpi=300)
ax2 = fig.add_subplot(122, projection='3d')
ax2.set_xlabel(r'$\alpha_{u_0}$')
ax2.set_ylabel(r'$\alpha_{u_4}$')
ax2.set_zlabel(r'$\alpha_{l_1}$')
ax2.set_title(r'Normalized Drag Coefficient, $\overline C_d$')

for (n, (x, airfoil, color)) in enumerate(zip(xs, airfoils, bowcolors)):
    
    # Plot airfoil
    ax.plot(airfoil[:,0], airfoil[:,1], 
            n, 
            '-', markersize=1, linewidth=0.5, 
            label=f'{n}', color=color)

    # this method is called for each point
    # ax.text(airfoil[0,0], airfoil[0,1], n, f'{n}', size=5, zorder=3)

    # Plot noise
    p2 = ax2.scatter(*x, c=color, cmap='viridis') 
    ax2.text(*x, f'{n}', size=6, zorder=1,)

clb = fig.colorbar(p2)
clb.ax.set_title(r'$\overline C_d$')

fig.tight_layout()
fig.savefig('profiles.png', format='png', dpi=300)
# fig2.legend()
# fig2.savefig('design_space.png', format='png', dpi=300)

fig3 = plt.figure(3, dpi=300)
ax = fig3.add_subplot(111,
                    # projection='3d',
                    )
ax.set_aspect('equal')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
[ ax.plot(airfoil[:,0], airfoil[:,1], '-', markersize=1, linewidth=0.5, 
            label=f'{n}') for airfoil in airfoils ]
fig3.show()

# %%
