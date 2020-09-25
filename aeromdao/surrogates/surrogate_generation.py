import numpy as np
# from modules import xfoil_noise
import os, importlib
import itertools
import matplotlib.pyplot as plt
from smt.sampling_methods import LHS
from smt.surrogate_models import RMTB,KRG

# importlib.reload(xfoil_noise)

# Sample generation via Latin Hypercube Sampling (LHS)
xlimits = np.array([[0.0, 8.0], [0.0, 0.5],[2.0,7.0],[2.0,7.0]])
sampling = LHS(xlimits=xlimits)


num = 200
x = sampling(num)

for i,j in enumerate(x[:,2:3]):
    x[:,2:3][i] = np.round(j,2)
    
for i,j in enumerate(x[:,3:4]):
    x[:,3:4][i] = np.round(j,2)

# plt.plot(x[:, 2], x[:, 3], "o")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.show()

# Analysis script for NACA xx16 foils
def nacaXX16Analysis(aoa,mach,m,p):
    airfoil = f'{int(m)}{int(p)}16'
    re = xfoil_noise.reynoldsNumber(1.,mach,'air')
    results = xfoil_noise.polar(airfoil,
                                re,
                                mach,
                                aoa,
                                'results',
                                ref_pressure=1e5, rho=1.225, sound_speed=330.)
    return list(results.values())

# dataset = []
# for i in range(len(x)):
#     aoa = np.round(x[i][0], 4)
#     mach = np.round(x[i][1], 4)
#     m = x[i][2]
#     p = x[i][3]
#     print(f"Mach Number: {mach}, Angle of Attack: {aoa}, m: {m}, p: {p}\n")
#     vals = nacaXX16Analysis(aoa,mach,m,p)
#     print(f"CL: {vals[2]}, CD: {vals[3]}, CM: {vals[4]}, Noise: {vals[5]}")
#     lolz = [ vals[0], vals[1], m, p, vals[2], vals[3], vals[4], vals[5] ]
#     dataset.append(lolz)

# # Filename to save dataset
# filesave =f'results/lolz-dataset.txt'
# # Saves dataset into a file
# np.savetxt(filesave, np.array(dataset), fmt='%1.6f', delimiter=', ', header='Angle, Mach, m, p, Cl, Cd, Cm, Noise:')

# Filename to load dataset
fileload = f'results/lolz-dataset.txt'
# Loads previous data into NumPy array
loader = np.loadtxt(fileload, delimiter=', ')


# Partitions specifically for the NACA parametrisation
def partitionNACAData(raw):

    deg2rad = np.pi / 180

    xt = np.array(raw[:,:4])
    yt = np.array(raw[:,4:])
    xlimits = np.array([
                        [0,8.0],
                        [0.0,0.5],
                        [2.0, 7.0],
                        [2.0, 7.0]
                        ])

    xt[:, 0] *= deg2rad
    xlimits[0, :] *= deg2rad

    return xt, yt, xlimits

xt, yt, xlimits = partitionNACAData(loader)

def calculateSurrogateNoise(xt, yt, tol=1e-5):
    trainer = KRG(theta0=[tol, tol, tol, tol])
    trainer.set_training_values(xt, yt)
    trainer.train()
    return trainer

surrogateNoise = calculateSurrogateNoise(xt,yt)

def SurrogateN(x): return surrogateNoise.predict_values(x)[0]

print(SurrogateN(np.array([[3.0*np.pi/180, 0.2, 3.5, 4.5]])))

