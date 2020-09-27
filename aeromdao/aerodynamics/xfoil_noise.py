import numpy as np
import xfoil
import subprocess as sp

class XFoil(xfoil.XFoil):

    def __init__(self):
        super().__init__()

    def xfoil_noise(self, sound_speed = 330, ref_pressure = 101325):
        # Noise computations
        xs, cps = self.get_cp_distribution()
        noise = np.average(noise_level(pressure_from_cp(cps, self.M * sound_speed, ref_pressure, sound_speed)))

        return noise

def xfoil_obj(polar_obj, airfoil, mach, reynolds, max_iter = 20):
    
    polar_obj.max_iter = max_iter
    polar_obj.airfoil = airfoil
    polar_obj.M = mach
    polar_obj.Re = reynolds
    
    return polar_obj

def pressure_from_cp(cp, speed, ref_pressure=101325, rho=1.225): return (cp*0.5*rho*speed**2)

def noise_level(press, ref_press=2e-5): return 20*np.log10(np.abs(press)/ref_press)       

# Reynolds number calculator
def reynolds_number(chord_length, freestream_velocity, rho = 1.225, mu = 1.7894e-5, water=False):
    if water:
        rho = 998.2
        mu = 8.899e-4
    return rho * freestream_velocity * chord_length / mu