"""This file generates a representative sample of muons which will hit a detector at sea level.
The muons have an origin energy, the pathlength they must travel to the detector, and an an
azimuthal angle relative to the detector."""

import numpy as np  
import matplotlib.pyplot as plt


#CONSTANTS
R_earth_SB = 5256 #KM


#HELPER FUNCTIONS. Group 3 just look at the bottom main function gen_muons.py

#Generates normal sample of muon heights with mean 15km and stdev 2km.
def gen_heights(n_muons, mean = 15, stdev = 2):
    heights_list = np.random.normal(mean, stdev, size=n_muons)
    return heights_list
   
#Path length for flat Earth. From https://arxiv.org/pdf/1606.06907v3.pdf (5)
def gen_pathLengthFlat(theta,d):
    return d/np.cos(np.deg2rad(theta))

#Path length for curved Earth. From https://arxiv.org/pdf/1606.06907v3.pdf (7)
def gen_pathLength(theta, d, R):
    return d * ((R**2/d**2 * np.cos(np.deg2rad(theta))**2 + 2 * R / d + 1)**(.5) - R / d * np.cos(np.deg2rad(theta)))

def monte_carlo_sample(f, bounds, n_samples):
    """Generates a sample from a continuous probability distribution."""
    samples = []
    pmax = f(bounds[0])
    tries_per_run = int(n_samples*1/pmax)
    while len(samples) < n_samples:
        x = np.random.rand(tries_per_run)*(bounds[1]-bounds[0])+bounds[0]
        y = np.random.rand(tries_per_run)*pmax
        good = x[y <= f(x)]
        samples = samples + [i for i in x[y <= f(x)]]
    return np.array(np.array(samples))[:n_samples]

def gen_energies(n_muons):
    """Generates the muon initial energies."""
    pdist, bounds = fit_energylaw()
    samples = monte_carlo_sample(pdist, bounds, n_muons)
    return samples


def fit_energylaw(showplots = False):
    """ Fits some data to a power law energy vs intensity distribution.
    Returns
    ----------
    f: lambda x
        Intensity, as a function of energy 
    xbounds: list
        min, max
    """
    #Data is from Cosmlc Ray Muon Spectrum In the Atmoephere M. Circella et al 1993 Fig 4
    #(at 15KM. conversion from depth to altitude using https://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html)
    #Units are GeV/c vs (cm^2 s sr Gev / c) ^ -1
    data = np.array([[.4, .025], [.5, .017], [.7, .01], [1, .008], [1.25, .004], [1.8, .003], [2.5, .0015], [5,.00035], [18, .00001]])
    xbounds = [1, 100]
    #Fit data to ax^b
    data_log = np.log(data)
    fits = np.polyfit(data_log[:,0], data_log[:,1], 1)
    a = np.exp(fits[1])
    b = fits[0]
    if(showplots):
        fitdata = np.polyfit(data_log[:,0], data_log[:,1], 1,cov=True)
        print(fitdata[1])
        x = np.linspace(.4, 50, 1000)
        plt.scatter(data[:,0], data[:,1], label="Data from Circella")
        plt.loglog(x,  a * x **b, color="green", label="ax^b fit")
        plt.xlabel("Muon Energy (GeV/c)")
        plt.ylabel("Differential Intensity (cm^2 s sr Gev / c)^-1")
        plt.title("Fitting Flux vs Energy at 15km from Circella et al.")
        plt.legend()
        plt.show()
    f = lambda x: a * x**b
    return f, xbounds
#fit_energylaw(True)

#MAIN GENERATION FUNCTION

def gen_muons(n_muons: int):
    """ Generates a representative sample of muons.

    Parameters
    ----------
    n_muons: int
        Number of muons to generate. 1E5 is a max reasonable sample size with the current sampling method.
    flat: bool
        Whether or not to use a flat earth model for pathlength. Set true if you're in flat gang.

    Returns
    ----------
    muons: np.array(n_muons, 4)
        Array of muons in the form
            [[muon_1_pathlength_flat (M), muon_1_pathlength_round (M), muon_1_energy (MeV/c), muon_1_angle (deg)],
             [muon_2_pathlength_flat (M), muon_2_pathlength_round (M), muon_2_energy (MeV/c), muon_2_angle (deg)] 
             ...]
    """
    thetas = np.random.uniform(0,90, size=n_muons)
    heights = gen_heights(n_muons)
    pathlengths_flat = gen_pathLengthFlat(thetas, heights).reshape(n_muons,1)
    pathlengths_round = gen_pathLength(thetas, heights, R_earth_SB).reshape(n_muons,1)

    heights = heights.reshape(n_muons,1)
    energies = gen_energies(n_muons).reshape(n_muons,1)
    thetas = thetas.reshape(n_muons,1)

    muons = np.concatenate((np.abs(pathlengths_flat * 1000), np.abs(pathlengths_round * 1000), energies * 1000, thetas), axis=1)
    return muons

"""
muons = gen_muons(1000)
print(muons)
x0_flat_initial = muons[:,0] 
x0_round_initial = muons[:,1] 

E_initial = muons[:,2] 
zenithAngle_final = muons[:,3]

plt.figure(3)
plt.scatter(zenithAngle_final, x0_flat_initial)
#plt.figure(4)
plt.scatter(zenithAngle_final, x0_round_initial)

plt.figure(5)
x = np.linspace(-89.9,89.9,100000)
plt.plot(x, gen_pathLength(x, 15, R_earth_SB))
plt.plot(x, gen_pathLengthFlat(x, 15))
plt.show()
"""