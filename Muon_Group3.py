import numpy as np
import matplotlib.pyplot as plt
import random
import generate_muons
import vectorized_group_one_electron_density
import copy

#Constants: Not subject to change
c=3e8 #Speed of light, m/s
m=1.88e-28/1.6605e-27 #mass of muon, AU
p0=101325 #sea level standard atmospheric pressure, Pa
T0= 288.15 #sea level standard temperature, K
g= 9.807 #m/s^2
L= .0065 #temperature lapse rate, K/m
R= 8.31446 #ideal gas constant, J/(mol K)
M= 0.0289652 #molar mass of dry air, kg/mol
e = 1 #electron charge, e

#changes to the constants, I'm overwriting down here in case we want to go back to the definitions above
m = 105.658 #MeV/c^2

#Constants: Subject to change based on info we get from previous groups
A=13 #atomic weight of the stopping medium
Z=7.32 #atomic number of N
I=A*Z #approximation of mean excitation energy of the stopping medium of atoms
z=1*e #charge of muon relative to electron charge

def get_C0(gamma,Beta):
    """Given the current gamma and beta, returns the updated C0 value
    C0 is defined as the right side of the Bethe formula divided by rho"""
    C0 =  100 * 0.00003071*Z*z*z/(A*Beta*Beta)*(np.log(2*m*c**2*10**6*Beta*Beta*gamma*gamma/I)-Beta*Beta)
    C0[C0 > .2] = 0
    return C0

def get_rho(x):
    """Given the height from the ground, returns the mass density, rho of the air at that height IN Kg/m^3."""
    #rho=1000*p0*M/(R*T0)*(1-L*x/T0)**(g*M/(R*L)-1)
    return vectorized_group_one_electron_density.density_of_atmosphere(x)

def get_dt_prime(C0,rho,gamma1,gamma2):
	dg = gamma1-gamma2
	dt_prime = (m/(c*rho*C0)) * dg * 1./(np.sqrt(gamma1**2-1))
	return dt_prime


def findDecayProbability(t_prime):
    """Given the the time it takes for the muon to reach the detector in its own reference frame, returns the probability that it decays. """
    tau = 2.2*10**-6. #seconds
    return np.exp(-1*t_prime/tau)
    
"""
E = np.linspace(1e-2, 3000, 2000)
gamma = E/(m) + 1
print(gamma)
beta = np.sqrt(1 - 1/(gamma*gamma))
print(beta)
c0 = get_C0(gamma, beta)
dEdx = c0 * 1
plt.figure(7)
plt.loglog(E, dEdx)
plt.show()
"""

#Initial Conditions, will be provided from the previous groups as a list for each condition of all the muons
muons = generate_muons.gen_muons(10000)
#print(muons)
x0_flat_initial = muons[:,0] 
x0_round_initial = muons[:,1] 

E_initial = muons[:,2] 
zenithAngle_final = muons[:,3]

"""
E_initial = np.array([1000, 2000, 3000], dtype = 'f') #units of MeV 
x0_flat_initial = np.array([15000, 15000, 15000], dtype = 'f') #height of troposphere in m
x0_round_initial = np.array([15000, 15000, 15000], dtype = 'f')
zenithAngle_final = np.array([0, 0, 0])
"""

gamma_initial = E_initial/(m)
Beta_initial = np.sqrt(1-1/(gamma_initial**2))

t_prime_flat = 0
t_prime_round = 0


E_flat_final = np.zeros(len(x0_flat_initial))
E_round_final = np.zeros(len(x0_flat_initial))
t_prime_flat_final = np.zeros(len(x0_round_initial))
t_prime_round_final = np.zeros(len(x0_round_initial))

number_of_steps = 500. #dummy amount, will change
dx_flat = x0_flat_initial/number_of_steps
dx_round = x0_round_initial/number_of_steps

#Find final energies and t_prime for flat earth
Beta = Beta_initial
x = x0_flat_initial
gamma = gamma_initial
E_flat = copy.deepcopy(E_initial)

for i in range(int(number_of_steps)):
    C0=get_C0(gamma,Beta)
    rho=get_rho(x)
    dE=C0*rho*dx_flat
    E_flat -= dE
    #print(E_flat)
    #E2 = E1 - dE #feel like this should be adding dE but that results in increasing energy
    gamma1 = gamma
    gamma2 = E_flat/m + 1
    Beta1=Beta
    Beta2=np.sqrt(1-1/(gamma2*gamma2))
    t_prime_flat -= get_dt_prime(C0,rho,gamma1,gamma2)

    gamma = gamma2
    Beta = Beta2
    x -= dx_flat

E_flat_final = E_flat
t_prime_flat_final = t_prime_flat

print("ROUND: ")

#Find the final energies and t_prime for round earth
Beta = Beta_initial
x = x0_round_initial
gamma = gamma_initial
E_round = copy.deepcopy(E_initial)

for i in range(int(number_of_steps)):
    C0=get_C0(gamma,Beta)
    rho=get_rho(x)
    dE=C0*rho*dx_round
    #E1 = E + 0
    #E2 = E1 - dE #feel like this should be adding dE but that results in increasing energy
    E_round -= dE
    #print(E_round)
    gamma1 = gamma
    gamma2 = E_round/m + 1
    Beta1=Beta
    Beta2=np.sqrt(1-1/(gamma2*gamma2))
    t_prime_round -= get_dt_prime(C0,rho,gamma1,gamma2)

    gamma = gamma2
    Beta = Beta2
    x -= dx_round

E_round_final = E_round
t_prime_round_final = t_prime_round

"""
plt.figure(10)
plt.scatter(zenithAngle_final, E_flat_final - E_round_final)
print("EDIFF")
ediff = E_flat_final - E_round_final
print(E_round_final)
print(E_flat_final)
"""

#Simulating how many muons hit the detector:
    #Checks for muons in the proper energy range
    #Makes a probablistic decision whether a paticular muon has decayed
    #If muon has the proper energy and is found not to decay, add to histogram
        
energyDetected = 160 #Mev

energyAllowance = 40 #Mev

 
anglesOfDetectedMuons_flatEarth = []
anglesOfDetectedMuons_roundEarth = []

for muon in range(len(x0_flat_initial)):
    #Check detection for flat Earth
    if (energyDetected - energyAllowance <= E_flat_final[muon] <= energyDetected + energyAllowance):
        #Decide if it decayed
        if (random.random() < findDecayProbability(t_prime_flat_final[muon])):
            anglesOfDetectedMuons_flatEarth.append(zenithAngle_final[muon])
        
    #Check detection for round Earth
    if (energyDetected - energyAllowance <= E_round_final[muon] <= energyDetected + energyAllowance):
        #Decide if it decayed
        if (random.random() < findDecayProbability(t_prime_round_final[muon])):
            anglesOfDetectedMuons_roundEarth.append(zenithAngle_final[muon])
        
bin_number = 10    

#print(anglesOfDetectedMuons_flatEarth)

plt.figure(1)
plt.hist(anglesOfDetectedMuons_flatEarth, bins = bin_number, label="flat", alpha=.5, color="blue")
plt.title('Distribution of Detected Muons as a Function of Zenith Angle')
plt.xlabel('Zenith Angle (Radians)')
plt.ylabel('Number of Detected Muons')
plt.hist(anglesOfDetectedMuons_roundEarth, bins = bin_number, label="round", alpha=.5, color="green")
plt.legend()

"""
plt.figure(5)
plt.scatter(zenithAngle_final, t_prime_round_final)

plt.figure(2)
t_prime_round_final = np.array(t_prime_round_final)
t_prime_round_final = t_prime_round_final[t_prime_round_final < 1]
t_prime_round_final = t_prime_round_final[t_prime_round_final > 0]
t_prime_flat_final = np.array(t_prime_flat_final)
t_prime_flat_final = t_prime_flat_final[t_prime_flat_final < 1]
t_prime_flat_final = t_prime_flat_final[t_prime_flat_final > 0]
plt.hist(t_prime_round_final, bins = bin_number, label="flat", alpha=.5, color="blue")
plt.title('Distribution of Muons as a Function of muon-frame lifetime')
plt.xlabel('Lifetime (s)')
plt.ylabel('Number of Detected Muons')
plt.hist(t_prime_flat_final, bins = bin_number, label="round", alpha=.5, color="green")
plt.legend()
"""
#plt.figure(3)
#plt.hist(E_initial, bins = 50)


print("Detected flat muons: ")
print(len(anglesOfDetectedMuons_flatEarth))
print("Detected round muons: ")
print(len(anglesOfDetectedMuons_roundEarth))
plt.show()
