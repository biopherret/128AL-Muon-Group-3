import numpy as np
import random
import matplotlib.pyplot as plt

def find_decay_probability(t_prime):
    """Given the the time it takes for the muon to reach the detector in its own reference frame, returns the probability that it decays. """
    tau = 2.2*10**-6 #seconds
    lam = 1 / tau
    decay_prob = np.exp(-1 * lam * t_prime) 
    
    return decay_prob

def is_detectable_energy(energy):
    """Given the final energy of a single muon, returns True if the energy is in the detectable range"""
    detection_energy = 160 #MeV
    energy_allowance = 20 #MeV

    if (energy <= detection_energy + energy_allowance) and (energy >= detection_energy - energy_allowance):
        return True
    else:
        return False

def particle_does_not_decay(t_prime):
    """Given the time the muon travels in its own reference frame, returns True if the muon does not decay
    probabilistic in nature, uses random, and will vary in output"""
    random_num = random.random()
    decay_prob = find_decay_probability(t_prime)

    if random_num > decay_prob:
        return True
    else:
        return False

def is_detected(energy, t_prime):
    """Given a muon's final energy and total time traveled in its own reference frame, returns True if the muon is detected by the detector"""
    if is_detectable_energy(energy) and particle_does_not_decay(t_prime):
        return True
    else:
        return False

#dummy lists for now for testing purposes
zenith_angles = np.linspace(0,90,100)

detected_angles_flat = []
detected_angles_round = []
for muon in range(len(E_flat_final)):
    #find the angles of detected muons for flat earth
    if is_detected(E_flat_final[muon], t_prime_flat_final[muon]):
        detected_angles_flat.append(zenith_angles[muon])
    #find the angles of detected muons for round earth
    if is_detected(E_round_final[muon], t_prime_round_final[muon]):
        detected_angles_round.append(zenith_angles[muon])

bin_number = 20       

plt.hist(detected_angles_flat, bins = bin_number)
plt.title('Flat Earth Distribution of Detected Muons as a Function of Zenith Angle')
plt.xlabel('Zenith Angle (Radians)')
plt.ylabel('Number of Detected Muons')
plt.show()

plt.hist(detected_angles_round, bins = bin_number)
plt.title('Round Earth Distribution of Detected Muons as a Function of Zenith Angle')
plt.xlabel('Zenith Angle (Radians)')
plt.ylabel('Number of Detected Muons')
plt.show()