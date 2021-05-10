# -*- coding: utf-8 -*-
"""
Created on Tue May  4 15:40:02 2021

@author: marti
"""
import numpy as np
import matplotlib.pyplot as plt

# takes an input of height in meters and returns the electron density of the atmosphere
# at that height in number of electrons per cubic meter
def density_of_atmosphere(height):
    temperature = np.zeros(len(height))
    pressure = np.zeros(len(height))
    temperature[height < 11000] = 15.04 - 0.00649 * height[height < 11000]
    pressure[height < 11000] = 101.29 * ((temperature[height < 11000] + 273.1) / 288.08) ** 5.256
    temperature[height > 11000] = -56.46
    pressure[height > 11000] = 22.65 * (2.71828 ** (1.73 - 0.000157 * height[height > 11000]))
    density = pressure / (0.2869 * (temperature + 273.1))
    return density

#plt.plot(np.linspace(1000,20000,10000), density_of_atmosphere(np.linspace(1000,20000,10000)))
#plt.show()