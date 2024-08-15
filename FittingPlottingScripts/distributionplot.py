# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 19:39:53 2024

@author: Arran
"""

import numpy as np
import matplotlib.pyplot as plt

kappa = 10

x = np.linspace(-1, 1, 100)


def Kappa(x, kappa):
    #defining the variance-like variable indicated in Jeffreys. N. S et al (2016)
    variance_k = 1 * (1 - 3/(2*kappa))
    #argument used in kappa function 
    arg = (1 + ((x - 0)**2)/(2 * (0.5**2) * kappa))
    
    return 1 * (arg) ** (-1 *kappa + 1)



fig = plt.figure()

plt.title("Kappa Distributions,\n Intensity about a central wavelength $\lambda_{0}$")

plt.plot(x, Kappa(x, 2), label = "$\kappa$ = 2")

plt.plot(x, Kappa(x, 5), label = "$\kappa$ = 5")

plt.plot(x, Kappa(x, 10), label = "$\kappa$ = 10")

plt.plot(x, Kappa(x, 20), label = "$\kappa$ = 20")

# plt.plot(x, Kappa(x, 100), label = "$\kappa$ = 100")

plt.plot(x, Kappa(x, 10000), label = r'$\kappa$ $\rightarrow$ $\infty $')

plt.xlabel("$\lambda - \lambda_{0}$")
plt.ylabel("I($\lambda$)")

plt.legend()

plt.savefig("KappaDistrib.pdf")

plt.show()