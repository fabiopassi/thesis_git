""" This module contains all the functions for ... script """

# Importing modules
import numpy as np
import matplotlib.pyplot as plt
from numba import jit, njit, prange

# Functions

@njit
def U(x) :
    a = 0.1
    omega = 4
    return np.sin(omega*x) + a * np.power(x,2)
