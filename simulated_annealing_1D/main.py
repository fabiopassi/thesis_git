# *     1-D simulated annealing     *

# Importing modules
from utils.functions import *
import numpy as np
import matplotlib.pyplot as plt
import time
import sys



# Variables for the system boundaries and move
x_min = -10                                                 # Lower extremum
x_max = 10                                                  # Upper extremum
max_move = x_max / 10                                       # Maximum possible move

# Variables for the animated plot
dx_plot = 0.01                                              # Spacing on the axis
x_axis = np.arange(x_min, x_max, dx_plot)                   # x axis
U_plot = U(x_axis)                                          # Potential curve for the plot

# Variables for SA procedure
MC_steps = 40000                                            # Number of MC steps
tau = - MC_steps / np.log(0.001)                            # Decay time

# Lists for static plot
T_static_plot = np.zeros(MC_steps)
x_static_plot = np.arange(0, MC_steps, 1)

# Define initial position
x_0 = (x_max - x_min) * np.random.rand() + x_min            # Initial position picked at random
U_0 = U(x_0)


# Estimate initial temperature (i.e. accept a move with a probability of 75%)
init_acc_prob = 0.75
num_moves = 1000
deltas = np.zeros(num_moves)

for i in range(num_moves) :
    # Generate a random point
    rand_cand = (x_max - x_min) * np.random.rand() + x_min
    U_cand = U(rand_cand)
    # Generate a random neighbour
    rand_neigh = rand_cand + 2 * max_move * np.random.rand() - max_move
    U_neigh = U(rand_neigh)
    # Calculate the uphill transition cost
    deltas[i] = np.fabs(U_cand - U_neigh)    

T_0 = - np.mean(deltas) / np.log(init_acc_prob)

# Status message
print("Starting SA procedure with the following parameters :\n")
print(f"\t-> MC_steps = {MC_steps}")
print(f"\t-> tau = {tau:.2f}")
print(f"\t-> T_0 = {T_0:.3f}")
print(f"\t-> acc_prob = {init_acc_prob:.3f}")
print(f"\t-> x_0 = {x_0:.3f}\n")


# SA procedure
x_SA = [x_0]
U_SA = [U_0]

for step in range(MC_steps) :

    # Save potential and temperature
    T_static_plot[step] = T_0 * np.exp(-step/tau)

    # Propose a move
    x_new = x_0 + 2 * max_move * np.random.rand() - max_move

    # Calculate the potential
    U_new = U(x_new)

    # Metropolis rule
    if U_new < U_0 :
        x_0 = x_new
        U_0 = U_new
    else :
        T = T_0 * np.exp(-step/tau)
        acc_prob = np.exp(-(U_new - U_0) / T)
        if np.random.rand() < acc_prob :
            x_0 = x_new
            U_0 = U_new
    
    # Append data to list of positions
    x_SA.append(x_0)
    U_SA.append(U_0)


# Animated plot
# to run GUI event loop
plt.ion()
 
# here we are creating sub plots
figure, ax = plt.subplots(figsize=(10, 8))
line1, = ax.plot(x_axis, U_plot)
point, = ax.plot(x_SA[0], U_SA[0], marker = "o", color = "red", markersize = 9, label = f"T = {T_0:.3f}")
ax.legend(handles = [point])

 
# setting title
plt.title("Simulated Annealing", fontsize=20)
 
# setting x-axis label and y-axis label
plt.xlabel("x")
plt.ylabel("U")
 
# Loop
step_plot = 30
for i in range(0, MC_steps, step_plot):

    # updating data values
    line1.set_xdata(x_axis)
    line1.set_ydata(U_plot)
    point.set_xdata(x_SA[i])
    point.set_ydata(U_SA[i])
    point.set_label(f"T = {T_0 * np.exp(-i/tau):.3f}")
    ax.legend(handles = [point])
 
    figure.canvas.draw()

    figure.canvas.flush_events()
 
    time.sleep(0.001)

# Static plot of temperature + cost function
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))

ax1.plot(x_static_plot, T_static_plot)
ax1.set_xlabel("MC step", fontsize=20)
ax1.set_ylabel("T", fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=18)

ax2.plot(x_static_plot, U_SA[:-1], color="red")
ax2.set_xlabel("MC step", fontsize=20)
ax2.set_ylabel("S", fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=18)

plt.show()