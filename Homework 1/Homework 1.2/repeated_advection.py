#Exercise 2: This program implements one of two possible schemes for the advection equation, repeating the simulation for 
# different values of the courant factor. The output is a plot of the evolution of the norm for the different values of the courant factor
import numpy as np
import matplotlib.pyplot as plt
#PARAMETERS OF THE SIMULATION
#velocity
a = 1.0
#simulated domain
L = 10
T = 20
#Number of points for space (J) and time (N)
J = 101

#The variable "method" selects which method will be applied during the simulation. The available values for this variable are:
#"LAX_FRIEDRICHS","LAX_WENDROFF". For any other value the code will display an error message and the execution
#will be terminated
method = "LAX_FRIEDRICHS"
if(not(method in ["LAX_FRIEDRICHS","LAX_WENDROFF"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()

#Function defining the step function. np.where acts as an "elementwise if": the first argument is an array of boolean values
#(representing the regions where the step-function returns 1 or 0) and the output is an array of the same size as the original one 
# where if the boolean value is true (i.e.: the condition is satisfied) for a given index, then the output array will have in that 
# index the value of the second argument (1), otherwise the value of the third (0)
def step_function(arg):
    return np.where((arg>4) & (arg<6),1,0)

for c_f in [0.3,0.5,0.8]:
    #Each method results in a different scheme. Here we set the function that implements the scheme for each of the possible choices
    if(method == "LAX_FRIEDRICHS"):
        def advance(space_arr):
            return 0.5*(np.roll(space_arr,-1) + np.roll(space_arr,+1)) - (a*dt/(2*dx))*(np.roll(space_arr, -1) - np.roll(space_arr, +1))
    elif(method == "LAX_WENDROFF"):
        def advance(space_arr):
            return space_arr - 0.5*(a*dt/dx)*(np.roll(space_arr,-1) - np.roll(space_arr,+1)) + 0.5*((a*dt/dx)**2)*(np.roll(space_arr,-1) - 2*space_arr + np.roll(space_arr,+1))
        

    #Lattice spacing
    dx = L/(J-1)
    dt = c_f*dx/a
    #Define frames number (used to stop the simulation)
    frames_number = int(T/dt)

    x = np.linspace(0,10,J)

    u_init = step_function(x)
    u_old = u_init.copy()
    u_curr = []

    time = 0.0
    frame = 0
    times = []
    norms = []

    while(frame<=frames_number):
        #For all the schemes different form Leapfrog we just apply the standard "advance" function, which implements the whole method
        u_curr = advance(u_old)
        u_old = u_curr.copy()
        norms.append(np.linalg.norm(u_curr))
        frame+=1
        time+=dt
        print("timestep #{}, time: {}".format(frame, time))
        times.append(time)
    
    plt.plot(times,norms, label = r'$\Delta x$' + "={:.3f}".format(dx) + ", " + r'$C_f$' + "={:.3f}".format(c_f))
    plt.legend()

plt.grid()
plt.title("Evolution of the L2 norm")
plt.savefig("norm_{}_repeated.png".format(method))
plt.close()