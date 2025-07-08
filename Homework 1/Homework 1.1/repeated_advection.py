#Exercise 1: This program implements one of four possible schemes for the advection equation, repeating the simulation for 
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
N = 101#Relevant only for FTCS

#The variable "method" selects which method will be applied during the simulation. The available values for this variable are:
#"FTCS","LAX_FRIEDRICHS","LAX_WENDROFF","LEAPFROG". For any other value the code will display an error message and the execution
#will be terminated
method = "LAX_WENDROFF"
if(not(method in ["FTCS","LAX_FRIEDRICHS","LEAPFROG","LAX_WENDROFF"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()

for c_f in [0.3,0.5,0.8]:#VALUES OF THE COURANT FACTOR TO USE
    #Each method results in a different scheme. Here we set the function that implements the scheme for each of the possible choices
    if(method == "FTCS"):
        def advance(space_arr):
            return space_arr - (a*dt/(2*dx))*(np.roll(space_arr, -1) - np.roll(space_arr, +1))
    elif(method == "LAX_FRIEDRICHS"):
        def advance(space_arr):
            return 0.5*(np.roll(space_arr,-1) + np.roll(space_arr,+1)) - (a*dt/(2*dx))*(np.roll(space_arr, -1) - np.roll(space_arr, +1))
    elif(method == "LEAPFROG"):
        #For the leapfrog methods two functions are defined. Since this scheme requires the information about three subsequent instants,
        #the first step is done with the Lax-Friedrichs scheme, while all the others are done with the Leapfrog scheme
        def advance_LF(space_arr):
            return 0.5*(np.roll(space_arr,-1) + np.roll(space_arr,+1)) - (a*dt/(2*dx))*(np.roll(space_arr, -1) - np.roll(space_arr, +1))

        def advance_Leapfrog(space_arr,space_arr_older):
            return space_arr_older - (a*dt/dx)*(np.roll(space_arr, -1) - np.roll(space_arr, +1))
    elif(method == "LAX_WENDROFF"):
        def advance(space_arr):
            return space_arr - 0.5*(a*dt/dx)*(np.roll(space_arr,-1) - np.roll(space_arr,+1)) + 0.5*((a*dt/dx)**2)*(np.roll(space_arr,-1) - 2*space_arr + np.roll(space_arr,+1))
        

    #Lattice spacing
    dx = L/(J-1)
    #The dime discretization is defined in two different ways depending on if the chosen method is FTCS or not. In the latter case
    #We use the courant factor to define dt
    if(method == "FTCS"):
        dt = T/(N-1)
    else:
        dt = c_f*dx/a
    #Define frames number (used to stop the simulation)
    frames_number = int(T/dt)

    x = np.linspace(0,10,J)
    x0 = 5
    u_init = np.exp(-(x-x0)**2)
    u_old = u_init.copy()
    u_curr = []

    time = 0.0
    frame = 0
    times = []
    norms = []

    while(frame<=frames_number):
        if(method == "LEAPFROG"):
            #In the case of the Leapfrog method we need to define two separate steps: For the first frame we apply the Lax-Friedrichs scheme
            #while for all the other frames we apply the ordinary Leapfrog method
            if(frame == 0):
                print("first step\n")
                u_curr = advance_LF(u_old)
                norms.append(np.linalg.norm(u_curr))
            else:
                u_new = advance_Leapfrog(u_curr,u_old)
                u_old = u_curr.copy()
                u_curr = u_new.copy()
                norms.append(np.linalg.norm(u_new))
        else:
            #For all the schemes different form Leapfrog we just apply the standard "advance" function, which implements the whole method
            u_curr = advance(u_old)
            u_old = u_curr.copy()
            norms.append(np.linalg.norm(u_curr))
        frame+=1
        time+=dt
        print("timestep #{}, time: {}".format(frame, time))
        times.append(time)
    if(method=="FTCS"):
        plt.plot(times,np.log10(norms), label = r'$\Delta x$' + "={:.3f}".format(dx) + ", " + r'$\Delta t$' + "={:.3f}".format(dt))
    else:
        plt.plot(times,norms, label = r'$\Delta x$' + "={:.3f}".format(dx) + ", " + r'$C_f$' + "={:.3f}".format(c_f))
    plt.legend()

plt.grid()
plt.title("Evolution of the L2 norm")
plt.xlabel("t")
plt.ylabel("L2 norm of u(x,t)")
plt.savefig("norm_{}_repeated.png".format(method))
plt.close()