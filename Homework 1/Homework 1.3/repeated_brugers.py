#Exercise 3: This program implements one of two possible schemes for the Burgers' equation, repeating the simulation for 
# different values of the spacing dx. The output is a plot of the evolution of the 1-norm in which each line represents 
# the evolution of the norm for a different value of dx
import numpy as np
import matplotlib.pyplot as plt

#PARAMETERS OF THE SIMULATION
#simulated domain
L = 10
T = 0.5
#Courant factor
c_f = 0.5

#The variable "method" selects which method will be applied during the simulation. The available values for this variable are:
#"FLUX","NON_FLUX". For any other value the code will display an error message and the execution
#will be terminated
method = "FLUX"
if(not(method in ["FLUX","NON_FLUX"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()

#flux function
def flux(u):
    return 0.5*np.multiply(u,u)

x0 = 5

for J in [101,501,1001]:#VALUES OF THE NUMBER OF POINTS TO USE
    #Each method results in a different scheme. Here we set the function that implements the scheme for each of the possible choices
    if(method == "FLUX"):
        def advance(space_arr):
            return space_arr - (dt/dx)*(flux(space_arr) - flux(np.roll(space_arr,+1)))
    elif(method == "NON_FLUX"):
        def advance(space_arr):
            return space_arr - (dt/dx)*np.multiply(space_arr,space_arr-np.roll(space_arr,+1))

    #Lattice spacing
    dx = L/(J-1)
    dt = c_f*dx/10
    #Define frames number (used to stop the simulation)
    frames_number = int(T/dt)

    x = np.linspace(0,10,J)
    u_init = 10*np.exp(-(x-x0)**2)
    u_old = u_init.copy()
    u_curr = []

    time = 0.0
    frame = 0
    times = []
    norms = []

    while(frame<=frames_number):
        u_curr = advance(u_old)
        u_old = u_curr.copy()
        norms.append(np.linalg.norm(u_curr, ord = 1)/J)
        frame+=1
        time+=dt
        #Here we print also the value of the flux at the boundary
        print("timestep #{}, time: {}, f(0) = {}, f(10) = {}".format(frame, time, flux(u_curr)[0], flux(u_curr)[-1]))
        times.append(time)
    
    plt.plot(times,norms, label = r'$\Delta x$' + "={:.3f}".format(dx) + ", " + r'$C_f$' + "={:.3f}".format(c_f))
    plt.legend()

plt.grid()
plt.title("Evolution of the L1 norm")
plt.xlabel("t")
plt.ylabel("L1 norm of u(x,t)")
plt.savefig("norm_{}_repeated.png".format(method))
plt.close()