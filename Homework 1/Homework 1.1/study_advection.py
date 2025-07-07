#Exercise 1: code which implements FTCS, Lax-Friedrichs, Lax-wendroff, Leapfrog methods. 
# The outputs are two plots
#- One contains four subplots that represent the evolution of the profile at some given instants
#- The other plot shows the evolution of the norm of the solution over time
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
c_f=0.5#Relevant for Lax-Friedrichs, Leapfrog, Lax-Wendroff
col = "C0"

#The variable "method" selects which method will be applied during the simulation. The available values for this variable are:
#"FTCS","LAX_FRIEDRICHS","LAX_WENDROFF","LEAPFROG". For any other value the code will display an error message and the execution
#will be terminated
method = "LEAPFROG"
if(not(method in ["FTCS","LAX_FRIEDRICHS","LEAPFROG","LAX_WENDROFF"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()

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

#Define the moments in which to save the evolved profile. The "save_times" array contains three values, which are integers smaller then
# the total amount of frames. When the frame counter reaches one of these values a plot is drawn.
if(method == "FTCS"):
    save_times = [int(frames_number/10),int(frames_number*2/10),int(frames_number*3/10)]
else:
    save_freq = int(5/dt)
    save_times = [save_freq,2*save_freq,3*save_freq]

x = np.linspace(0,10,J)
#Initial gaussian profile
x0 = 5
u_init = np.exp(-(x-x0)**2)
#arrays used to store information about the profile at times t^(n+1) and t^n
u_old = u_init.copy()
u_curr = []

time = 0.0
frame = 0
#arrays to save different time instants and the corresponding values of the norm, to plot the norm evolution after the simulation
times = []
norms = []

#Create a figure with 4 subplots when T=20
plt.figure(1)
fig, ax = plt.subplots(1,int(T/5), figsize=(24,6))
plotnum = 0

#Simulation loop
#NOTE: sometimes the condition may be changed to "frame<=frames_number" if the simulation time does not reach T=20
while(frame<frames_number):
    if(method == "LEAPFROG"):
        #In the case of the Leapfrog method we need to define two separate steps: For the first frame we apply the Lax-Friedrichs scheme
        #while for all the other frames we apply the ordinary Leapfrog method
        if(frame == 0):
            print("first step\n")
            u_curr = advance_LF(u_old)
            norms.append(np.log10(np.linalg.norm(u_curr)))
        else:
            u_new = advance_Leapfrog(u_curr,u_old)
            u_old = u_curr.copy()
            u_curr = u_new.copy()
            norms.append(np.log10(np.linalg.norm(u_new)))
    else:
        #For all the schemes different form Leapfrog we just apply the standard "advance" function, which implements the whole method
        u_curr = advance(u_old)
        u_old = u_curr.copy()
        norms.append(np.linalg.norm(u_curr))
    frame+=1
    time+=dt
    print("timestep #{}, time: {}".format(frame, time))
    times.append(time)
    if(frame in save_times):
        #Make a plot at a intermaediate time
        ax[plotnum].plot(x,u_curr,label = "t={:.3f}".format(time), color = col)
        ax[plotnum].set_ylim(-0.2,2)
        ax[plotnum].grid()
        ax[plotnum].plot(x,u_init,label = "t=0",color = "black", linestyle = "dashed")#Initial profile for Comparison
        ax[plotnum].set_xlabel("x")
        ax[plotnum].set_ylabel("u(x,t)")
        ax[plotnum].legend(fontsize = 10)
        plotnum+=1

#Plot the final profile
ax[plotnum].plot(x,u_curr,label = "t={:.3f}".format(time), color = col)
ax[plotnum].set_ylim(-0.2,2)
ax[plotnum].grid()
ax[plotnum].plot(x,u_init,label = "t=0",color = "black", linestyle = "dashed")
ax[plotnum].set_xlabel("x")
ax[plotnum].set_ylabel("u(x,t)")
ax[plotnum].legend(fontsize = 10)
#AUTOMATIC NAMING CONVENTION
if(method == "FTCS"):
    #use dx and dt in the title
    fig.suptitle("Evolution of the profile, " + r'$\Delta x $=' + "{:.3f}".format(dx) + ", " + r'$\Delta t $=' + "{:.3f}".format(dt), fontsize = 20)
    #use number of points J and N in file names
    plt.savefig("Comparison_{}_J{}_N{}.png".format(method,J,N))
else:
    #use dx and courant factor in the title
    fig.suptitle("Evolution of the profile, " + r'$\Delta x $=' + "{:.3f}".format(dx) + ", " + r'$C_f $=' + "{}".format(c_f), fontsize = 20)
    #use J and courant factor in the title. The courant factor is multiplied by 100 and then the result is converted to an integer: 
    #this is just to have no "." in the file name, which may corrupt the file (By doing this the name encodes the courant factor up to 
    #the second decimal place)
    plt.savefig("Comparison_{}_J{}_cf{}.png".format(method,J,int(c_f*100)))   

plt.close()

#Plot the evolution of the norm
plt.grid()
plt.xlabel("t")
plt.ylabel("L2 norm of u(x,t)")
if(method == "FTCS"):
    plt.plot(times,np.log10(norms))
    plt.title("Evolution of the L2 norm, " + r'$\Delta x $=' + "{:.3f}".format(dx) + ", " + r'$\Delta t $=' + "{:.3f}".format(dt))
    plt.savefig("Norm_{}_J{}.png".format(method,J))
else:
    plt.plot(times,norms)
    plt.title("Evolution of the L2 norm, " + r'$\Delta x $=' + "{:.3f}".format(dx) + ", " + r'$C_f $=' + "{}".format(c_f))
    plt.savefig("Norm_{}_J{}_cf{}.png".format(method,J,int(c_f*100)))
plt.close()