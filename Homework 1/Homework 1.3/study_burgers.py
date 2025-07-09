#Exercise 3: code which implements the flux and non-flux conservative schemes for a gaussian initial profile for the Burger's equation. 
# The output is a single plot containing the evolution of the profile at some given instants
import numpy as np
import numpy as np
import matplotlib.pyplot as plt

#PARAMETERS OF THE SIMULATION
#simulated domain
L = 10
T = 0.5
#Number of points for space and courant factor
J = 101
c_f=0.5

#The variable "method" selects which method will be applied during the simulation. The available values for this variable are:
#"FLUX","NON_FLUX". For any other value the code will display an error message and the execution
#will be terminated
method = "NON_FLUX"
if(not(method in ["FLUX","NON_FLUX"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()

#Lattice spacings
dx = L/(J-1)
dt = c_f * dx/10.0 #I take 10 as maximum propagation speed as the initial data is a gaussian with 10 as max height

frames_number = int(T/dt)
save_freq = int(0.1/dt)

#flux function
def flux(u):
    return 0.5*np.multiply(u,u)

#Each method results in a different scheme. Here we set the function that implements the scheme for each of the possible choices
if(method == "FLUX"):
    def advance(space_arr):
        return space_arr - (dt/dx)*(flux(space_arr) - flux(np.roll(space_arr,+1)))
elif(method == "NON_FLUX"):
    def advance(space_arr):
        return space_arr - (dt/dx)*np.multiply(space_arr,space_arr-np.roll(space_arr,+1))

x = np.linspace(0,10,J)
x0 = 5
u_init = 10*np.exp(-(x-x0)**2)
u_old = u_init.copy()
u_curr = []

time = 0.0
times = []
frame=0

plt.figure(1)

#Simulation loop
#NOTE: sometimes the condition may be changed to "frame<=frames_number" if the simulation time does not reach T=20
while(frame<frames_number):
    frame+=1
    u_curr = advance(u_old)
    u_old = u_curr.copy()
    time+=dt
    print("timestep #{}, time: {}".format(int(time/dt), time))
    times.append(time)
    if(frame%save_freq==0):
        plt.plot(x,u_curr,label = "t={:.3f}".format(time))
    
plt.ylim(-1,11)
plt.grid()
plt.plot(x,u_init,label = "t=0")
plt.title("Evolution of the profile, " + r'$\Delta x $=' + "{:.3f}".format(dx) + ", " + r'$C_f $=' + "{}".format(c_f))
plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.legend()
#AUTOMATIC NAMING CONVENTION
#use J and courant factor in the title. The courant factor is multiplied by 100 and then the result is converted to an integer: 
#this is just to have no "." in the file name, which may corrupt the file (By doing this the name encodes the courant factor up to 
#the second decimal place)
plt.savefig("final_Upwind_{}_J{}_cf{}.png".format(method,J,int(c_f*100)))
plt.close()
