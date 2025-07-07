#Exercise 3, Upwind flux-conservative method for burgers equation
import numpy as np
import matplotlib.pyplot as plt
#PARAMETERS OF THE SIMULATION
#simulated domain
L = 10
T = 0.5
#Number of points for space and courant factor
J = 101
c_f=0.5

method = "FLUX"
if(not(method in ["FLUX","NON_FLUX"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()

#Lattice spacings
dx = L/(J-1)
dt = c_f * dx/10.0 #I take 10 as maximum propagation speed as the initial data is a gaussian with 10 as max height

frames_number = int(0.5/dt)
save_freq = int(0.1/dt)

def flux(u):
    return 0.5*np.multiply(u,u)

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

while(time<=0.5):
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
plt.title("Comparison between different times, " + r'$\Delta x $=' + "{:.3f}".format(dx))
plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.legend()
plt.savefig("final_Upwind_{}_J{}_cf{}.png".format(method,J,int(c_f*100)))