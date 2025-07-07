#Exercise 1: This program produces an animation of the evolution of a given initial profile, using one of four possible schemes
#for the advection equation and a gaussian initial profile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#PARAMETERS OF THE SIMULATION
#velocity
a = 1.0
#simulated domain
L = 10
T = 20
#Number of points for space (J) and time (N)
J = 101
N = 101
c_f=0.5

method = "LEAPFROG"
if(not(method in ["FTCS","LAX_FRIEDRICHS","LEAPFROG","LAX_WENDROFF"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()
#Lattice spacings
dx = L/(J-1)
if(method == "FTCS"):
    dt = T/(N-1)
else:
    dt = c_f*dx/a
#Define frames number (used to stop the simulation)
frames_number = int(20/dt)

if(method == "FTCS"):
    def advance(space_arr):
        return space_arr - (a*dt/(2*dx))*(np.roll(space_arr, -1) - np.roll(space_arr, +1))
elif(method == "LAX_FRIEDRICHS"):
    def advance(space_arr):
        return 0.5*(np.roll(space_arr,-1) + np.roll(space_arr,+1)) - (a*dt/(2*dx))*(np.roll(space_arr, -1) - np.roll(space_arr, +1))
elif(method == "LEAPFROG"):
    def advance_LF(space_arr):
        return 0.5*(np.roll(space_arr,-1) + np.roll(space_arr,+1)) - (a*dt/(2*dx))*(np.roll(space_arr, -1) - np.roll(space_arr, +1))

    def advance_Leapfrog(space_arr,space_arr_older):
        return space_arr_older - (a*dt/dx)*(np.roll(space_arr, -1) - np.roll(space_arr, +1))
elif(method == "LAX_WENDROFF"):
    def advance(space_arr):
        return space_arr - 0.5*(a*dt/dx)*(np.roll(space_arr,-1) - np.roll(space_arr,+1)) + 0.5*((a*dt/dx)**2)*(np.roll(space_arr,-1) - 2*space_arr + np.roll(space_arr,+1))

x = np.linspace(0,10,J)
#Initial gaussian profile
x0 = 5
u_init = np.exp(-(x-x0)**2)
#arrays used to store information about the profile at times t^(n+1) and t^n
u_old = u_init.copy()
u_curr = []

#create figure with a single subplot
fig, ax = plt.subplots()
plot, = ax.plot(x,u_init)
plt.grid()
ax.set_xlabel("x")
ax.set_ylabel("u(x,t)")
plt.title(method + ", " + r'$C_f = $' + "{}".format(c_f))

time = 0.0
#arrays to save different time instants and the corresponding values of the norm, to plot the norm evolution after the simulation
times = []
norms = []

#For animation purposes: intialize plot
def init():
    plot.set_ydata([])
    return plot,

#Function defining how to update the plot
def update(frame):
    global u_curr,u_old,time

    if(method == "LEAPFROG"):
        if(frame == 0):
            print("first step\n")
            u_curr = advance_LF(u_old)
            plot.set_ydata(u_curr)
            norms.append(np.log10(np.linalg.norm(u_curr)))
        else:
            u_new = advance_Leapfrog(u_curr,u_old)
            u_old = u_curr.copy()
            u_curr = u_new.copy()
            plot.set_ydata(u_new)
            norms.append(np.log10(np.linalg.norm(u_new)))
    else:
        #apply the scheme
        u_curr = advance(u_old)
        u_old = u_curr.copy()
        norms.append(np.linalg.norm(u_curr))

    plot.set_ydata(u_curr)    

    time += dt
    print("frame #{}/{}\ttime: {}".format(frame+1,frames_number,time))
    times.append(time)
    
    plt.title(method + ", " + r'$C_f = $' + "{}".format(c_f),  + ", t={:.3f}".format(time))
    return plot,

ani = animation.FuncAnimation(fig=fig, func=update, interval=1,init_func=init, frames=frames_number, repeat = False)
#OPTION TO SAVE THE ANIMATION TO A FILE
ani.save(filename="{}.gif".format(method), writer="imagemagick")
# #OPTION TO SHOW THE ANIMATION LIVE
# plt.show()
