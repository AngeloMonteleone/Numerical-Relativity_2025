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
J = 501
N = 101#Relevant only for FTCS
c_f=0.5#Relevant for all methods besides FTCS

method = "LAX_FRIEDRICHS"
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
# plt.title(method + ", " + r'$C_f = $' + "{}".format(c_f))

time = 0.0

#For animation purposes: intialize plot
def init():
    plot.set_ydata(u_init)
    return plot,

#parameter used to speed up the animation. Each frame corresponds to making a numner of steps equal to "repeat_integration"
repeat_integration = 5
#Function defining how to update the plot
def update(frame):
    global u_curr,u_old,time

    if(method == "LEAPFROG"):
        if(frame == 0):
            print("first step\n")
            #apply the scheme
            u_curr = advance_LF(u_old)
            plot.set_ydata(u_curr)
            time += dt
        else:
            for i in range(0,repeat_integration):
                #apply the scheme
                u_new = advance_Leapfrog(u_curr,u_old)
                u_old = u_curr.copy()
                u_curr = u_new.copy()
                time += dt
            plot.set_ydata(u_new)
    else:
        for i in range(0,repeat_integration):
            #apply the scheme
            u_curr = advance(u_old)
            u_old = u_curr.copy()
            time += dt
        plot.set_ydata(u_curr)    

    print("frame #{}\ttime: {:.5f}\tnorm: {:.5f}".format(frame+1,time,np.linalg.norm(u_curr)))
    if(method=="FTCS"):
        plt.title(method + ", " + r'$\Delta x $=' + "{:.3f}".format(dx) + ", " + r'$\Delta t $=' + "{:.3f}".format(dt)  + ", t={:.3f}".format(time))
    else:
        plt.title(method + ", " + r'$\Delta x $=' + "{:.3f}".format(dx) + ", " + r'$C_f $=' + "{}".format(c_f)  + ", t={:.3f}".format(time))
    return plot,

ani = animation.FuncAnimation(fig=fig, func=update, interval=1, init_func=init, frames=int(frames_number/repeat_integration), repeat = False)
#OPTION TO SAVE THE ANIMATION TO A FILE
if(method=="FTCS"):
    ani.save(filename="{}_J{}_N{}.gif".format(method,J,N), writer="imagemagick")
else:
    ani.save(filename="{}_J{}_cf{}.gif".format(method,J,int(c_f*100)), writer="imagemagick")
#OPTION TO SHOW THE ANIMATION LIVE
# plt.show()
