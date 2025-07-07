#Exercise 2: This program produces an animation of the evolution of a given initial profile, using one of two possible schemes
#for the advection equation and the stepfunction as an initial profile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#PARAMETERS OF THE SIMULATION
#velocity
a = 1.0
#simulated domain
L = 10
T = 20
#Number of points for space and courant factor
J = 101
c_f=0.5

method = "LAX_WENDROFF"
if(not(method in ["LAX_FRIEDRICHS","LAX_WENDROFF"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()
#Lattice spacings
dx = L/(J-1)
dt = c_f*dx/a
#Define frames number (used to stop the simulation)
frames_number = int(20/dt)

if(method == "LAX_FRIEDRICHS"):
    def advance(space_arr):
        return 0.5*(np.roll(space_arr,-1) + np.roll(space_arr,+1)) - (a*dt/(2*dx))*(np.roll(space_arr, -1) - np.roll(space_arr, +1))
elif(method == "LAX_WENDROFF"):
    def advance(space_arr):
        return space_arr - 0.5*(a*dt/dx)*(np.roll(space_arr,-1) - np.roll(space_arr,+1)) + 0.5*((a*dt/dx)**2)*(np.roll(space_arr,-1) - 2*space_arr + np.roll(space_arr,+1))

x = np.linspace(0,10,J)
#Function defining the step function. np.where acts as an "elementwise if": the first argument is an array of boolean values
#(representing the regions where the step-function returns 1 or 0) and the output is an array of the same size as the original one 
# where if the boolean value is true (i.e.: the condition is satisfied) for a given index, then the output array will have in that 
# index the value of the second argument (1), otherwise the value of the third (0)
def step_function(arg):
    return np.where((arg>4) & (arg<6),1,0)

u_init = step_function(x)
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
