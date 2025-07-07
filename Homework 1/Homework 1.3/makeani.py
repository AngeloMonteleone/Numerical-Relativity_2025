#Exercise 3, Upwind flux-conservative method for burgers equation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#PARAMETERS OF THE SIMULATION
#simulated domain
L = 10
T = 0.5
#Number of points for space and courant factor
J = 101
c_f=0.5

method = "NON_FLUX"
if(not(method in ["FLUX","NON_FLUX"])):
    print("ERROR: method \"{}\" not available!".format(method))
    quit()

#Lattice spacings
dx = L/(J-1)
dt = c_f * dx/10.0 #I take 10 as maximum propagation speed as the initial data is a gaussian with 10 as max height

frames_number = int(0.5/dt)
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

#create figure with a single subplot
fig, ax = plt.subplots()
plot, = ax.plot(x,u_init)
ax.set_xlabel("x")
ax.set_ylabel("u(x,t)")
plt.grid()
plt.ylim(-1,11)
    
time = 0.0
times = []
norms = []

#For animation purposes: intialize plot
def init():
    plot.set_ydata([])
    return plot,

def update(frame):
    global u_curr,u_old,time

    u_curr = advance(u_old)
    u_old = u_curr.copy()
    
    if(method == "FLUX"):
        plt.title("Upwind Flux Conservative " + r'$C_f = $' + "{} ".format(c_f) + r'$\Delta x = $' + "{:.3f}".format(dx) + ", t={:.3f}".format(time))
    else:
        plt.title("Upwind Non-Flux Conservative " + r'$C_f = $' + "{} ".format(c_f) + r'$\Delta x = $' + "{:.3f}".format(dx) + ", t={:.3f}".format(time))
    plot.set_ydata(u_curr)

    time += dt
    print("frame #{}\ttime: {}".format(frame+1,time))
    times.append(time)
    norms.append(np.linalg.norm(u_curr))
    return plot,

ani = animation.FuncAnimation(fig=fig, func=update,init_func=init, frames = frames_number, interval=1, repeat = False)
#OPTION TO SAVE THE ANIMATION TO A FILE
ani.save(filename="Upwind_{}.gif".format(method), writer="imagemagick")
# #OPTION TO SHOW THE ANIMATION LIVE
# plt.show()