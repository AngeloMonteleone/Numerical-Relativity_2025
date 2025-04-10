#Exercise 1, FCTS method for advection equation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#Define parameters of the simulation
a = 1.0
L = 10
T = 20
J = 101
N = 101

#Lattice spacings
dx = L/(J-1)
dt = T/(N-1)

def advance(space_arr):
    # ret = []
    # for j in range(0,len(space_arr)):
    #     ret.append(space_arr[j%J]-(a*dt/(2*dx))*(space_arr[(j+1)%J]-space_arr[(j-1)%J]))
    ret = space_arr - (a*dt/(2*dx))*(np.roll(space_arr, +1) - np.roll(space_arr, -1))

    return ret

x = np.linspace(0,10,J)
x0 = 0.5
u_init = np.exp(-(x-x0)**2)
u_old = u_init.copy()
u_curr = []

#create figure with a single subplot
fig, ax = plt.subplots()
plot, = ax.plot(x,u_init)
ax.set_xlabel("x")
ax.set_ylabel("u(x,t)")

# for n in range(0,N):
#     u_curr = advance(u_init)
#     u_init = u_curr.copy()

time = 0.0
times = []
norms = []
def update(frame):
    global u_curr,u_old,time

    if(time > 20.0):
        ani.event_source.stop()
        #close previous plot
        plt.close()
        #plot initial and final state
        plt.plot(x,u_init,label = "t=0")
        plt.plot(x,u_curr,label = "t=20")
        plt.title("Comparison between t=0 and t=20")
        plt.xlabel("x")
        plt.ylabel("u(x,t)")
        plt.ylim(0,1.5)
        plt.legend()
        plt.savefig("final_FCTS.png")
        #close the plot (otherwise animation crashes)
        plt.close()

        plt.plot(times,norms)
        plt.title("Evolution of the L2 norm")
        plt.xlabel("t")
        plt.ylabel("L2 norm of u(x,t)")
        plt.savefig("norm_FCTS.png")
        plt.close()

    u_curr = advance(u_old)
    u_old = u_curr.copy()

    plot.set_ydata(u_curr)

    time += dt
    print("frame #{}/{}\ttime: {}".format(frame+1,N,time))
    times.append(time)
    norms.append(np.linalg.norm(u_curr))
    return plot,

ani = animation.FuncAnimation(fig=fig, func=update, interval=10)
plt.show()