# Write a short simulation tool, that computes the trajectories of gravitational systems with addaptive Stepsize (2D)

# import packages
import numpy as np
import matplotlib.pyplot as plt

# Constants
M_sun = 2*10**30
G = 6.67*10**-11

# initial conditions
x_i = 4e12
y_i = 0

v_xi = 0
v_yi = 500

N = 100
t_i = 0
t_f = 60 * 60 * 24* 365 *100
h = (t_f- t_i)/N


# define underlying forces/ functions
def f(r,t):
    x, y, v_x, v_y = r
    
    fx = v_x
    fy = v_y

    fv_x = - G*M_sun *x / (x**2 + y**2)**1.5
    fv_y = - G*M_sun *y/ (x**2 + y**2)**1.5

    return np.array([fx, fy, fv_x, fv_y] ,float)



# RK4 without adapt. stepsize
def RK4(t_i, t_f ,h, x_i, y_i, v_xi, v_yi):
    tpoints = np.arange(t_i, t_f ,h)
    xpoints = []
    ypoints = []

    r = np.array([x_i, y_i, v_xi, v_yi] ,float)
    for t in tpoints:
        xpoints.append(r[0])
        ypoints.append(r[1])
        k1 = h*f(r,t)
        k2 = h*f(r+0.5*k1,t+0.5*h)
        k3 = h*f(r+0.5*k2,t+0.5*h)
        k4 = h*f(r+k3,t+h)
        r += (k1+2*k2+2*k3+k4)/6

    return xpoints, ypoints
x, y = RK4(t_i, t_f ,h, x_i, y_i, v_xi, v_yi)




# RK4 WITH adapt. stepsize
def RK4_adapt(t_i, t_f ,h, x_i, y_i, v_xi, v_yi):
    xpoints = []
    ypoints = []

    r = np.array([x_i, y_i, v_xi, v_yi] ,float)
    xpoints.append(r[0])
    ypoints.append(r[1])
    time=[]
    time.append(0)
    t= t_i
    delta = 1000 /(60*60*24*365)      # one kilometer per year

    while t <= t_f:

        r_temporary = r.copy()
        # One step, big h
        h_b =2*h
        k1 = h_b*f(r,t)
        k2 = h_b*f(r+0.5*k1,t+0.5*h_b)
        k3 = h_b*f(r+0.5*k2,t+0.5*h_b)
        k4 = h_b*f(r+k3,t+h_b)
        r_b  = r+ (k1+2*k2+2*k3+k4)/6


        # Two steps, smaller h
        x_temporary = []
        y_temporary = []
        for _ in (1,2):
            k1 = h*f(r,t)
            k2 = h*f(r+0.5*k1,t+0.5*h)
            k3 = h*f(r+0.5*k2,t+0.5*h)
            k4 = h*f(r+k3,t+h)
            r += (k1+2*k2+2*k3+k4)/6
            x_temporary.append(r[0])
            y_temporary.append(r[1])
            


        x1 = r[0]   # estimate from smaller step
        y1 = r[1]   # estimate from smaller step
        x2 = r_b[0] # estimate from big step
        y2 = r_b[1] # estimate from big step

        epsilon_x   = 1/30* (x1 -x2)
        epsilon_y   = 1/30* (y1 -y2)
        epsilon     = np.sqrt(epsilon_x**2 + epsilon_y**2)
        rho = (30 * h * delta)/(epsilon)


        if rho >= 1:     # better then threshold
            xpoints.append(x_temporary[0])  # after first small step
            xpoints.append(x_temporary[1])  # after second small step

            ypoints.append(y_temporary[0])
            ypoints.append(y_temporary[1])
            time.append(time[-1]+h/2)
            time.append(time[-1]+h/2)
            h_new = h * rho**0.25
            if h_new >= 2*h:
                h_new = 2*h
            t = t+2*h
            h= h_new
        


        elif rho < 1:
            r= r_temporary.copy()
            h = h * rho**0.25

    print(time)
    return xpoints, ypoints, time


x_adap, y_adap, t= RK4_adapt(t_i, t_f ,h, x_i, y_i, v_xi, v_yi)
#x, y = RK4(t_i, t_f ,h, x_i, y_i, v_xi, v_yi)

print(len(x_adap))
print(len(t))

#plt.scatter(x, y , marker='o')
plt.plot(x_adap, y_adap, linewidth= 1)
plt.plot(x_adap, y_adap , marker='x')
#plt.plot(t, x_adap)
plt.plot(2, 3, marker='o', color='yellow', markersize=10)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
#plt.axis("equal")
plt.show()