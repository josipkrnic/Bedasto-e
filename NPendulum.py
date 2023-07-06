import numpy as np
from math import sin
from math import cos
import matplotlib.pyplot as plt
import matplotlib.animation as animation
def N_Pendulum(N,t_uk,I,M=1,L=1,g=9.81):
    m = M/N
    l = L/N
    def angular_acceleration0(N):
        theta0 = np.zeros(N)
        omega0 = np.zeros(N)
        alfa0 = np.zeros(N)
        for i in range(N):
            theta0[i] = np.random.uniform(-np.pi/4,np.pi/4)
            omega0[i] = np.random.uniform(-5,5)
        A = np.zeros((N,N))
        B = np.zeros(N)
        suma = 0
        for i in range(N):
            for j in range(N):
                if j <= i:
                    suma += -(N-i)*m*l*omega0[j]**2*sin(theta0[i]-theta0[j])
                    A[i][j] = (N-i)*m*l*cos(theta0[i]-theta0[j])
            for k in range(N):
                if k > i:
                    suma += -(N-k)*m*l*omega0[k]**2*sin(theta0[i]-theta0[k])
                    A[i][k] = (N-k)*m*l*cos(theta0[i]-theta0[k])
            B[i]=(-(N-i)*m*g*sin(theta0[i])+suma)
        alfa0 = np.linalg.inv(A).dot(B)
        return theta0,omega0,alfa0
    def angular_acceleration(N,theta0,omega0,M=1,L=1,g=9.81):
        alfa = np.zeros(N)
        A = np.zeros((N,N))
        B = np.zeros(N)
        for i in range(N):
            suma = 0
            for j in range(N):
                if j <= i:
                    suma += -(N-i)*m*l*omega0[j]**2*sin(theta0[i]-theta0[j])
                    A[i][j] = (N-i)*m*l*cos(theta0[i]-theta0[j])
            for k in range(N):
                if k > i:
                    suma += -(N-k)*m*l*omega0[k]**2*sin(theta0[i]-theta0[k])
                    A[i][k] = (N-k)*m*l*cos(theta0[i]-theta0[k])
            B[i]=(-(N-i)*m*g*sin(theta0[i])+suma)
        alfa = np.linalg.inv(A).dot(B)
        return alfa
    def N_Pendulum_RK4(N,t_uk,I):
        dt = t_uk/I
        t = []
        theta0,omega0,alfa0 = angular_acceleration0(N)
        theta1 = theta0
        omega1 = omega0
        alfa1 = alfa0
        THETA = np.zeros((I,N))
        OMEGA = np.zeros((I,N))
        for i in range(I):
            t.append(dt*i)
            kv1 = np.zeros(N)
            kv2 = np.zeros(N)
            kv3 = np.zeros(N)
            kv4 = np.zeros(N)
            kx1 = np.zeros(N)
            kx2 = np.zeros(N)
            kx3 = np.zeros(N)
            kx4 = np.zeros(N)
            theta2 = np.zeros(N)
            omega2 = np.zeros(N)
            theta3 = np.zeros(N)
            omega3 = np.zeros(N)
            theta4 = np.zeros(N)
            omega4 = np.zeros(N)
            for j in range(N):
                kv1[j] = alfa1[j]
                kx1[j] = omega1[j]
                theta2[j] = theta1[j]+kx1[j]*dt/2
                omega2[j] = omega1[j]+kv1[j]*dt/2
            alfa2 = angular_acceleration(N,theta2,omega2)
            for j in range(N):
                kv2[j] = alfa2[j]
                kx2[j] = omega2[j]
                theta3[j] = theta2[j]+kx2[j]*dt/2
                omega3[j] = omega2[j]+kv2[j]*dt/2
            alfa3 = angular_acceleration(N,theta3,omega3)
            for j in range(N):
                kv3[j] = alfa3[j]
                kx3[j] = omega3[j]
                theta4[j] = theta3[j]+kx3[j]*dt
                omega4[j] = omega3[j]+kv3[j]*dt
            alfa4 = angular_acceleration(N,theta4,omega4)
            for j in range(N):
                kv4[j] = alfa4[j]
                kx4[j] = omega4[j]
                OMEGA[i][j] = omega1[j] + dt*(kv1[j]+kv2[j]/2+kv3[j]/2+kv4[j])/6
                THETA[i][j] = theta1[j] + dt*(kx1[j]+kx2[j]/2+kx3[j]/2+kx4[j])/6
            theta1 = np.zeros(N)
            omega1 = np.zeros(N)
            theta1 = THETA[i,:]
            omega1 = OMEGA[i,:]
            alfa1 = angular_acceleration(N,theta1,omega1)
        return THETA,OMEGA,t
    T,O,vrijeme = N_Pendulum_RK4(N,t_uk,I)
    return T,O,vrijeme
N = int(input("Unesite broj njihala: "))
T,O,t = N_Pendulum(N,100,10000)
L = 1
l = L/N
if N == 1:
    theta = T[:,0]
    x = []
    y = []
    for i in theta:
        x.append(l*sin(i))
        y.append(-l*cos(i))
    x0 = x[0]
    y0 = y[0]
    fig = plt.figure("N-Pendulum",figsize=(10,10))
    anim, = plt.plot([x0],[y0],"bo")
    line, = plt.plot([],[],"bo-",lw=2)
    def init_line():
        line.set_data([],[])
        return line,
    def init_anim():
        anim.set_data([],[])
        return anim,
    def animate_line(i):
        lx = [0,x[i]]
        ly = [0,y[i]]
        line.set_data(lx,ly)
        return line,
    def animate_anim(i):
        iks = x[i]
        ips = y[i]
        anim.set_data(iks,ips)
        return anim,
    myAnimation0 = animation.FuncAnimation(fig,animate_line,frames=10000,init_func=init_line,interval=0,blit=False)
    myAnimation1 = animation.FuncAnimation(fig,animate_anim,frames=10000,init_func=init_anim,interval=0,blit=False)
    plt.xlim(-L,L)
    plt.ylim(-L-0.1,L-0.5)
    plt.scatter(0,0,c="r")
    plt.show()
if N == 2:
    theta1 = T[:,0]
    theta2 = T[:,1]
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    for i in range(len(theta1)):
        x1.append(l*sin(theta1[i]))
        y1.append(-l*cos(theta1[i]))
        x2.append(x1[i]+l*sin(theta2[i]))
        y2.append(y1[i]-l*cos(theta2[i]))
    x02 = x2[0]
    y02 = y2[0]
    fig = plt.figure("N-Pendulum",figsize=(10,10))
    line1, = plt.plot([],[],"bo-",lw=2)
    anim2, = plt.plot([x02],[y02],"bo")
    line2, = plt.plot([],[],"bo-",lw=2)
    def init_line1():
        line1.set_data([],[])
        return line1,
    def init_line2():
        line2.set_data([],[])
        return line2,
    def init_anim2():
        anim2.set_data([],[])
        return anim2,
    def animate_line1(i):
        lx = [0,x1[i]]
        ly = [0,y1[i]]
        line1.set_data(lx,ly)
        return line1,
    def animate_line2(i):
        lx = [x1[i],x2[i]]
        ly = [y1[i],y2[i]]
        line2.set_data(lx,ly)
        return line2,
    def animate_anim2(i):
        iks = x2[i]
        ips = y2[i]
        anim2.set_data(iks,ips)
        return anim2,
    myAnimation0 = animation.FuncAnimation(fig,animate_line1,frames=10000,init_func=init_line1,interval=0,blit=False)
    myAnimation2 = animation.FuncAnimation(fig,animate_line2,frames=10000,init_func=init_line2,interval=0,blit=False)
    myAnimation3 = animation.FuncAnimation(fig,animate_anim2,frames=10000,init_func=init_anim2,interval=0,blit=False)
    plt.xlim(-L,L)
    plt.ylim(-L-0.1,L-0.5)
    plt.scatter(0,0,c="r")
    plt.show()
if N == 3:
    theta11 = T[:,0]
    theta22 = T[:,1]
    theta33 = T[:,2]
    x11 = []
    y11 = []
    x22 = []
    y22 = []
    x33 = []
    y33 = []
    for i in range(len(theta11)):
        x11.append(l*sin(theta11[i]))
        y11.append(-l*cos(theta11[i]))
        x22.append(x11[i]+l*sin(theta22[i]))
        y22.append(y11[i]-l*cos(theta22[i]))
        x33.append(x22[i]+l*sin(theta33[i]))
        y33.append(y22[i]-l*cos(theta33[i]))
    x033 = x33[0]
    y033 = y33[0]
    fig = plt.figure("N-Pendulum",figsize=(10,10))
    line11, = plt.plot([],[],"bo-",lw=2)
    line22, = plt.plot([],[],"bo-",lw=2)
    anim33, = plt.plot([x033],[y033],"bo")
    line33, = plt.plot([],[],"bo-",lw=2)
    def init_line11():
        line11.set_data([],[])
        return line11,
    def init_line22():
        line22.set_data([],[])
        return line22,
    def init_line33():
        line33.set_data([],[])
        return line33,
    def init_anim33():
        anim33.set_data([],[])
        return anim33,
    def animate_line11(i):
        lx = [0,x11[i]]
        ly = [0,y11[i]]
        line11.set_data(lx,ly)
        return line11,
    def animate_line22(i):
        lx = [x11[i],x22[i]]
        ly = [y11[i],y22[i]]
        line22.set_data(lx,ly)
        return line22,
    def animate_line33(i):
        lx = [x22[i],x33[i]]
        ly = [y22[i],y33[i]]
        line33.set_data(lx,ly)
        return line33,
    def animate_anim33(i):
        iks = x33[i]
        ips = y33[i]
        anim33.set_data(iks,ips)
        return anim33,
    myAnimation0 = animation.FuncAnimation(fig,animate_line11,frames=10000,init_func=init_line11,interval=0,blit=False)
    myAnimation2 = animation.FuncAnimation(fig,animate_line22,frames=10000,init_func=init_line22,interval=0,blit=False)
    myAnimation4 = animation.FuncAnimation(fig,animate_line33,frames=10000,init_func=init_line33,interval=0,blit=False)
    myAnimation5 = animation.FuncAnimation(fig,animate_anim33,frames=10000,init_func=init_anim33,interval=0,blit=False)
    plt.xlim(-L,L)
    plt.ylim(-L-0.1,L-0.5)
    plt.plot(0,0,c="r")
    plt.show()
if N == 4:
    theta11 = T[:,0]
    theta22 = T[:,1]
    theta33 = T[:,2]
    theta44 = T[:,3]
    x11 = []
    y11 = []
    x22 = []
    y22 = []
    x33 = []
    y33 = []
    x44 = []
    y44 = []
    for i in range(len(theta11)):
        x11.append(l*sin(theta11[i]))
        y11.append(-l*cos(theta11[i]))
        x22.append(x11[i]+l*sin(theta22[i]))
        y22.append(y11[i]-l*cos(theta22[i]))
        x33.append(x22[i]+l*sin(theta33[i]))
        y33.append(y22[i]-l*cos(theta33[i]))
        x44.append(x33[i]+l*sin(theta44[i]))
        y44.append(y33[i]-l*cos(theta44[i]))
    x044 = x44[0]
    y044 = y44[0]
    fig = plt.figure("N-Pendulum",figsize=(10,10))
    line11, = plt.plot([],[],"bo-",lw=2)
    line22, = plt.plot([],[],"bo-",lw=2)
    line33, = plt.plot([],[],"bo-",lw=2)
    anim44, = plt.plot([x044],[y044],"bo")
    line44, = plt.plot([],[],"bo-",lw=2)
    def init_line11():
        line11.set_data([],[])
        return line11,
    def init_line22():
        line22.set_data([],[])
        return line22,
    def init_line33():
        line33.set_data([],[])
    def init_line44():
        line44.set_data([],[])
        return line44,
    def init_anim44():
        anim44.set_data([],[])
        return anim44,
    def animate_line11(i):
        lx = [0,x11[i]]
        ly = [0,y11[i]]
        line11.set_data(lx,ly)
        return line11,
    def animate_line22(i):
        lx = [x11[i],x22[i]]
        ly = [y11[i],y22[i]]
        line22.set_data(lx,ly)
        return line22,
    def animate_line33(i):
        lx = [x22[i],x33[i]]
        ly = [y22[i],y33[i]]
        line33.set_data(lx,ly)
        return line33,
    def animate_line44(i):
        lx = [x33[i],x44[i]]
        ly = [y33[i],y44[i]]
        line44.set_data(lx,ly)
        return line44,
    def animate_anim44(i):
        iks = x44[i]
        ips = y44[i]
        anim44.set_data(iks,ips)
        return anim44,
    myAnimation0 = animation.FuncAnimation(fig,animate_line11,frames=10000,init_func=init_line11,interval=0,blit=False)
    myAnimation2 = animation.FuncAnimation(fig,animate_line22,frames=10000,init_func=init_line22,interval=0,blit=False)
    myAnimation4 = animation.FuncAnimation(fig,animate_line33,frames=10000,init_func=init_line33,interval=0,blit=False)
    myAnimation6 = animation.FuncAnimation(fig,animate_line44,frames=10000,init_func=init_line44,interval=0,blit=False)
    myAnimation7 = animation.FuncAnimation(fig,animate_anim44,frames=10000,init_func=init_anim44,interval=0,blit=False)
    plt.xlim(-L,L)
    plt.ylim(-L-0.1,L-0.5)
    plt.scatter(0,0,c="r")
    plt.show()
if N == 5:
    theta11 = T[:,0]
    theta22 = T[:,1]
    theta33 = T[:,2]
    theta44 = T[:,3]
    theta55 = T[:,4]
    x11 = []
    y11 = []
    x22 = []
    y22 = []
    x33 = []
    y33 = []
    x44 = []
    y44 = []
    x55 = []
    y55 = []
    for i in range(len(theta11)):
        x11.append(l*sin(theta11[i]))
        y11.append(-l*cos(theta11[i]))
        x22.append(x11[i]+l*sin(theta22[i]))
        y22.append(y11[i]-l*cos(theta22[i]))
        x33.append(x22[i]+l*sin(theta33[i]))
        y33.append(y22[i]-l*cos(theta33[i]))
        x44.append(x33[i]+l*sin(theta44[i]))
        y44.append(y33[i]-l*cos(theta44[i]))
        x55.append(x44[i]+l*sin(theta55[i]))
        y55.append(y44[i]-l*cos(theta55[i]))
    x055 = x55[0]
    y055 = y55[0]
    fig = plt.figure("N-Pendulum",figsize=(10,10))
    line11, = plt.plot([],[],"bo-",lw=2)
    line22, = plt.plot([],[],"bo-",lw=2)
    line33, = plt.plot([],[],"bo-",lw=2)
    anim55, = plt.plot([x055],[y055],"bo")
    line44, = plt.plot([],[],"bo-",lw=2)
    line55, = plt.plot([],[],"bo-",lw=2)
    def init_line11():
        line11.set_data([],[])
        return line11,
    def init_line22():
        line22.set_data([],[])
        return line22,
    def init_line33():
        line33.set_data([],[])
    def init_line44():
        line44.set_data([],[])
        return line44,
    def init_line55():
        line55.set_data([],[])
        return line55,
    def init_anim55():
        anim55.set_data([],[])
        return anim55,
    def animate_line11(i):
        lx = [0,x11[i]]
        ly = [0,y11[i]]
        line11.set_data(lx,ly)
        return line11,
    def animate_line22(i):
        lx = [x11[i],x22[i]]
        ly = [y11[i],y22[i]]
        line22.set_data(lx,ly)
        return line22,
    def animate_line33(i):
        lx = [x22[i],x33[i]]
        ly = [y22[i],y33[i]]
        line33.set_data(lx,ly)
        return line33,
    def animate_line44(i):
        lx = [x33[i],x44[i]]
        ly = [y33[i],y44[i]]
        line44.set_data(lx,ly)
        return line44,
    def animate_line55(i):
        lx = [x44[i],x55[i]]
        ly = [y44[i],y55[i]]
        line55.set_data(lx,ly)
        return line55,
    def animate_anim55(i):
        iks = x55[i]
        ips = y55[i]
        anim55.set_data(iks,ips)
        return anim55,
    myAnimation0 = animation.FuncAnimation(fig,animate_line11,frames=10000,init_func=init_line11,interval=0,blit=False)
    myAnimation2 = animation.FuncAnimation(fig,animate_line22,frames=10000,init_func=init_line22,interval=0,blit=False)
    myAnimation4 = animation.FuncAnimation(fig,animate_line33,frames=10000,init_func=init_line33,interval=0,blit=False)
    myAnimation6 = animation.FuncAnimation(fig,animate_line44,frames=10000,init_func=init_line44,interval=0,blit=False)
    myAnimation7 = animation.FuncAnimation(fig,animate_anim55,frames=10000,init_func=init_anim55,interval=0,blit=False)
    myAnimation8 = animation.FuncAnimation(fig,animate_line55,frames=10000,init_func=init_line55,interval=0,blit=False)
    plt.xlim(-L,L)
    plt.ylim(-L-0.1,L-0.5)
    plt.scatter(0,0,c="r")
    plt.show()
if N == 6:
    theta11 = T[:,0]
    theta22 = T[:,1]
    theta33 = T[:,2]
    theta44 = T[:,3]
    theta55 = T[:,4]
    theta66 = T[:,5]
    x11 = []
    y11 = []
    x22 = []
    y22 = []
    x33 = []
    y33 = []
    x44 = []
    y44 = []
    x55 = []
    y55 = []
    x66 = []
    y66 = []
    for i in range(len(theta11)):
        x11.append(l*sin(theta11[i]))
        y11.append(-l*cos(theta11[i]))
        x22.append(x11[i]+l*sin(theta22[i]))
        y22.append(y11[i]-l*cos(theta22[i]))
        x33.append(x22[i]+l*sin(theta33[i]))
        y33.append(y22[i]-l*cos(theta33[i]))
        x44.append(x33[i]+l*sin(theta44[i]))
        y44.append(y33[i]-l*cos(theta44[i]))
        x55.append(x44[i]+l*sin(theta55[i]))
        y55.append(y44[i]-l*cos(theta55[i]))
        x66.append(x55[i]+l*sin(theta66[i]))
        y66.append(y55[i]-l*cos(theta66[i]))
    x066 = x66[0]
    y066 = y66[0]
    fig = plt.figure("N-Pendulum",figsize=(10,10))
    line11, = plt.plot([],[],"bo-",lw=2)
    line22, = plt.plot([],[],"bo-",lw=2)
    line33, = plt.plot([],[],"bo-",lw=2)
    anim66, = plt.plot([x066],[y066],"bo")
    line44, = plt.plot([],[],"bo-",lw=2)
    line55, = plt.plot([],[],"bo-",lw=2)
    line66, = plt.plot([],[],"bo-",lw=2)
    def init_line11():
        line11.set_data([],[])
        return line11,
    def init_line22():
        line22.set_data([],[])
        return line22,
    def init_line33():
        line33.set_data([],[])
    def init_line44():
        line44.set_data([],[])
        return line44,
    def init_line55():
        line55.set_data([],[])
        return line55,
    def init_line66():
        line66.set_data([],[])
        return line66,
    def init_anim66():
        anim66.set_data([],[])
        return anim66,
    def animate_line11(i):
        lx = [0,x11[i]]
        ly = [0,y11[i]]
        line11.set_data(lx,ly)
        return line11,
    def animate_line22(i):
        lx = [x11[i],x22[i]]
        ly = [y11[i],y22[i]]
        line22.set_data(lx,ly)
        return line22,
    def animate_line33(i):
        lx = [x22[i],x33[i]]
        ly = [y22[i],y33[i]]
        line33.set_data(lx,ly)
        return line33,
    def animate_line44(i):
        lx = [x33[i],x44[i]]
        ly = [y33[i],y44[i]]
        line44.set_data(lx,ly)
        return line44,
    def animate_line55(i):
        lx = [x44[i],x55[i]]
        ly = [y44[i],y55[i]]
        line55.set_data(lx,ly)
        return line55,
    def animate_line66(i):
        lx = [x55[i],x66[i]]
        ly = [y55[i],y66[i]]
        line66.set_data(lx,ly)
        return line66,
    def animate_anim66(i):
        iks = x66[i]
        ips = y66[i]
        anim66.set_data(iks,ips)
        return anim66,
    myAnimation0 = animation.FuncAnimation(fig,animate_line11,frames=10000,init_func=init_line11,interval=0,blit=False)
    myAnimation2 = animation.FuncAnimation(fig,animate_line22,frames=10000,init_func=init_line22,interval=0,blit=False)
    myAnimation4 = animation.FuncAnimation(fig,animate_line33,frames=10000,init_func=init_line33,interval=0,blit=False)
    myAnimation6 = animation.FuncAnimation(fig,animate_line44,frames=10000,init_func=init_line44,interval=0,blit=False)
    myAnimation7 = animation.FuncAnimation(fig,animate_anim66,frames=10000,init_func=init_anim66,interval=0,blit=False)
    myAnimation8 = animation.FuncAnimation(fig,animate_line55,frames=10000,init_func=init_line55,interval=0,blit=False)
    myAnimation9 = animation.FuncAnimation(fig,animate_line66,frames=10000,init_func=init_line66,interval=0,blit=False)
    plt.xlim(-L,L)
    plt.ylim(-L-0.1,L-0.5)
    plt.scatter(0,0,c="r")
    plt.show()
if N == 7:
    theta11 = T[:,0]
    theta22 = T[:,1]
    theta33 = T[:,2]
    theta44 = T[:,3]
    theta55 = T[:,4]
    theta66 = T[:,5]
    theta77 = T[:,6]
    x11 = []
    y11 = []
    x22 = []
    y22 = []
    x33 = []
    y33 = []
    x44 = []
    y44 = []
    x55 = []
    y55 = []
    x66 = []
    y66 = []
    x77 = []
    y77 = []
    for i in range(len(theta11)):
        x11.append(l*sin(theta11[i]))
        y11.append(-l*cos(theta11[i]))
        x22.append(x11[i]+l*sin(theta22[i]))
        y22.append(y11[i]-l*cos(theta22[i]))
        x33.append(x22[i]+l*sin(theta33[i]))
        y33.append(y22[i]-l*cos(theta33[i]))
        x44.append(x33[i]+l*sin(theta44[i]))
        y44.append(y33[i]-l*cos(theta44[i]))
        x55.append(x44[i]+l*sin(theta55[i]))
        y55.append(y44[i]-l*cos(theta55[i]))
        x66.append(x55[i]+l*sin(theta66[i]))
        y66.append(y55[i]-l*cos(theta66[i]))
        x77.append(x66[i]+l*sin(theta77[i]))
        y77.append(y66[i]-l*cos(theta77[i]))
    x077 = x77[0]
    y077 = y77[0]
    fig = plt.figure("N-Pendulum",figsize=(10,10))
    line11, = plt.plot([],[],"bo-",lw=2)
    line22, = plt.plot([],[],"bo-",lw=2)
    line33, = plt.plot([],[],"bo-",lw=2)
    anim77, = plt.plot([x077],[y077],"bo")
    line44, = plt.plot([],[],"bo-",lw=2)
    line55, = plt.plot([],[],"bo-",lw=2)
    line66, = plt.plot([],[],"bo-",lw=2)
    line77, = plt.plot([],[],"bo-",lw=2)
    def init_line11():
        line11.set_data([],[])
        return line11,
    def init_line22():
        line22.set_data([],[])
        return line22,
    def init_line33():
        line33.set_data([],[])
    def init_line44():
        line44.set_data([],[])
        return line44,
    def init_line55():
        line55.set_data([],[])
        return line55,
    def init_line66():
        line66.set_data([],[])
        return line66,
    def init_line77():
        line77.set_data([],[])
        return line77,
    def init_anim77():
        anim77.set_data([],[])
        return anim77,
    def animate_line11(i):
        lx = [0,x11[i]]
        ly = [0,y11[i]]
        line11.set_data(lx,ly)
        return line11,
    def animate_line22(i):
        lx = [x11[i],x22[i]]
        ly = [y11[i],y22[i]]
        line22.set_data(lx,ly)
        return line22,
    def animate_line33(i):
        lx = [x22[i],x33[i]]
        ly = [y22[i],y33[i]]
        line33.set_data(lx,ly)
        return line33,
    def animate_line44(i):
        lx = [x33[i],x44[i]]
        ly = [y33[i],y44[i]]
        line44.set_data(lx,ly)
        return line44,
    def animate_line55(i):
        lx = [x44[i],x55[i]]
        ly = [y44[i],y55[i]]
        line55.set_data(lx,ly)
        return line55,
    def animate_line66(i):
        lx = [x55[i],x66[i]]
        ly = [y55[i],y66[i]]
        line66.set_data(lx,ly)
        return line66,
    def animate_line77(i):
        lx = [x66[i],x77[i]]
        ly = [y66[i],y77[i]]
        line77.set_data(lx,ly)
        return line77,
    def animate_anim77(i):
        iks = x77[i]
        ips = y77[i]
        anim77.set_data(iks,ips)
        return anim77,
    myAnimation0 = animation.FuncAnimation(fig,animate_line11,frames=10000,init_func=init_line11,interval=0,blit=False)
    myAnimation2 = animation.FuncAnimation(fig,animate_line22,frames=10000,init_func=init_line22,interval=0,blit=False)
    myAnimation4 = animation.FuncAnimation(fig,animate_line33,frames=10000,init_func=init_line33,interval=0,blit=False)
    myAnimation6 = animation.FuncAnimation(fig,animate_line44,frames=10000,init_func=init_line44,interval=0,blit=False)
    myAnimation7 = animation.FuncAnimation(fig,animate_anim77,frames=10000,init_func=init_anim77,interval=0,blit=False)
    myAnimation8 = animation.FuncAnimation(fig,animate_line55,frames=10000,init_func=init_line55,interval=0,blit=False)
    myAnimation9 = animation.FuncAnimation(fig,animate_line66,frames=10000,init_func=init_line66,interval=0,blit=False)
    myAnimation10 = animation.FuncAnimation(fig,animate_line77,frames=10000,init_func=init_line77,interval=0,blit=False)
    plt.xlim(-L,L)
    plt.ylim(-L-0.1,L-0.5)
    plt.scatter(0,0,c="r")
    plt.show()




        
            
            

        
