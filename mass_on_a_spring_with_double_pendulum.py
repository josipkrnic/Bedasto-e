import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
def rk4(N=40000,t0=0,tN=400,s0=0,v0=0,theta10=0,w10=0,theta20=0,w20=0,l0=20,l1=4,l2=4,y0=10,M=20,m1=1,m2=1,g=9.81,k=100,fi=(numpy.pi/6)):
    def funkcija(f,s,v,theta1,w1,theta2,w2,t):
        return f(s,v,theta1,w1,theta2,w2,t)
    def a(s,v,theta1,w1,theta2,w2,t):
        det1 = (M+m1+m2)*(m1+m2)+2*(m1+m2)*numpy.cos(fi+theta1)*numpy.cos(fi+theta2)*numpy.cos(theta1-theta2)-(m1+m2)*m2*(numpy.cos(fi+theta1)**2+numpy.cos(fi+theta2)**2)-(m1+m2)**2*numpy.cos(theta1-theta2)**2
        W1 = (m1+m2)*l1*w1**2*((m1+m2)*numpy.sin(fi+theta1)-m2*numpy.sin(fi+theta1)*numpy.cos(theta1-theta2)**2+m2*numpy.cos(fi+theta1)*numpy.cos(theta1-theta2)*numpy.sin(theta1-theta2)-m2*numpy.sin(theta1-theta2)*numpy.cos(fi+theta2))
        W2 = m2*l2*w2**2*((m1+m2)*numpy.sin(fi+theta2)-m2*numpy.sin(fi+theta2)*numpy.cos(theta1-theta2)**2-m2*numpy.cos(fi+theta2)*numpy.cos(theta1-theta2)*numpy.sin(theta1-theta2)-(m1+m2)*numpy.sin(theta1-theta2)*numpy.cos(fi+theta1))
        E = (-(M+m1+m2)*g*numpy.sin(fi)+k*s)*(m2*numpy.cos(theta1-theta2)**2-m1-m2)
        T1 = (m1+m2)*g*numpy.sin(theta1)*(-(m1+m2)*numpy.cos(fi+theta1)-m2*numpy.cos(fi+theta2)*numpy.cos(theta1-theta2))
        T2 = (m1+m2)*m2*g*numpy.sin(theta2)*(numpy.cos(fi+theta2)-numpy.cos(fi+theta1)*numpy.cos(theta1-theta2))
        return (W1+W2+E+T1+T2)/det1
    def alfa1(s,v,theta1,w1,theta2,w2,t):
        det2 = l1*((M+m1+m2)*(m1+m2)+2*(m1+m2)*numpy.cos(fi+theta1)*numpy.cos(fi+theta2)*numpy.cos(theta1-theta2)-(m1+m2)*m2*(numpy.cos(fi+theta1)**2+numpy.cos(fi+theta2)**2)-(m1+m2)**2*numpy.cos(theta1-theta2)**2)
        W1 = w1**2*l1*((m1+m2)*m2*numpy.sin(fi+theta1)*numpy.cos(fi+theta2)*numpy.cos(theta1-theta2)-(m1+m2)**2*numpy.cos(fi+theta2)*numpy.sin(fi+theta2)+(m1+m2)*m2*numpy.cos(fi+theta1)*numpy.cos(fi+theta2)*numpy.sin(theta1-theta2)-(M+m1+m2)*m2*numpy.cos(theta1-theta2)*numpy.sin(theta1-theta2))
        W2 = m2*l2*w2**2*(m2*numpy.sin(fi+theta1)*numpy.cos(fi+theta2)*numpy.sin(theta1-theta2)-(m1+m2)*numpy.cos(fi+theta1)*numpy.sin(fi+theta2)-(M+m1+m2)*numpy.sin(theta1-theta2)+m2*numpy.cos(fi+theta2)**2*numpy.sin(theta1-theta2))
        Eg = (m1+m2+M)*g*numpy.sin(fi)*((m1+m2)*numpy.cos(fi+theta1)-m2*numpy.cos(fi+theta2)*numpy.cos(theta1-theta2))
        Ek = k*s*((m1+m2)*numpy.cos(fi+theta1)-m2*numpy.cos(fi+theta2)*numpy.cos(theta1-theta2))
        T1 = (m1+m2)*g*numpy.sin(theta1)*(numpy.cos(fi+theta2)**2*m2-(M+m1+m2))
        T2 = m2*g*numpy.sin(theta2)*((M+m1+m2)*numpy.cos(theta1-theta2)-(m1+m2)*numpy.cos(fi+theta1)*numpy.cos(fi+theta2))
        return (W1+W2+Eg+Ek+T1+T2)/det2
    def alfa2(s,v,theta1,w1,theta2,w2,t):
        det3 = l2*((M+m1+m2)*(m1+m2)+2*(m1+m2)*numpy.cos(fi+theta1)*numpy.cos(fi+theta2)*numpy.cos(theta1-theta2)-(m1+m2)*m2*(numpy.cos(fi+theta1)**2+numpy.cos(fi+theta2)**2)-(m1+m2)**2*numpy.cos(theta1-theta2)**2)
        W1 = (m1+m2)*l1*w1**2*((m1+m2)*numpy.cos(fi+theta2)*numpy.sin(fi+theta1)*numpy.cos(theta1-theta2)-(m1+m2)*numpy.cos(fi+theta2)*numpy.sin(fi+theta1)+(M+m1+m2)*numpy.sin(theta1-theta2)-(m1+m2)*numpy.cos(fi+theta1)**2*numpy.sin(theta1-theta2))
        W2 = m2*l2*w2**2*((m1+m2)*numpy.cos(fi+theta1)*numpy.sin(fi+theta2)*numpy.cos(theta1-theta2)-(m1+m2)*numpy.cos(fi+theta2)*numpy.sin(fi+theta2)-(M+m1+m2)*numpy.cos(theta1-theta2)*numpy.sin(theta1-theta2)+(m1+m2)*numpy.cos(fi+theta1)*numpy.cos(fi+theta2)*numpy.sin(theta1-theta2))
        Eg = (M+m1+m2)*(m1+m2)*g*numpy.sin(fi)*(numpy.cos(fi+theta2)-numpy.cos(fi+theta1)*numpy.cos(theta1-theta2))
        Ek = (m1+m2)*k*s*(numpy.cos(fi+theta2)-numpy.cos(fi+theta1)*numpy.cos(theta1-theta2))
        T1 = (m1+m2)*g*numpy.sin(theta1)*((m1+m2)*numpy.cos(fi+theta1)*numpy.cos(fi+theta2)-(M+m1+m2)*numpy.cos(theta1-theta2))
        T2 = (m1+m2)*g*numpy.sin(theta2)*((m1+m2)*numpy.cos(fi+theta1)**2-(M+m1+m2))
        return (W1+W2+Eg+Ek+T1+T2)/det3
    l0 = 20
    h = (tN-t0)/N
    s = [s0]
    v = [v0]
    theta1 = [theta10]
    w1 = [w10]
    theta2 = [theta20]
    w2 = [w20]
    v = [v0]
    t = [t0]
    for i in range(1,N):
        t.append(t[i-1]+h)
        kv1_1 = a(s[i-1],v[i-1],theta1[i-1],w1[i-1],theta2[i-1],w2[i-1],t[i-1])
        kv1_2 = alfa1(s[i-1],v[i-1],theta1[i-1],w1[i-1],theta2[i-1],w2[i-1],t[i-1])
        kv1_3 = alfa2(s[i-1],v[i-1],theta1[i-1],w1[i-1],theta2[i-1],w2[i-1],t[i-1])
        ky1_1 = v[i-1]
        ky1_2 = w1[i-1]
        ky1_3 = w2[i-1]
        kv2_1 = a(s[i-1]+ky1_1*h/2,v[i-1]+kv1_1*h/2,theta1[i-1]+ky1_2*h/2,w1[i-1]+kv1_2*h/2,theta2[i-1]+ky1_3*h/2,w2[i-1]+kv1_3*h/2,t[i-1]+h/2)
        kv2_2 = alfa1(s[i-1]+ky1_1*h/2,v[i-1]+kv1_1*h/2,theta1[i-1]+ky1_2*h/2,w1[i-1]+kv1_2*h/2,theta2[i-1]+ky1_3*h/2,w2[i-1]+kv1_3*h/2,t[i-1]+h/2)
        kv2_3 = alfa2(s[i-1]+ky1_1*h/2,v[i-1]+kv1_1*h/2,theta1[i-1]+ky1_2*h/2,w1[i-1]+kv1_2*h/2,theta2[i-1]+ky1_3*h/2,w2[i-1]+kv1_3*h/2,t[i-1]+h/2)
        ky2_1 = (v[i-1]+kv1_1*h/2)
        ky2_2 = (w1[i-1]+kv1_2*h/2)
        ky2_3 = (w2[i-1]+kv1_3*h/2)
        kv3_1 = a(s[i-1]+ky2_1*h/2,v[i-1]+kv2_1*h/2,theta1[i-1]+ky2_2*h/2,w1[i-1]+kv2_2*h/2,theta2[i-1]+ky2_3*h/2,w2[i-1]+kv2_3*h/2,t[i-1]+h/2)
        kv3_2 = alfa1(s[i-1]+ky2_1*h/2,v[i-1]+kv2_1*h/2,theta1[i-1]+ky2_2*h/2,w1[i-1]+kv2_2*h/2,theta2[i-1]+ky2_3*h/2,w2[i-1]+kv2_3*h/2,t[i-1]+h/2)
        kv3_3 = alfa2(s[i-1]+ky2_1*h/2,v[i-1]+kv2_1*h/2,theta1[i-1]+ky2_2*h/2,w1[i-1]+kv2_2*h/2,theta2[i-1]+ky2_3*h/2,w2[i-1]+kv2_3*h/2,t[i-1]+h/2)
        ky3_1 = (v[i-1]+kv2_1*h/2)
        ky3_2 = (w1[i-1]+kv2_2*h/2)
        ky3_3 = (w2[i-1]+kv2_3*h/2)
        kv4_1 = a(s[i-1]+ky3_1*h,v[i-1]+kv3_1*h,theta1[i-1]+ky3_2*h,w1[i-1]+kv3_2*h,theta2[i-1]+ky3_3*h,w2[i-1]+kv3_3*h,t[i])
        kv4_2 = alfa1(s[i-1]+ky3_1*h,v[i-1]+kv3_1*h,theta1[i-1]+ky3_2*h,w1[i-1]+kv3_2*h,theta2[i-1]+ky3_3*h,w2[i-1]+kv3_3*h,t[i])
        kv4_3 = alfa2(s[i-1]+ky3_1*h,v[i-1]+kv3_1*h,theta1[i-1]+ky3_2*h,w1[i-1]+kv3_2*h,theta2[i-1]+ky3_3*h,w2[i-1]+kv3_3*h,t[i])
        ky4_1 = (v[i-1]+kv3_1*h)
        ky4_2 = (w1[i-1]+kv3_2*h)
        ky4_3 = (w2[i-1]+kv3_3*h)
        v.append(v[i-1]+h*(kv1_1+2*kv2_1+2*kv3_1+kv4_1)/6)
        w1.append(w1[i-1]+h*(kv1_2+2*kv2_2+2*kv3_2+kv4_2)/6)
        w2.append(w2[i-1]+h*(kv1_3+2*kv2_3+2*kv3_3+kv4_3)/6)
        s.append(s[i-1]+h*(ky1_1+2*ky2_1+2*ky3_1+ky4_1)/6)
        theta1.append(theta1[i-1]+h*(ky1_2+2*ky2_2+2*ky3_2+ky4_2)/6)
        theta2.append(theta2[i-1]+h*(ky1_3+2*ky2_3+2*ky3_3+ky4_3)/6)
    return t,s,theta1,v,w1,theta2,w2
t,s,theta1,v,w1,theta2,w2 = rk4()
min2pi = -2*numpy.pi
M = len(s)
theta_graf = []
novi = 0
nov = 0
for i in range(len(theta1)):
    f = theta1[i]
    if f >= 0 and f <= 2*numpy.pi:
        theta_graf.append(f)
    elif f > 2*numpy.pi:
        novi = f - int(f/(2*numpy.pi))*2*numpy.pi
        theta_graf.append(novi)
    elif f < 0:
        nov = f + abs(int(f/(2*numpy.pi)-1)*2*numpy.pi)
        theta_graf.append(nov)
theta_graf_finale = []
for j in theta_graf:
    if j >= 0 and j <= numpy.pi/2:
        theta_graf_finale.append(j)
    elif j > numpy.pi/2 and j <= 3*numpy.pi/2:
        theta_graf_finale.append(numpy.pi-j)
    elif j > 3*numpy.pi/2 and j <= 2*numpy.pi:
        theta_graf_finale.append(j-2*numpy.pi)
#crtanje grafova
plt.plot(t,s)
plt.title("s-t graf")
plt.xlabel("t [s]")
plt.ylabel("produljenost opruge - s [m]")
plt.show()
plt.plot(t,theta_graf)
plt.title("theta-t graf")
plt.xlabel("t [s]")
plt.ylabel("kut otklona njihala - theta [rad]")
plt.show()
plt.plot(t,theta_graf_finale)
plt.title("theta-t graf, theta=kut otklona od vertikale")
plt.xlabel("t [s]")
plt.ylabel("theta [rad]")
plt.show()
plt.plot(t,v)
plt.title("v-t graf")
plt.xlabel("t [s]")
plt.ylabel("v [m/s]")
plt.show()
plt.plot(t,w1)
plt.title("w-t graf")
plt.xlabel("t [s]")
plt.ylabel("w [rad/s]")
plt.show()
x1 = []
y1 = []
x2 = []
y2 = []
x3 = []
y3 = []
l1 = 4
l2 = 4
y0 = 15
l0 = 15
fi = numpy.pi/6
for i in s:
    x1.append((l0+i)*numpy.cos(fi))
    y1.append(y0-(l0+i)*numpy.sin(fi))
for i in range(len(s)):
    x2.append(x1[i] + l1*numpy.sin(theta1[i]))
    y2.append(y1[i] - l1*numpy.cos(theta1[i]))
for k in range(len(theta2)):
    x3.append(x2[k] + l2*numpy.sin(theta2[k]))
    y3.append(y2[k] - l2*numpy.cos(theta2[k]))
x1_1 = []
x1_2 = []
x1_3 = []
x1_4 = []
y1_1 = []
y1_2 = []
y1_3 = []
y1_4 = []
b = 1.7 #bočna stranica pravokutnika
alfa = numpy.pi/6 #kut tijela na kosini
hip = b/numpy.cos(numpy.pi/6)
a = b/numpy.tan(fi)
for i in range(len(x1)):
    x1_4.append(x1[i]+hip*numpy.cos(alfa+fi))
    y1_4.append(y1[i]-hip*numpy.sin(alfa+fi))
    x1_2.append(x1[i]-hip*numpy.cos(alfa+fi))
    y1_2.append(y1[i]+hip*numpy.sin(alfa+fi))
for j in range(len(x1)):
    x1_3.append(x1_4[j]-a*numpy.cos(fi))
    y1_3.append(y1_4[j]+a*numpy.sin(fi))
    x1_1.append(x1_4[j]+b*numpy.sin(fi))
    y1_1.append(y1_4[j]+b*numpy.cos(fi))
x0_1 = x1[0]
y0_1 = y1[0]
x0_2 = x2[0]
y0_2 = y2[0]
x0_3 = x3[0]
y0_3 = y3[0]
x01_1 = x1_1[0]
y01_1 = y1_1[0]
x01_2 = x1_2[0]
y01_2 = y1_2[0]
x01_3 = x1_3[0]
y01_3 = y1_3[0]
x01_4 = x1_4[0]
y01_4 = y1_4[0]
fig = plt.figure("PROJEKT",figsize=(10,10))
line1, = plt.plot([],[],"o-",c="b",lw=2)
line2, = plt.plot([],[],"o-",c="y",lw=2)
line3, = plt.plot([],[],"o-",c="y",lw=2)
line4, = plt.plot([],[],"o-",c="y",lw=2)
line5, = plt.plot([],[],"o-",c="y",lw=2)
line6, = plt.plot([],[],"o-",c="r",lw=2)
line78, = plt.plot([],[],"o-",c="b",lw=2)
anim1, = plt.plot([x0_1],[y0_1],"yo")
anim2, = plt.plot([x0_2],[y0_2],"bo")
anim78, = plt.plot([x0_3],[y0_3],"bo")
anim3, = plt.plot([x01_1],[y01_1],"yo")
anim4, = plt.plot([x01_2],[y01_2],"yo")
anim5, = plt.plot([x01_3],[y01_3],"yo")
anim6, = plt.plot([x01_4],[y01_4],"yo")
def init_line1():
    line1.set_data([],[])
    return line1,
def init_line2():
    line2.set_data([],[])
    return line2,
def init_line3():
    line3.set_data([],[])
    return line3,
def init_line4():
    line4.set_data([],[])
    return line4,
def init_line5():
    line5.set_data([],[])
    return line5,
def init_line6():
    line6.set_data([],[])
    return line6,
def init_line78():
    line78.set_data([],[])
    return line78,
def init_anim1():
    anim1.set_data([],[])
    return anim1,
def init_anim2():
    anim2.set_data([],[])
    return anim2,
def init_anim3():
    anim3.set_data([],[])
    return anim3,
def init_anim4():
    anim4.set_data([],[])
    return anim4,
def init_anim5():
    anim5.set_data([],[])
    return anim5,
def init_anim6():
    anim6.set_data([],[])
    return anim6,
def init_anim78():
    anim78.set_data([],[])
    return anim78,
def animate_line1(i):
    linijax = [x1[i],x2[i]]
    linijay = [y1[i],y2[i]]
    line1.set_data(linijax,linijay)
    return line1,
def animate_line2(i):
    linijax = [x1_1[i],x1_2[i]]
    linijay = [y1_1[i],y1_2[i]]
    line2.set_data(linijax,linijay)
    return line2,
def animate_line3(i):
    linijax = [x1_2[i],x1_3[i]]
    linijay = [y1_2[i],y1_3[i]]
    line3.set_data(linijax,linijay)
    return line3,
def animate_line4(i):
    linijax = [x1_3[i],x1_4[i]]
    linijay = [y1_3[i],y1_4[i]]
    line4.set_data(linijax,linijay)
    return line4,
def animate_line5(i):
    linijax = [x1_4[i],x1_1[i]]
    linijay = [y1_4[i],y1_1[i]]
    line5.set_data(linijax,linijay)
    return line5,
def animate_line6(i):
    linijax = [0,x1[i]]
    linijay = [y0+0.75,y1[i]]
    line6.set_data(linijax,linijay)
    return line6,
def animate_line78(i):
    linijax = [x2[i],x3[i]]
    linijay = [y2[i],y3[i]]
    line78.set_data(linijax,linijay)
    return line78,
def animate_anim1(i):
    iks1 = x1[i]
    ips1 = y1[i]
    anim1.set_data(iks1,ips1)
    return anim1,
def animate_anim2(i):
    iks2 = x2[i]
    ips2 = y2[i]
    anim2.set_data(iks2,ips2)
    return anim2,
def animate_anim78(i):
    iks78 = x3[i]
    ips78 = y3[i]
    anim2.set_data(iks78,ips78)
    return anim78,
def animate_anim3(i):
    iks3 = x1_1[i]
    ips3 = y1_1[i]
    anim3.set_data(iks3,ips3)
    return anim3,
def animate_anim4(i):
    iks4 = x1_2[i]
    ips4 = y1_2[i]
    anim4.set_data(iks4,ips4)
    return anim4,
def animate_anim5(i):
    iks5 = x1_3[i]
    ips5 = y1_3[i]
    anim5.set_data(iks5,ips5)
    return anim5,
def animate_anim6(i):
    iks6 = x1_4[i]
    ips6 = y1_4[i]
    anim6.set_data(iks6,ips6)
    return anim6,
myAnimation0 = animation.FuncAnimation(fig,animate_anim1,frames=M,init_func=init_anim1,interval=0,blit=False)
myAnimation1 = animation.FuncAnimation(fig,animate_anim2,frames=M,init_func=init_anim2,interval=0,blit=False)
myAnimation2 = animation.FuncAnimation(fig,animate_line1,frames=M,init_func=init_line1,interval=0,blit=False)
myAnimation3 = animation.FuncAnimation(fig,animate_anim3,frames=M,init_func=init_anim3,interval=0,blit=False)
myAnimation4 = animation.FuncAnimation(fig,animate_anim4,frames=M,init_func=init_anim4,interval=0,blit=False)
myAnimation5 = animation.FuncAnimation(fig,animate_anim5,frames=M,init_func=init_anim5,interval=0,blit=False)
myAnimation6 = animation.FuncAnimation(fig,animate_anim6,frames=M,init_func=init_anim6,interval=0,blit=False)
myAnimation7 = animation.FuncAnimation(fig,animate_line2,frames=M,init_func=init_line2,interval=0,blit=False)
myAnimation8 = animation.FuncAnimation(fig,animate_line3,frames=M,init_func=init_line3,interval=0,blit=False)
myAnimation9 = animation.FuncAnimation(fig,animate_line4,frames=M,init_func=init_line4,interval=0,blit=False)
myAnimation10 = animation.FuncAnimation(fig,animate_line5,frames=M,init_func=init_line5,interval=0,blit=False)
myAnimation11 = animation.FuncAnimation(fig,animate_line6,frames=M,init_func=init_line6,interval=0,blit=False)
myAnimation12 = animation.FuncAnimation(fig,animate_line78,frames=M,init_func=init_line78,interval=0,blit=False)
myAnimation13 = animation.FuncAnimation(fig,animate_anim78,frames=M,init_func=init_anim78,interval=0,blit=False)
plt.grid()
plt.xlim(-5,27)
plt.ylim(-5,27)
nasupr_x = []
nasupr_y = []
prilez_x = []
prilez_y = []
hipot_x = []
hipot_y = []
y0 = 14
x00 = y0/numpy.tan(fi)
range_x = int(x00*100)
for i in range(y0*100+150):
    nasupr_x.append(0)
    nasupr_y.append(i*0.01)
for j in range(range_x):
    prilez_x.append(j*0.01)
    prilez_y.append(0)
range_hip = int(numpy.sqrt(x00**2+y0**2))*100
for k in range(range_hip):
    hipot_x.append(numpy.cos(fi)*k*0.01)
    hipot_y.append(numpy.sin(fi)*k*0.01)
hipot_y.reverse()
plt.plot(nasupr_x,nasupr_y,"k")
plt.plot(prilez_x,prilez_y,"k")
plt.plot(hipot_x,hipot_y,"k")
plt.show()
