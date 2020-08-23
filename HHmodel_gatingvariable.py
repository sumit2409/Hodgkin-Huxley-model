import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
tmin = 0.0
tmax = 50
gK = 36.0
gNa = 120.0
gL = 0.3
Cm = 1.0
VK = -12.0
VNa = 115.0
Vl = 10.613
T = np.linspace(tmin,tmax,10000)
def alpha_n(Vm):
    return (0.01 * (10.0 - Vm)) / (np.exp(1.0 - (0.1 * Vm)) - 1.0)
def beta_n(Vm):
    return 0.125 * np.exp(-Vm / 80.0)
def alpha_m(Vm):
    return (0.1 * (25.0 - Vm)) / (np.exp(2.5 - (0.1 * Vm)) - 1.0)
def beta_m(Vm):
    return 4.0 * np.exp(-Vm / 18.0)
def alpha_h(Vm):
    return 0.07 * np.exp(-Vm / 20.0)
def beta_h(Vm):
    return 1.0 / (np.exp(3.0 - (0.1 * Vm)) + 1.0)
def n_inf(Vm=0.0):
    return alpha_n(Vm) / (alpha_n(Vm) + beta_n(Vm))
def m_inf(Vm=0.0):
    return alpha_m(Vm) / (alpha_m(Vm) + beta_m(Vm))
def h_inf(Vm=0.0):
    return alpha_h(Vm) / (alpha_h(Vm) + beta_h(Vm))

def Id(t):
    if 0.0 < t < 10:
        return 0
    elif 10.0 < t < 30:
        return 50
    return 0
def compute_derivatives(y, t):
    dy = np.zeros((4,))
    
    Vm = y[0]
    n = y[1]
    m = y[2]
    h = y[3]
    GK = (gK / Cm) * np.power(n, 4.0)
    GNa = (gNa / Cm) * np.power(m, 3.0) * h
    GL = gL / Cm
    dy[0] = (Id(t)/ Cm) - (GK * (Vm - VK)) - (GNa * (Vm - VNa)) - (GL * (Vm - Vl))
    dy[1] = (alpha_n(Vm) * (1.0 - n)) - (beta_n(Vm) * n)
    dy[2] = (alpha_m(Vm) * (1.0 - m)) - (beta_m(Vm) * m)
    dy[3] = (alpha_h(Vm) * (1.0 - h)) - (beta_h(Vm) * h)
    
    return dy
t=np.linspace(0,60,1000)
Y = np.array([0.0, n_inf(), m_inf(), h_inf()])
Vy = odeint(compute_derivatives, Y, T)
mylist1=[]
mylist2=[]
mylist3=[]
for i in range(0,10000):
    mylist1.append(Vy[i][1])
for i in range(0,10000):
    mylist2.append(Vy[i][2])
for i in range(0,10000):
    mylist3.append(Vy[i][3])
plt.plot(T,mylist1)
plt.plot(T,mylist2)  
plt.plot(T,mylist3)
plt.xlabel('time')
plt.ylabel('x(t)or gating variable')
plt.show()
    
