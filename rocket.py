import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

## constants
deg = np.pi/180
g0 = 9.81                   # sea-level acceleration of gravity (m/s)
Re = 6378E3                 # ...Radius of the earth (m)
h_scale = 7.5e3             # ...Density scale height (m)
rho0 = 1.225                # ...Sea level density of atmosphere (kg/m^3)
    
## inputs
diam = 5.0                  # ...Vehicle diameter (m)
cd = 0.5                    # ...Drag coefficient (assumed constant)
m0 = 68000                  # ...Lift-off mass (kg)
thrust = 933910             # ...Thrust (kN)
isp = 390                   # ...Specific impulse (s)
h_turn  = 130               # ...Height at which pitchover begins (m)
t_bo = 260                  # ...burn out (s)
    
area = np.pi/4*np.power(diam,2) # ...Frontal area (m^2)
m_dot = thrust/isp/g0       # ...Propellant mass flow rate (kg/s)
mprop = m0 - m_dot*t_bo     # ...Propellant mass (kg)
    
# ...Initial conditions:
v0 = 0                      # ...Initial velocity (m/s)
gamma0 = 89.85*deg          # ...Initial flight path angle (rad)
x0 = 0                      # ...Initial downrange distance (km)
h0 = 0                      # ...Initial altitude (km)
vD0	= 0;                    # ...Initial value of velocity loss to drag (m/s)
vG0 = 0;                    # ...Initial value of velocity loss due to gravity (m/s) 
y0 = [v0, gamma0, x0, h0, vD0, vG0]
dydt = [None]*len(y0)
t = np.linspace(0,t_bo)

def rates(t, y):
    v, gamma, x, h, vD, vG = y
    v *= 1E3
    x *= 1E3
    h *= 1E3
    vD *= 1E3
    vG *= 1E3

    if t < t_bo:
        m = m0 - m_dot*t    # ...Current vehicle mass
        T = thrust
    else:
        m = m0 - m_dot*t_bo# ...Current vehicle mass
        T = 0          # ...Current thrust

    g = g0/(1 + h/Re)**2
    rho = rho0 * np.exp(-h/h_scale)
    D = 1/2 * rho * v**2 * area * cd

    if h <= h_turn:
        gamma_dot = 0
        v_dot     = T/m - D/m - g
        x_dot     = 0
        h_dot     = v
        vG_dot    = -g
    else:
        v_dot = T/m - D/m - g*np.sin(gamma) # ...Equation 11.6
        gamma_dot = -1/v * (g - v**2/(Re + h))*np.cos(gamma)# ...Equation 11.7
        x_dot    = Re/(Re + h)* v *np.cos(gamma)   # ...Equation 11.8(1)
        h_dot    = v*np.sin(gamma)       # ...Equation 11.8(2)
        vG_dot	 = -g*np.sin(gamma)      # ...Gravity loss rate

    vD_dot = -D/m   # ...Drag loss rate

    dydt[0]  = v_dot*1E-3
    dydt[1]  = gamma_dot
    dydt[2]  = x_dot*1E-3
    dydt[3]  = h_dot*1E-3
    dydt[4]  = vD_dot*1E-3
    dydt[5]  = vG_dot*1E-3
    return dydt

sol = integrate.solve_ivp(rates, [min(t),max(t)], y0)
print(sol)
v, gamma, x, h, vD, vG = sol.y
    
    
print('Initial flight path angle = {} [deg] '.format(gamma0/deg))
print('Pitchover altitude= {} g [m]'.format(h_turn))
print('Burn time = {} [s]'.format(t_bo))
print('Final speed = {} [km/s]'.format(v[-1]))
print('Final flight path angle = {} [deg] '.format(gamma[-1]/deg))
print('Altitude = {} [km]  '.format(h[-1]))
print('Downrange distance = {} [km] '.format(x[-1]))
print('Drag loss = {} [km/s]'.format(-vD[-1]))
print('Gravity loss = {} [km/s]'.format(-vG[-1]))

plt.plot(x,h)
plt.ylabel('altitude')
plt.xlabel('downrange distance')
plt.show()
