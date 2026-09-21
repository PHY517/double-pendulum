import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters
g = 9.81
L1 = 0.20
L2 = 0.15
m1 = 1.0
m2 = 1.0

# Double-pendulum equations
def rhs(variables, t):

    theta1, omega1, theta2, omega2 = variables
    delta = theta1 - theta2

    den1 = L1 * (2*m1 + m2 - m2*np.cos(2*delta))
    den2 = L2 * (2*m1 + m2 - m2*np.cos(2*delta))

    domega1 = (
        -g*(2*m1 + m2)*np.sin(theta1)
        -m2*g*np.sin(theta1 - 2*theta2)
        -2*m2*np.sin(delta) *
        (omega2**2*L2 + omega1**2*L1*np.cos(delta))
    ) / den1

    domega2 = (
        2*np.sin(delta) *
        (omega1**2*L1*(m1 + m2)
         + g*(m1 + m2)*np.cos(theta1)
         + omega2**2*L2*m2*np.cos(delta))
    ) / den2

    return omega1, domega1, omega2, domega2


# Time array
t_init = 0.0
t_final = 15.0
dt = 0.01

N = int((t_final - t_init)/dt) + 1
t = np.linspace(t_init, t_final, N)

# Initial conditions
theta1_0 = 181
theta2_0 = 120
omega1_0 = 0.0
omega2_0 = 0.0

variables_init = [
    np.radians(theta1_0),
    omega1_0,
    np.radians(theta2_0),
    omega2_0
]

# Solve equations
solution = integrate.odeint(rhs, variables_init, t)

theta1 = solution[:, 0]
theta2 = solution[:, 2]
omega1 = solution[:, 1]
omega2 = solution[:, 3]


def energy(theta1, omega1, theta2, omega2):
    T = (
        0.5*(m1 + m2)*L1**2*omega1**2
        + 0.5*m2*L2**2*omega2**2
        + m2*L1*L2*omega1*omega2*np.cos(theta1 - theta2)
    )

    V = (
        -(m1 + m2)*g*L1*np.cos(theta1)
        - m2*g*L2*np.cos(theta2)
    )

    return T + V



E = energy(theta1, omega1, theta2, omega2)

E_avg = np.mean(E)
E_std = np.std(E)

print(f"Average energy = {E_avg:.8f} J")
print(f"Energy dispersion = {E_std:.8e} J")
print(f"Relative dispersion = {E_std/abs(E_avg):.8e}")

E0 = E[0]

drift = np.max(np.abs(E - E0))
relative_drift = drift / abs(E0)

print(f"Maximum energy drift = {drift:.8e} J")
print(f"Maximum relative drift = {relative_drift:.8e}")





# Set up animation
fig, ax = plt.subplots()

ax.set_aspect("equal")
ax.set_xlim(-(L1 + L2), L1 + L2)
ax.set_ylim(-(L1 + L2), L1 + L2)

pivot, = ax.plot(0, 0, "ks", ms=8)
rods, = ax.plot([], [], "k-")
bob1, = ax.plot([], [], "bo", ms=10)
bob2, = ax.plot([], [], "ro", ms=10)


def animate(i):

    x1 = L1*np.sin(theta1[i])
    y1 = -L1*np.cos(theta1[i])

    x2 = x1 + L2*np.sin(theta2[i])
    y2 = y1 - L2*np.cos(theta2[i])

    rods.set_data([0, x1, x2], [0, y1, y2])
    bob1.set_data([x1], [y1])
    bob2.set_data([x2], [y2])

    return rods, bob1, bob2, pivot


mt = 1000*(t[1] - t[0])

ani = FuncAnimation(
    fig,
    animate,
    frames=len(t),
    interval=mt
)

plt.show()



