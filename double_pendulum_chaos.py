import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from matplotlib.animation import FuncAnimation

# --------------------------------------------------
# Physical parameters
# --------------------------------------------------

g = 9.81
L1 = 0.20
L2 = 0.2
m1 = 1.0
m2 = 2.0


# --------------------------------------------------
# Double-pendulum equations
# --------------------------------------------------

def rhs(variables, t):

    theta1, omega1, theta2, omega2 = variables
    delta = theta1 - theta2

    denominator1 = L1 * (
        2*m1 + m2 - m2*np.cos(2*delta)
    )

    denominator2 = L2 * (
        2*m1 + m2 - m2*np.cos(2*delta)
    )

    domega1 = (
        -g*(2*m1 + m2)*np.sin(theta1)
        -m2*g*np.sin(theta1 - 2*theta2)
        -2*m2*np.sin(delta)
        * (
            omega2**2*L2
            + omega1**2*L1*np.cos(delta)
        )
    ) / denominator1

    domega2 = (
        2*np.sin(delta)
        * (
            omega1**2*L1*(m1 + m2)
            + g*(m1 + m2)*np.cos(theta1)
            + omega2**2*L2*m2*np.cos(delta)
        )
    ) / denominator2

    return omega1, domega1, omega2, domega2


# --------------------------------------------------
# Time array
# --------------------------------------------------

t_init = 0.0
t_final = 15.0
dt = 0.01

N = int((t_final - t_init)/dt) + 1
t = np.linspace(t_init, t_final, N)


# --------------------------------------------------
# Nearly identical initial conditions
# --------------------------------------------------

initial1 = [
    np.radians(180.0000),
    0.0,
    np.radians(120.0),
    0.0
]

initial2 = [
    np.radians(180.0001),   # very small difference
    0.0,
    np.radians(120.0),
    0.0
]


# --------------------------------------------------
# Solve both systems
# --------------------------------------------------

solution1 = integrate.odeint(rhs, initial1, t)
solution2 = integrate.odeint(rhs, initial2, t)

theta1a = solution1[:, 0]
theta2a = solution1[:, 2]

theta1b = solution2[:, 0]
theta2b = solution2[:, 2]


# --------------------------------------------------
# Calculate the bob positions
# --------------------------------------------------

x1a = L1*np.sin(theta1a)
y1a = -L1*np.cos(theta1a)

x2a = x1a + L2*np.sin(theta2a)
y2a = y1a - L2*np.cos(theta2a)

x1b = L1*np.sin(theta1b)
y1b = -L1*np.cos(theta1b)

x2b = x1b + L2*np.sin(theta2b)
y2b = y1b - L2*np.cos(theta2b)


# Separation between the second bobs
d = np.sqrt((x2a - x2b)**2 + (y2a - y2b)**2)


# --------------------------------------------------
# Create plots
# --------------------------------------------------

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5))

limit = 1.1*(L1 + L2)

for ax in (ax1, ax2):
    ax.set_aspect("equal")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.plot(0, 0, "ks", ms=8)

ax1.set_title("Pendulum A")
ax2.set_title("Pendulum B")

rod_a, = ax1.plot([], [], "k-")
bob1_a, = ax1.plot([], [], "bo", ms=8)
bob2_a, = ax1.plot([], [], "ro", ms=10)

rod_b, = ax2.plot([], [], "k-")
bob1_b, = ax2.plot([], [], "bo", ms=8)
bob2_b, = ax2.plot([], [], "ro", ms=10)


# Semilog plot of the separation
ax3.semilogy(t, d)
marker, = ax3.semilogy([], [], "ro")

ax3.set_xlabel("Time (s)")
ax3.set_ylabel("Separation d (m)")
ax3.set_title("Divergence of trajectories")
ax3.grid()


# --------------------------------------------------
# Animation function
# --------------------------------------------------

def animate(i):

    # Update pendulum A
    rod_a.set_data(
        [0, x1a[i], x2a[i]],
        [0, y1a[i], y2a[i]]
    )

    bob1_a.set_data([x1a[i]], [y1a[i]])
    bob2_a.set_data([x2a[i]], [y2a[i]])

    # Update pendulum B
    rod_b.set_data(
        [0, x1b[i], x2b[i]],
        [0, y1b[i], y2b[i]]
    )

    bob1_b.set_data([x1b[i]], [y1b[i]])
    bob2_b.set_data([x2b[i]], [y2b[i]])

    # Move the marker along the separation curve
    marker.set_data([t[i]], [d[i]])

    return (
        rod_a, bob1_a, bob2_a,
        rod_b, bob1_b, bob2_b,
        marker
    )


# Animation interval in milliseconds
mt = 1000*(t[1] - t[0])

ani = FuncAnimation(
    fig,
    animate,
    frames=len(t),
    interval=mt
)

plt.tight_layout()
plt.show()


