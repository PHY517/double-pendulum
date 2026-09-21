# Double Pendulum, Chaos, and Lyapunov Exponents

## Overview

The double pendulum consists of two pendulums connected end to end. Despite its simple construction and deterministic equations of motion, it can exhibit highly complex and chaotic behavior.

This project uses Python to:

- Solve the equations of motion of a double pendulum.
- Animate its dynamics.
- Monitor conservation of mechanical energy.
- Compare two systems with nearly identical initial conditions.
- Identify chaotic behavior from the divergence of nearby trajectories.
- Estimate the largest Lyapunov exponent from a semilogarithmic plot.

## Project Files

### [`double_pendulum.py`](double_pendulum.py)

This is the main double-pendulum simulation. It:

- Defines the coupled equations of motion.
- Solves them using `scipy.integrate.odeint`.
- Calculates the total mechanical energy at every time.
- Reports the average energy, energy dispersion, and maximum energy drift.
- Animates the motion of a single double pendulum.

Use this file to study the basic dynamics and verify that the numerical solution approximately conserves mechanical energy.

### [`double_pendulum_chaos.py`](double_pendulum_chaos.py)

This script demonstrates sensitive dependence on initial conditions. It:

- Simulates two double pendulums whose initial angles differ by only a very small amount.
- Displays the two animations next to each other.
- Calculates the distance `d` between the second bobs.
- Shows `d(t)` on a semilogarithmic plot.
- Moves a marker along the separation curve while the pendulums are animated.

Use this file to visualize how two initially similar trajectories diverge and to identify the exponential-growth region associated with a positive Lyapunov exponent.

## Physical Model

The system contains two point masses, $m_1$ and $m_2$, connected by massless rods of lengths $L_1$ and $L_2$.

The angular coordinates are:

- $\theta_1$: angle of the first rod measured from the downward vertical.
- $\theta_2$: angle of the second rod measured from the downward vertical.

The corresponding angular velocities are:

- $\omega_1 = d\theta_1/dt$
- $\omega_2 = d\theta_2/dt$

The Cartesian positions of the masses are

```math
x_1 = L_1\sin\theta_1,
\qquad
y_1 = -L_1\cos\theta_1,
```

and

```math
x_2 = x_1 + L_2\sin\theta_2,
\qquad
y_2 = y_1 - L_2\cos\theta_2.
```

The coupled nonlinear equations of motion are integrated numerically using `scipy.integrate.odeint`.

## Mechanical Energy

The total mechanical energy is the sum of the kinetic and potential energies:

```math
E = T + V.
```

For the double pendulum, the kinetic energy is

```math
T = \frac{1}{2}(m_1+m_2)L_1^2\omega_1^2
  + \frac{1}{2}m_2L_2^2\omega_2^2
  + m_2L_1L_2\omega_1\omega_2\cos(\theta_1-\theta_2),
```

and the potential energy is

```math
V = -(m_1+m_2)gL_1\cos\theta_1
    -m_2gL_2\cos\theta_2.
```

A useful Python function is:

```python
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
```

The mean energy and its standard deviation can be used to monitor numerical conservation:

```python
E = energy(theta1, omega1, theta2, omega2)

E_average = np.mean(E)
E_dispersion = np.std(E)

print("Average energy:", E_average)
print("Energy dispersion:", E_dispersion)
```

For an ideal isolated system, the total energy is constant. In a numerical solution, small changes may appear because of integration error. A smaller energy dispersion generally indicates better numerical conservation.

## What Is Chaos?

Chaos is complex behavior produced by deterministic equations. The future motion is completely determined by the initial conditions, but tiny uncertainties in those initial conditions can grow rapidly.

For example, consider two double pendulums whose initial angles differ by only

```text
180.0000 degrees
180.0001 degrees
```

Their motions initially appear identical. After sufficient time, however, their trajectories may become very different.

This phenomenon is called **sensitive dependence on initial conditions**.

Chaos does not mean that the equations contain randomness. It means that long-term prediction becomes difficult because no physical or numerical initial condition can be specified with infinite precision.

## Measuring Trajectory Separation

To test for sensitive dependence on initial conditions, solve the equations for two nearly identical initial states.

Let $(x_{2a},y_{2a})$ and $(x_{2b},y_{2b})$ be the positions of the second bob in the two simulations. Their separation is

```math
d(t)=\sqrt{[x_{2a}(t)-x_{2b}(t)]^2+[y_{2a}(t)-y_{2b}(t)]^2}.
```

In Python:

```python
d = np.sqrt((x2a - x2b)**2 + (y2a - y2b)**2)
```

A rapid increase in $d(t)$ shows that the two trajectories are diverging.

## Lyapunov Exponent

During an interval of exponential divergence, the separation approximately follows

```math
d(t) \approx d_0 e^{\lambda t},
```

where:

- $d_0$ is the initial separation.
- $\lambda$ is the largest Lyapunov exponent.

Taking the natural logarithm gives

```math
\ln d(t) = \ln d_0 + \lambda t.
```

Therefore, a plot of $\ln d$ versus time should contain an approximately straight region. The slope of that region estimates the Lyapunov exponent:

```math
\lambda \approx \frac{d\ln d}{dt}.
```

Equivalently, exponential divergence appears approximately linear when $d(t)$ is shown on a semilogarithmic plot:

```python
plt.semilogy(t, d)
plt.xlabel("Time (s)")
plt.ylabel("Separation d (m)")
plt.show()
```

A positive Lyapunov exponent indicates that nearby trajectories separate exponentially and is an important signature of chaos.

## Estimating the Lyapunov Exponent

Select an interval in which $\ln d(t)$ is approximately linear, before the separation saturates:

```python
from scipy.stats import linregress

log_d = np.log(d)
mask = (t > 2.0) & (t < 6.0)

fit = linregress(t[mask], log_d[mask])
lyapunov = fit.slope

print("Estimated Lyapunov exponent:", lyapunov)
```

The selected fitting interval depends on the trajectory. It should exclude:

1. The initial transient region.
2. Points affected strongly by numerical roundoff.
3. The late-time saturation region.

This two-trajectory fit provides a simple classroom estimate of the largest Lyapunov exponent. A more rigorous calculation repeatedly renormalizes the perturbation so that the two trajectories remain close in phase space.

## Signatures of Chaos

The main signatures to examine are:

1. **Sensitive dependence on initial conditions**  
   Nearly identical starting conditions eventually produce visibly different trajectories.

2. **Exponential separation**  
   The distance between nearby trajectories grows approximately as $e^{\lambda t}$ during an intermediate time interval.

3. **Linear region on a semilog plot**  
   Exponential growth appears approximately linear when $d(t)$ is plotted with a logarithmic vertical axis.

4. **Positive Lyapunov exponent**  
   A positive fitted slope indicates exponential divergence.

5. **Aperiodic motion**  
   The pendulum motion does not repeat with a simple, fixed period.

Irregular motion alone is not sufficient evidence of chaos. The divergence of nearby trajectories provides a clearer quantitative test.

## Important Numerical Considerations

- Use radians in the equations of motion.
- Use a sufficiently small time step.
- Monitor total mechanical energy.
- Verify that the observed divergence is not caused primarily by poor numerical integration.
- Compare results using smaller time steps or tighter integration tolerances.
- Avoid fitting the Lyapunov exponent after the separation has saturated.
- Do not evaluate `np.log(d)` at points where `d` is zero.

One safe approach is

```python
mask = (d > 0) & (t > 2.0) & (t < 6.0)
fit = linregress(t[mask], np.log(d[mask]))
```

## Requirements

```bash
pip install numpy scipy matplotlib
```

## Summary

The double pendulum is a useful example of classical chaos because it combines simple physical components with nonlinear, strongly coupled motion. By comparing two nearly identical simulations, one can observe sensitive dependence on initial conditions. If their separation grows approximately exponentially, the slope of $\ln d(t)$ provides an estimate of the largest Lyapunov exponent.

> A chaotic system is deterministic, but tiny uncertainties in its initial state can grow exponentially and make long-term prediction impractical.
