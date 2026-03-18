# Lyapunov Spectrum of a Hamiltonian Mean-Field System

This repository contains the code developed during my scientific initiation project in computational physics. The main goal is to study dynamical chaos in a Hamiltonian system with long-range interactions by computing its Lyapunov exponents via tangent-space dynamics.

The project deepened my understanding of scientific computing, symplectic integration, and chaos theory.

## Overview

The code simulates the time evolution of $N$ classical particles with positions $(x_i, y_i)$ and momenta $(p_{x,i}, p_{y,i})$, governed by a mean-field Hamiltonian. Using a 4th-order symplectic Suzuki–Trotter integrator, we preserve the symplectic structure of phase space, which is crucial for long-time simulations of Hamiltonian systems. To quantify chaos, we simultaneously evolve a set of tangent vectors in the linearized phase space. The growth rates of these vectors yield the Lyapunov exponents—the hallmark of sensitivity to initial conditions. 

## Mathematical Model

The Hamiltonian is given by:

$$
H = \sum_{i=1}^{N} \frac{p_{x,i}^2 + p_{y,i}^2}{2} + V(x_i, y_i)
$$

where the potential energy contains mean-field couplings:

$$
V = \frac{c}{2}(2 - M_1^2 - M_2^2) + \frac{d}{4}(2 - P_+^2 - P_-^2)
$$

The collective variables are defined as:

$$
M_1 = \frac{1}{N} \sum_{j} e^{i x_j}
$$

$$
M_2 = \frac{1}{N} \sum_{j} e^{i y_j}
$$

$$
P_+ = \frac{1}{N} \sum_{j} e^{i (x_j + y_j)}
$$

$$
P_- = \frac{1}{N} \sum_{j} e^{i (x_j - y_j)}
$$

The parameters $c$ and $d$ control the coupling strengths (set to $c_1 = c_2 = 1$, $c = -1$, $d = 1$ in the code).

Periodic boundary conditions confine the positions to the interval:

$$
[-2n\pi, \, 2n\pi], \quad \text{with } n = 4
$$

## Numerical Methods

### Symplectic Integration (Suzuki–Trotter 4th Order)

Hamiltonian dynamics are integrated using a composition of position and momentum updates (drift–kick steps).

The 4th-order symplectic scheme ensures good energy conservation over long times, which is essential for reliable chaos indicators.

The coefficients are defined according to the Suzuki fractal decomposition:

$$
s = \frac{1}{4 - 4^{1/3}}
$$

$$
\delta_1 = h s, \quad
\delta_2 = \frac{h s}{2}, \quad
\delta_3 = \frac{h(1 - 3s)}{2}, \quad
\delta_4 = h(1 - 4s)
$$

where $h = 0.5$ is the time step.

---

### Tangent Space Dynamics & Lyapunov Exponents

The evolution of an infinitesimal perturbation $(\delta x, \delta y, \delta p_x, \delta p_y)$ is governed by the linearized equations of motion (i.e., the Jacobian of the flow). These equations are integrated alongside the main trajectory using the same symplectic scheme.

At each time step, we:

- Evolve the tangent vector  
- Accumulate the logarithm of its norm:
  
  $$
  \lambda(t) = \frac{1}{t} \log \frac{|\mathbf{w}(t)|}{|\mathbf{w}(0)|}
  $$

- Renormalize the vector to prevent overflow/underflow  

The final value after long integration approximates the **largest Lyapunov exponent**.

- A **positive exponent** indicates chaotic motion  
- A **zero exponent** indicates regular (quasiperiodic) motion  

---

## Code Structure

| Function / Section        | Description |
|--------------------------|------------|
| `evolveQ`, `evolveP`     | Drift and kick steps for the original system |
| `evolvedQ`, `evolvedP`   | Linearized versions for tangent dynamics |
| `Suzuki4`, `Suzuki4_linear` | 4th-order integrators for original and tangent systems |
| `waterbag`               | Initializes particles in a rectangular region of phase space |
| Energy functions         | Compute kinetic, potential, and total energy per particle |

## Usage

### Requirements

You can install the required dependencies using:

```bash
pip install numpy matplotlib
```

### Running the Simulation

Simply execute the script. It will:

- Set the parameters $(N,\; h,\; n = 4)$  
- Generate an initial condition
- Evolve the system for $t$ time units  
- Compute the largest Lyapunov exponent and store its time evolution in `L`  
- Print the final exponent  

## Results & Discussion

For the chosen parameters, the system typically exhibits **chaotic behavior**, indicated by a positive Lyapunov exponent.

The finite-time Lyapunov exponent converges to a positive value, confirming the **exponential divergence of nearby trajectories**.

---

## What I Learned

Working on this project allowed me to develop skills in:

- Implementing high-order symplectic integrators  
- Understanding the importance of preserving geometric structures in long-time simulations  
- Deriving and coding linearized equations for tangent-space dynamics  
- Computing Lyapunov exponents (renormalization, Gram–Schmidt methods)  
- Debugging and optimizing scientific Python code  

---

## References

- Suzuki, M. (1990). *Fractal decomposition of exponential operators with applications to many-body theories and Monte Carlo simulations*. Physics Letters A, 146(6), 319–323.

- Benettin, G., Galgani, L., Giorgilli, A., & Strelcyn, J.-M. (1980). *Lyapunov characteristic exponents for smooth dynamical systems and for Hamiltonian systems; a method for computing all of them*. Meccanica, 15(1), 9–20.

- Antoni, M., & Ruffo, S. (1995). *Clustering and relaxation in Hamiltonian long-range dynamics*. Physical Review E, 52(3), 2361.

- Maciel, J. M., Firpo, M.-C., & Amato, M. A. (2015). *Some statistical equilibrium mechanics and stability properties of a class of two-dimensional Hamiltonian mean-field models*.  
  CNRS–École Polytechnique (France) & Universidade de Brasília (Brazil).

## Final Remarks

This project was essential for my growth as a computational physicist. It combined theoretical concepts from classical mechanics and chaos theory with hands-on implementation, providing a solid foundation for future research in complex systems.
