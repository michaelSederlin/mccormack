# McCormack 

Implements the McCormack Scheme for solving a PW type traffic flow model

This repository provides several predefined traffic scenarios for testing and demonstration of numerical simulation methods (e.g., McCormack scheme). Each scenario returns an initial density (`rho`), velocity (`v0`), and simulation parameters (`params`).

Please notify me about any bugs and/or other mistakes.

---

## Scenarios

### 1. `backward_wave_example()`
- Simulates a backward-propagating wave.
- Parameters: low congestion $(\rho_0 \approx 30)$, single Gaussian bump.

### 2. `gaussian_density_example()`
- Base case with a Gaussian bump in density.
- $\rho_0 \approx 30$, peak added at center.

### 3. `with_disturbance()`
- Adds a disturbance during the simulation at defined extent in space and time.
- Based on `gaussian_density_example()`.
- Can reduce jam density ($\rho_{jam}$) or decrease density.

### 4. `numerical_oscilations()`
- Tuned to demonstrate numerical instabilities.
- Uses $\Delta t$ slightly below CFL condition $\Delta x / V_0$ and small $\tau$
- Purpose: highlight need for stability conditions.

### 5. `lecture_high_res()`
- High-resolution grid for lecture demos.
- Sharp localized density spike.

### 6. `lecture_low_res()`
- Low-resolution version of the above.

Scenarios are packaged in functions named as above and returns parameters 

```python
rho, v0, params = scenario()
or 
rho, v0, params, disturbances = scenario_with_disturbance()
```
