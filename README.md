# Atmospheric Dynamics Lab

**[Open the live laboratory](https://genlier.github.io/atmospheric-dynamics-lab/)** · [中文版](README.zh-CN.md) · [中文网页](https://genlier.github.io/atmospheric-dynamics-lab/?lang=zh)

Explore atmospheric flow using Fourier spectral models. English is the default website language; Chinese is available from the language switch.

## Interactive laboratory

Paint positive or negative vorticity, erase part of the field, enter a mathematical function, or import a JSON grid. Choose the Reynolds number, beta, integration duration, maximum time step and output density, then compute an animation in your browser.

- Periodic domain: $[-\pi,\pi)\times[-\pi,\pi)$.
- Grid: 32×32 or 64×64; 61, 121 or 241 computed frames.
- Pause, scrub the timeline, change playback speed, inspect grid values, and view kinetic energy and enstrophy.
- Export the initial field or computed results as JSON. Imported grids use the same initial-field format.
- Computation and imported data stay in the browser.

The live solver integrates the **barotropic vorticity equation**. The nine reference experiments also include shallow-water adjustment, two-layer quasi-geostrophy and a rising thermal bubble; these supplementary models are presented as saved results.

Reference playback interpolates the existing snapshots. The live laboratory numerically integrates a new trajectory and interpolates adjacent computed frames for playback. Use **“Use this initial field in the live lab”** on the five BVE reference cases to compute a denser trajectory. The live calculation uses its displayed grid and settings; the 256×256 turbulence reference is sampled to 64×64 for display.

## Run locally

Serve the website from the repository root:

~~~bash
python -m http.server 8000 --directory docs
~~~

Open http://localhost:8000. All scripts are included locally; no build step or external resources are required. Worker computation runs from HTTP or HTTPS pages.

### MATLAB

~~~matlab
run('models/matlab/bve.m')
~~~

Set p.model_case to 1–5 for the BVE benchmarks. models/matlab/shallow_water.m uses p.exp_case to select initial height or vorticity perturbations.

### Python

~~~bash
python -m pip install -r requirements.txt
python models/python/bve.py --case 1 --output outputs/case1.npz
python models/python/experiments.py
python models/python/supplementary.py
~~~

The first script saves numerical arrays and, for the decaying vortex grid, numerical and exact energy histories. The auxiliary scripts explore asymmetric vortices, shallow water, two-layer QG and thermals. Their configurations are defined in the source.

## Mathematical model

The live model uses

$$
\frac{\partial\zeta}{\partial t}+J(\psi,\zeta)+\beta\frac{\partial\psi}{\partial x}
=\frac{1}{\mathrm{Re}}\nabla^2\zeta,\qquad
\nabla^2\psi=\zeta,
$$

with $u=-\partial_y\psi$, $v=\partial_x\psi$ and $J(\psi,\zeta)=\partial_x\psi\,\partial_y\zeta-\partial_y\psi\,\partial_x\zeta$.

Fourier differentiation represents spatial derivatives through wavenumber factors. Poisson inversion uses the zero-mean vorticity field, with the streamfunction's zero mode set to zero. Nonlinear products use 2/3 truncation. Advection advances with RK4 and viscous diffusion with an exponential factor. The browser solver bounds its time step using the selected maximum, the velocity CFL estimate and the beta-wave frequency.

Initial JSON fields contain n and values, ordered with x increasing fastest, followed by y increasing from −π. Results include the initial field, parameters, mean removed, completion state, time, vorticity, energy and enstrophy for each computed frame.

## Verification

~~~bash
node tests/test_browser_solver.js
python -m unittest discover -s tests
python tools/audit_publication.py
~~~

The JavaScript tests check FFT inversion, exact viscous decay, beta-wave phase, nonlinear energy and enstrophy conservation, zero-mean inversion and expression parsing. Browser checks cover both languages, drawing, computation, downloads, imports, cancellation and mobile layout.

## Repository

| Directory | Contents |
| --- | --- |
| models/matlab/ | Original BVE and shallow-water solvers |
| models/python/ | NumPy solver and supplementary experiments |
| docs/ | Bilingual site, reference snapshots and browser solver |
| tests/ | Analytic and conservation checks |
| tools/ | Publication audit |

The publication contains selected code and numerical snapshots. Personal reports, lecture files, credentials, private paths and raw datasets are excluded. Commits use the owner's GitHub privacy email.

