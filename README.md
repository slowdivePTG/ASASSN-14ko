# ASASSN-14ko
New project with Enrico and Jamie at UCSC.

## Sep. 17, 2020

[ASASSN-14ko is a Periodic Nuclear Transient in ESO 253−G003](https://arxiv.org/abs/2009.03321)

The orbital period is about 100 days. The idea will that this star was in an eccentric bound orbit, had a stellar interaction that slightly altered its orbit and redirected towards the black hole. 

[A Brief Review of TDEs](TDE.md)

## Oct. 15, 2020

1. [Toy models for a TDE event](./TDE_Analytical) (parabolic initial orbit)

2. FLASH & MESA documents from Jamie

   >https://github.com/jamielaw-smith/STARS-sims-doc are the main instructions
   >
   >https://github.com/jamielaw-smith/multitidal is the FLASH setup source code
   >
   >https://github.com/jamielaw-smith/jYT and https://github.com/jamielaw-smith/realistic-structure are analysis code.
   >
   >https://github.com/jamielaw-smith/FLASH-mentoring-instructions is a catch-all, copy and paste of misc instructions

3. [MESA tutorial](http://mesa.sourceforge.net/index.html)

## Nov. 4, 2020

1. Sink particles ([Federrath et al. 2010](https://iopscience.iop.org/article/10.1088/0004-637X/713/1/269/pdf))

   If the gas has reached a given density, a sink particle is introduced, which can accrete the gas exceeding the threshold, without al- tering the thermal physics.

   For a collapsing system, the free-fall timescale $t_\text{ff}$ decreases with increasing density
   $$
   t_\text{ff}=\sqrt{\frac{3\pi}{32G\rho}}
   $$
   Following the free-fall collapse/accretion, the change in density could cover over ten orders of magnitude. Such high dynamical range is beyond the modern numerical schemes and supercomputers. In `FLASH`, the `SinkParticle` unit can identify the dense yet gravitationally bound gas and turn it into sink particles.

   Once created, sink particles are free to move within the Cartesian computational domain, independent of the underlying grid, i.e., they move in the Lagrangian frame of reference, while the grid points are fixed in space (Eulerian frame of reference). The outer boundary conditions, i.e., *outflow*, *reflecting*, or *periodic* also apply to the sink particles in the simulation box.