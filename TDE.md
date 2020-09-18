# Tidal Disruption Event (TDE)

> The gravity of the black hole (BH) introduces strong tidal forces that can deform, mangle, and potentially destroy objects that approach it too closely.
>
> There are a number of phases of a tidal disruption that are potentially observable by Earth- and space-based observatories.
>
> 1. **Scattering of a star into the tidal radius.**
> 2. **The moment the star is destroyed or damaged: Crossing the tidal radius.**
> 3. **The fallback of debris back onto the supermassive black hole.**
>
> ----By James Guillochon

## Scattering of A Star

### Hill's Mechanism

### Lidov-Kozai Effect

To be continued...

## Crossing the Tidal Radius

- Tidal Radius $r_T$

  Essentially, $r_T$ is the distance from the hole at which $M_h/ r^3$ equals the mean internal density of the passing star. Here $M_h$ is the mass of the BH.
  $$
  r_T=r_*\left(\frac{M_h}{M_*}\right)^{1/3}
  $$
  In order of magnitude, $r_T$ is the same as the Roche radius: the latter is however only applicable to a star in a circular orbit with synchronized spin.

  When a star reaches the tidal radius of a BH, chances are that it could either lose its outer layers or suffers from complete disruption. Even if a solar-type star loses no material during each pericentric passage, it could get so distorted that the subsequent internal dissipation might exceed its incoming kinetic energy. Thus the star is captured by the BH.

- The destruction of a star

  The binding energy of a star with radius $r_*$ and mass $m_*$ is
  $$
  U=\frac{3Gm_*^2}{5r_*}
  $$
  For a sun-type star at $r_T$ of a SMBH in a parabolic orbit, the kinetic energy is
  $$
  K=\frac{GM_hm_*}{r_T}
  $$
  So
  $$
  \frac{K}{U}\sim \frac{r_*}{r_T}\frac{M_h}{m_*}=\left(\frac{M_h}{m_*}\right)^{2/3}\gg 1
  $$
  The energy required to tear the star apart (that is, the star's self-binding energy $U$) is supplied at the expense of the orbital kinetic energy $K$.

  Unless there were some explosive energy input, on average the debris would be bound to the hole, unless the star was initially on a hyperbolic orbit with asymptotic velocity over the escape velocity $v_*$ at the surface of the star, which is 1,000 km/s for solar-type stars.

  There are some other effects. The dominant effect is that while falling inwards towards the hole, the star would develop a **quadrupole distortion** which attains an amplitude of order unity by the time of disruption. Then the gravitational torque would spin it up and by the time it gets disrupted, it could be **spinning at close to its break-up angular velocity**.

  At $r_T$, the orbital velocity is approximately
  $$
  v_{orb}\sim\sqrt\frac{2GM_h}{r_T}\sim c\sqrt{\frac{r_g}{r_T}}
  $$
  where $r_g$ is the Schwarzschild radius, while the escape velocity
  $$
  v_*\sim\sqrt{\frac{Gm_*}{r_*}}\sim\left(\frac{m_*}{M_h}\right)^{1/3}v_{orb}
  $$
   In this way, the gas on the 'outer track' furthest from the BH moves $\sim\left({m_*}/{M_h}\right)^{1/3}v_{orb}$ faster than the 'inside track'.

  Morever, the slower-moving 'inside track' lies deeper in the BH's potential by an amount
  $$
  \Delta\Phi\sim -r_*\nabla_r\left(\frac{GM_h}{r}\right)_{r_T}=\frac{GM_h r_*}{r_T^2}=\frac{Gm_*}{r_*}\left(\frac{M_h}{m_*}\right)^{1/3}\sim v_{orb}v_*
  $$
  Even though the mean specific binding energy of the debris to the hole would be positive, and comparable with the self-binding energy ($Gm_*/r_*$) of the original star (**an assumption**, with which the initial orbit is slightly elliptical), the spread about this mean is larger by $(M_h/ m_*)^{1/3}$.
  
  In general, the original star moves in a parabolic, and thus zero-energy orbit. With this potential dispersion, however, part of the star is bound to the BH, while the other would be left unbound. For the outermost gas, **the specific orbital energy** is approximately
  $$
  -v_*^2+v_{orb}v_*=v_*(-v_*+v_{orb})=\frac{Gm_*}{r_*}\left[\left(\frac{M_h}{m_*}\right)^{1/3}-1\right]
  $$
  for the innermost gas, this turns to
  $$
  -v_*^2-v_{orb}v_*=-v_*(v_*+v_{orb})=-\frac{Gm_*}{r_*}\left[\left(\frac{M_h}{m_*}\right)^{1/3}+1\right]
  $$

## Fallback

**The average specific binding energy** of debris bound to the BH is
$$
\sim\frac{1}{2}\frac{Gm_*}{r_*}\left[\left(\frac{M_h}{m_*}\right)^{1/3}+1\right]
$$
The bound orbits are close to parabolics, and are thus very eccentric. In a rough estimation, since
$$
E=-\frac{GM_hm_*}{2a}
$$
where $a$ is the semi-major axis, we have
$$
a\sim\left(\frac{M_h}{m_*}\right)^{2/3}r_*
$$
For a sun-type star, $a\sim10^3\ M_6^{-1/3}r_g$, where $M_6=M_h/(10^6 M_\odot)$.

Now that $r_{min}=r_T$, 
$$
\left(\frac{M_h}{m_*}\right)^{2/3}r_*(1-e)=r_T=\left(\frac{M_h}{m_*}\right)^{1/3}r_*
$$
thus
$$
e=1-\left(\frac{M_h}{m_*}\right)^{-1/3}\sim1
$$
The orbital period is
$$
T_{debris}=\sqrt\frac{a^3}{GM_h}\sim\sqrt{\frac{M_Hr_*^3}{Gm_*^2}}\sim0.05\ M_6^{1/2}\text{ yr}
$$
while the initial orbital period is $\sim 20\ M_6\text{ yr}$, a factor of $400\ M_6^{1/2}$ longer.

This orbital period is much shorter compared to the estimated interval between each TDE ($\sim 10^{3-4}\text{ yr}$). Unless it takes many orbital periods to swallow the bound gas, the debris from each star would be digested separately.