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

- Tidal Radius $r_\text T$

  Essentially, $r_\text T$ is the distance from the hole at which $M_\text H/ r^3$ equals the mean internal density of the passing star. Here $M_\text H$ is the mass of the BH.
  $$
  r_\text T=r_*\left(\frac{M_\text H}{M_*}\right)^{1/3}
  $$
  In order of magnitude, $r_\text T$ is the same as the Roche radius: the latter is however only applicable to a star in a circular orbit with synchronized spin.

  When a star reaches the tidal radius of a BH, it could either lose its outer layers or suffers from complete disruption. Even if a solar-type star loses no material during each pericentric passage, it could get so distorted that the subsequent internal dissipation might exceed its incoming kinetic energy. Thus the star is captured by the BH.

- The destruction of a star

  The binding energy of a star with radius $r_*$ and mass $m_*$ is
  $$
  U=\frac{3Gm_*^2}{5r_*}
  $$
  For a sun-type star at $r_\text T$ of a SMBH in a parabolic orbit, the kinetic energy is
  $$
  K=\frac{GM_\text Hm_*}{r_\text T}
  $$
  So
  $$
  \frac{K}{U}\sim \frac{r_*}{r_\text T}\frac{M_\text H}{m_*}=\left(\frac{M_\text H}{m_*}\right)^{2/3}\gg 1
  $$
  The energy required to tear the star apart (that is, the star's self-binding energy $U$) is supplied at the expense of the orbital kinetic energy $K$.

  Unless there were some explosive energy input, on average the debris would be bound to the hole, unless the star was initially on a hyperbolic orbit with asymptotic velocity over the escape velocity $v_*$ at the surface of the star, which is 1,000 km/s for solar-type stars.

  There are some other effects. The dominant effect is that while falling inwards towards the hole, the star would develop a **quadrupole distortion** which attains an amplitude of order unity by the time of disruption. Then the gravitational torque would spin it up and by the time it gets disrupted, it could be **spinning at close to its break-up angular velocity**.

  At $r_\text T$, the orbital velocity is approximately
  $$
  v_{orb}\sim\sqrt\frac{2GM_\text H}{r_\text T}\sim c\sqrt{\frac{r_g}{r_\text T}}
  $$
  where $r_g$ is the Schwarzschild radius, while the escape velocity
  $$
  v_*\sim\sqrt{\frac{Gm_*}{r_*}}\sim\left(\frac{m_*}{M_\text H}\right)^{1/3}v_{orb}
  $$
  In this way, the gas on the 'outer track' furthest from the BH moves $\sim\left({m_*}/{M_\text H}\right)^{1/3}v_{orb}$ faster than the 'inside track'.

  Morever, the slower-moving 'inside track' lies deeper in the BH's potential by an amount
  $$
  \Delta\Phi\sim -r_*\nabla_r\left(\frac{GM_\text H}{r}\right)_{r_\text T}=\frac{GM_\text H r_*}{r_\text T^2}=\frac{Gm_*}{r_*}\left(\frac{M_\text H}{m_*}\right)^{1/3}\sim v_{orb}v_*
  $$
  Even though the mean specific binding energy of the debris to the hole would be positive, and comparable with the self-binding energy ($Gm_*/r_*$) of the original star (**an assumption**, with which the initial orbit is slightly elliptical), the spread about this mean is larger by $(M_\text H/ m_*)^{1/3}$.
  
  In general, the original star moves in a parabolic, and thus zero-energy orbit. With this potential dispersion, however, part of the star is bound to the BH, while the other would be left unbound. For the outermost gas, **the specific orbital energy** is approximately
  $$
  -v_*^2+v_{orb}v_*=v_*(-v_*+v_{orb})=\frac{Gm_*}{r_*}\left[\left(\frac{M_\text H}{m_*}\right)^{1/3}-1\right]
  $$
  for the innermost gas, this turns to
  $$
  -v_*^2-v_{orb}v_*=-v_*(v_*+v_{orb})=-\frac{Gm_*}{r_*}\left[\left(\frac{M_\text H}{m_*}\right)^{1/3}+1\right]
  $$
  
  Here, $-v_*^2$ accounts for the energy loss in the distruption of the star.

## Fallback

**The average specific binding energy** of debris bound to the BH is
$$
\sim\frac{1}{2}\frac{Gm_*}{r_*}\left[\left(\frac{M_\text H}{m_*}\right)^{1/3}+1\right]
$$
The bound orbits are close to parabolics, and are thus very eccentric. In a rough estimation, since
$$
E=-\frac{GM_\text Hm_*}{2a}
$$
where $a$ is the semi-major axis, so for the most tightly bound debris, we have
$$
a\sim\left(\frac{M_\text H}{m_*}\right)^{2/3}r_*
$$
For a sun-type star, $a\sim10^3\ M_6^{-1/3}r_g$, where $M_6=M_\text H/(10^6 M_\odot)$.

The orbital period is
$$
P=\sqrt\frac{a^3}{GM_\text H}\sim\sqrt{\frac{M_\text Hr_*^3}{Gm_*^2}}\sim0.05\ M_6^{1/2}\left(\frac{m_*}{M_\odot}\right)^{-1}\left(\frac{R_*}{R_\odot}\right)^{3/2}\text{ yr}
$$
while the initial orbital period is $\sim 20\ M_6\text{ yr}$, a factor of $400\ M_6^{1/2}$ longer for a solar-type star.

This orbital period is much shorter compared to the estimated interval between each TDE ($\sim 10^{3-4}\text{ yr}$). Unless it takes many orbital periods to swallow the bound gas, the debris from each star would be digested separately.

In order to form an accretion flow, the bound stellar debris must lose a significant amount of energy by viscous dissipation. If the viscosity is large enough to allow accretion onto the SMBH on a timescale shorter than $T_t$, the luminosity of the flare is expected to follow the rate of mass fallback
$$
\dot M=\frac{\text dM}{\text de}\frac{\text de}{\text da}\frac{\text da}{\text dT}\sim\frac{GM_\text H}{a^2}\cdot \frac{a}{T}\frac{\text dM}{\text de}\sim\frac{\text dM}{\text de}\left(GM_\text H\right)^{2/3}T^{-5/3}
$$
by assuming $\text dM/\text de\sim m_*/e_t$ is almost a constant (here $e$ refers to the specific energy).



## How Many Parameters Do We Need

Ignoring general relativistic effects and stellar rotation, and assuming stars to be polytropes, it may seem that a complete study of tidal disruptions would require an exhaustive study of the various combinations of six parameters: $M_∗$, $R_∗$, $M_\text h$, the orbital eccentricity $e$, the polytropic index $\gamma$, and $\beta\equiv  r_\text T/r_\text p$.

With $M_*$ and $\gamma$, we can without too much effort derive $R_*$.

We define a crossing timescale $t_\text p\equiv {r_\text p}/{v_\text{p}}$, where $\text p$ denotes for the pericenter. It reflects how long the encounter is. Since most of the stars that are scattered into disruptive orbits originate from the sphere of influence or beyond, $e\sim1$, thus the orbit is nearly parabolic. Fix $\beta$ and $M_*$, the definition of $\beta$ directly gives
$$
r_\text p=\frac{r_\text T}{\beta}=\frac1\beta R_*\left(\frac{M_\text H}{M_*}\right)^{1/3}\propto M_\text H^{1/3}
$$

And for a nearly parabolic orbit,
$$
v_\text p=\sqrt{\frac{2GM_\text H}{r_\text p}}\propto M_\text H^{1/3}
$$
As a result, $t_\text p$ is independent of $M_\text H$ for fixed $\beta$ and $M_*$.

How the encounter influence the star depends on the ratio of $t_\text p$ and the star's dynamic timescale
$$
t_\text{dyn}=\sqrt{\frac{R_*^3}{GM_*}}
$$
Surprisingly,
$$
\frac{t_\text p}{t_\text{dyn}}=\sqrt{\frac{r_\text p^3}{2GM_\text H}}\cdot\sqrt{\frac{GM_*}{R_*^3}}=\sqrt{\frac{\beta^3}{2}}
$$
is even independent of the stellar mass and radius.

Additionally, as the mass ratio approaches infinity, the asymmetry of the tidal field becomes progressively less important as $R_∗ \ll r_\text T$, with the difference in the strength of the tidal field at pericenter between the near side and the far side for a $10^6:1$ encounter being $\simeq3\%$ (Guillochon et al. 2011).

With all these facts above, it is claimed that the vast majority of stellar disruptions by SMBHs can be described by just two parameters: $\beta$ and $\gamma$.



## Eccentric Orbit

In an eccentric orbit, the object is more bound to the black hole than that in a parabolic orbit. The orbital energy for the most bound material is
$$
\epsilon=-\frac{GM_\text H}{2a_*}-\frac{GM_\text H}{r_\text T^2}R_*=-\frac{GM_\text H}{2a}
$$

$$
\begin{align*}
\Rightarrow P&=\frac{\pi}{\sqrt 2}GM_\text H(-\epsilon)^{-3/2}=\frac{\pi}{\sqrt 2}GM_\text H\left(\frac{GM_\text H}{2a_*}+\frac{GM_\text H}{r_\text T^2}R_*\right)^{-3/2}\\
&=\frac{\pi}{\sqrt 2}(GM_\text H)^{-1/2}\left(\frac{\beta_*(1-e_*)}{2r_\text T}+\frac{1}{r_\text T^2}R_*\right)^{-3/2}\\
&=\frac{\pi}{\sqrt 2}\left(\frac{GM_\text H}{R_*^3}\right)^{-1/2}\left(\frac{\beta_*(1-e_*)}{2}\left(\frac{M_\text H}{M_*}\right)^{-1/3}+\left(\frac{M_\text H}{M_*}\right)^{-2/3}\right)^{-3/2}\\
&=\frac{\pi}{\sqrt 2}\left(\frac{GM_*}{R_*^3}\right)^{-1/2}\left(\frac{M_\text H}{M_*}\right)^{1/2}\left(\frac{\beta_*(1-e_*)}{2}\left(\frac{M_\text H}{M_*}\right)^{1/3}+1\right)^{-3/2}\\
&=\frac{\pi}{\sqrt 2}\left(\frac{GM_*}{R_*^3}\right)^{-1/2}\beta_*^{-3/2}(1-e_*)^{-3/2}\left(\frac{1}{2}+\left(\frac{M_\text H}{M_*}\right)^{-1/3}\beta_*^{-1}(1-e_*)^{-1}\right)^{-3/2}
\end{align*}
$$

while the orbital period for the CoM is
$$
\begin{align*}P_*&=\frac{\pi}{\sqrt 2}\left(\frac{GM_*\beta_*^3}{R_*^3}\right)^{-1/2}\left(\frac{1-e_*}{2}\right)^{-3/2}\\
&=2\pi\left(\frac{r_\text p^3(1-e_*)^3}{GM_\text H}\right)^{1/2}\\
&=2\pi\left(\frac{a_*^3}{GM_\text H}\right)^{1/2}
\end{align*}
$$
This is simply Kepler's III law. Interestingly we have
$$
\left(P^{-2/3}-P_*^{-2/3}\right)^{-3/2}=\frac{\pi}{\sqrt 2}\left(\frac{GM_*}{R_*^3}\right)^{-1/2}\left(\frac{M_\text H}{M_*}\right)^{1/2}=\frac{\pi}{\sqrt 2}\left(\frac{M_\text H}{M_*}\right)^{1/2}t_\text{dyn}
$$
which is the $t_\text{peak}$ in the parabolic case.

Similarly, the outermost trajectory has an orbital energy of
$$
\epsilon=-\frac{GM_\text H\beta_*(1-e_*)}{2r_\text T}+\frac{GM_\text H}{r_\text T^2}R_*
$$
If it is also bound to the black hole, we have
$$
e_*<e_\text{crit}\equiv1-\frac2{\beta_*}\left(\frac{M_\text H}{M_*}\right)^{-1/3}
$$
The duration time ofmass fallback for eccentric TDEs with $e_*<e_\text{crit}$ is finite, and if $e_*\lesssim e_\text{crit}$,
$$
\begin{align*}
\Delta t&=\frac{\pi}{\sqrt 2}\left(\frac{GM_*}{R_*^3}\right)^{-1/2}\beta_*^{-3/2}(1-e_*)^{-3/2}\left[\left(\frac{1}{2}-\left(\frac{M_\text H}{M_*}\right)^{-1/3}\beta_*^{-1}(1-e_*)^{-1}\right)^{-3/2}-1\right]
\end{align*}
$$