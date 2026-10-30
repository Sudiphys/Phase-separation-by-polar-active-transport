# Documentation: Phase separation by polar active transport

The program integrates Eqs. (S10a–c) on a two-dimensional periodic square grid using centered spatial differences and explicit forward Euler time integration. This document describes the supplied source without modifying it.

## Fields and coefficients

| Code variable | Meaning |
|---|---|
| `x(i,j,iold)` | Microtubule area fraction, $\psi$ |
| `cp(i,j,iold)` | Plus-motor area fraction, $\hat c_+$ |
| `cm(i,j,iold)` | Minus-motor area fraction, $\hat c_-$ |
| `plx(i,j,iold)`, `ply(i,j,iold)` | Polarity components, $p_x,p_y$ |
| `cdif(i,j,iold)` | Motor difference, $\delta c=\hat c_+-\hat c_-$ |
| `atot` | Total area fraction, $s=\psi+\hat c_++\hat c_-$ |
| `polOP` | Squared polarity magnitude, $|\mathbf p|^2$ |
| `iold`, `inew` | Indices of the old and updated time buffers |

The coefficients entering the evolution equations are mapped as follows:

| Code variable | Definition in the code | Corresponding coefficient in S10 |
|---|---|---|
| `Dmt` | `Dpsi/(Dr*lp*lp)` | $\widetilde D_{mt}$ |
| `Dm` | `Dc/(Dr*lp*lp)` | $\widetilde D$ |
| `vp` | `v0/(Dr*lp)` | $\widetilde v_0$ |
| `ara` | `Dc*DLOG(Le/d)/(pi*d*v0)` | $1/\mathfrak R$ |
| `anup` | `anu*arho/lp` | Effective coupling $\widetilde\nu$ appearing in S10 |
| `azetap` | `azeta*arho/(lp*lp)` | $2\widetilde\zeta/(\hat c_0^++\hat c_0^-)$ |
| `akap` | Directly assigned | $\widetilde\kappa/\psi_0$ |
| `bet`, `alpha` | Directly assigned | $\beta,\alpha$ |
| `lsq` | `2*Le` | $2L$ |
| `arho` | `1/(pi*0.02D0*0.02D0)` | Reference concentration $\rho_0$ for the specified radius |

`lp` is the length used to scale speed and diffusion; `Le` is the microtubule length in the geometric factors. Both equal 5 in this source. `d` is the motor radius. `Dr` is the rotational diffusion coefficient; `Dpsi` and `Dc` are the input microtubule and motor diffusion coefficients.

The local velocity variable has the convention

$$\texttt{vm}=\texttt{vp}\frac{\hat c_+-\hat c_-}{\hat c_++\hat c_-}=-\widetilde v_{mt}.$$

It is calculated at every grid point and stored in `vmfield`. If the total motor fraction is zero, `vmfield` is set to zero.

## Spatial derivatives

`p(i)` and `m(i)` give the forward and backward periodic neighbors. They are index arrays, distinct from the polarity field $\mathbf p$.

For any scalar field $f$,

$$D_xf_{i,j}=\frac{f_{p(i),j}-f_{m(i),j}}{2\,\texttt{dx}},\qquad D_yf_{i,j}=\frac{f_{i,p(j)}-f_{i,m(j)}}{2\,\texttt{dx}}.$$

The Laplacian is

$$\Delta_hf_{i,j}=\frac{f_{p(i),j}+f_{m(i),j}+f_{i,p(j)}+f_{i,m(j)}-4f_{i,j}}{\texttt{dx}^2}.$$

`dx2=2*dx` is the denominator for first derivatives; `dxx=dx**2` is the denominator for second derivatives. Gradients and divergences of intermediate fields are obtained by applying the centered first differences successively. In the mappings below, $\nabla$ and $\nabla\cdot$ denote these centered operators, and $\Delta_h$ denotes the five-point Laplacian.

## Motor equations: S10a

The active transport vectors inside the right-hand-side divergences are

$$\mathbf V_+=\texttt{ara}(-\texttt{vp}+\texttt{vm})\psi\hat c_+\mathbf p,\qquad \mathbf V_-=\texttt{ara}(\texttt{vp}+\texttt{vm})\psi\hat c_-\mathbf p.$$

The gradient-energy and torque contributions to the motor chemical potentials are

$$\mu_+^{g}=-\texttt{azetap}\,\Delta_h\hat c_++\texttt{anup}\frac{\pi d}{\texttt{lsq}}\nabla\cdot(\psi\mathbf p),$$
$$\mu_-^{g}=-\texttt{azetap}\,\Delta_h\hat c_--\texttt{anup}\frac{\pi d}{\texttt{lsq}}\nabla\cdot(\psi\mathbf p).$$

| Code expression | Calculated quantity |
|---|---|
| `conpxp`, `conpyp` | Components of $\mathbf V_+$ |
| `conmxm`, `conmym` | Components of $\mathbf V_-$ |
| `conpxxp+conpyyp` | $\nabla\cdot\mathbf V_+$ |
| `conmxxm+conmyym` | $\nabla\cdot\mathbf V_-$ |
| `acpxx+acpyy` | $\Delta_h\hat c_+$ |
| `acmxx+acmyy` | $\Delta_h\hat c_-$ |
| `agdcpxx`, `agdcpyy` | Components of $\hat c_+\nabla\hat c_+/(1-s)$ |
| `agdcmxx`, `agdcmyy` | Components of $\hat c_-\nabla\hat c_-/(1-s)$ |
| `agdcpsxxx+agdcpsyyy` | $\nabla\cdot[\hat c_+\nabla\hat c_+/(1-s)]$ |
| `agdcmsxxx+agdcmsyyy` | $\nabla\cdot[\hat c_-\nabla\hat c_-/(1-s)]$ |
| `bpsx`, `bpsy` | Components of $\psi\mathbf p$ |
| `bpsxx+bpsyy` | $\nabla\cdot(\psi\mathbf p)$ |
| `amucp`, `amucm` | $\mu_+^{g},\mu_-^{g}$ |
| `ajcpx`, `ajcpy` | Components of $\hat c_+\nabla\mu_+^{g}$ |
| `ajcmx`, `ajcmy` | Components of $\hat c_-\nabla\mu_-^{g}$ |
| `divcpx+divcpy` | $\nabla\cdot(\hat c_+\nabla\mu_+^{g})$ |
| `divcmx+divcmy` | $\nabla\cdot(\hat c_-\nabla\mu_-^{g})$ |

The complete motor right-hand sides are

$$\texttt{addtotcp}=\nabla\cdot\mathbf V_++\texttt{Dm}\left[\Delta_h\hat c_++\nabla\cdot\left(\frac{\hat c_+}{1-s}\nabla\hat c_+\right)+\nabla\cdot(\hat c_+\nabla\mu_+^{g})\right],$$

$$\texttt{addtotcm}=\nabla\cdot\mathbf V_-+\texttt{Dm}\left[\Delta_h\hat c_-+\nabla\cdot\left(\frac{\hat c_-}{1-s}\nabla\hat c_-\right)+\nabla\cdot(\hat c_-\nabla\mu_-^{g})\right].$$

These contain active transport, ordinary diffusion, excluded-volume diffusion, motor gradient energy, and torque coupling. Vectors inside a positive right-hand-side divergence have the opposite sign to the physical current in $\partial_tc=-\nabla\cdot\mathbf j$.

## Microtubule equation: S10b

Define

$$q=\nabla\cdot\mathbf p,\qquad A=|\mathbf p|^2\left(1+\frac\alpha2|\mathbf p|^2\right),\qquad B=\mathbf p\cdot\nabla\delta c,\qquad H=\psi q^2.$$

| Code expression | Calculated quantity |
|---|---|
| `conpx`, `conpy` | Components of $\texttt{vm}\,\psi\mathbf p$ |
| `addcompsi` | $\nabla\cdot(\texttt{vm}\,\psi\mathbf p)$: active transport |
| `agrp2` | $A$ |
| `divpsx`, `divpsy` | Components of $\psi\nabla A$ |
| `dpcdf` | $B$ |
| `dpcdfpsi1`, `dpcdfpsi2` | Components of $\psi\nabla B$ |
| `aggp2` | $H$ |
| `aggpx2`, `aggpy2` | Components of $\nabla H$ |
| `aggpx2s`, `aggpy2s` | Components of $\psi\nabla H$ |
| `agdpx`, `agdpy` | Components of $(\texttt{lsq}/\pi d)\psi\nabla\psi/(1-s)$ |

The individual density contributions are

| Code term | Mathematical expression | Meaning |
|---|---|---|
| `tlapsi` | $\Delta_h\psi+\nabla\cdot[(2L/\pi d)\psi\nabla\psi/(1-s)]$ | Ordinary and excluded-volume diffusion |
| `tbetpsi` | $(\texttt{bet}/2)\nabla\cdot(\psi\nabla A)$ | Local polarity-energy contribution |
| `tnupsict` | $-\texttt{anup}\nabla\cdot(\psi\nabla B)$ | Torque-coupling contribution |
| `tkapsi` | $\texttt{bet}\,\texttt{akap}\nabla\cdot(\psi\nabla H)$ | Bending contribution |

The code assembles

```fortran
addpsitotal=(tlapsi+tbetpsi+tnupsict+tkapsi)*Dmt
```

The complete microtubule right-hand side is `addpsitotal+addcompsi`. The active transport term is added outside the `Dmt` factor.

## Polarity equation: S10c

The polarity right-hand side is calculated in the compact form

$$\partial_t\mathbf p=\texttt{vm}(\mathbf p\cdot\nabla)\mathbf p-\mathbf p(1+\alpha|\mathbf p|^2)+\frac{\texttt{anup}}{\texttt{bet}}\nabla\delta c+\frac{\texttt{akap}}{\psi}\nabla(\psi^2q).$$

At the continuum level, the compact bending term equals

$$\frac{\texttt{akap}}\psi\nabla(\psi^2q)=\texttt{akap}\left[2q\nabla\psi+\psi\nabla q\right],$$

which is the form in S10c with `akap` corresponding to $\widetilde\kappa/\psi_0$. The code differences the compact product directly.

| Code variable | Calculated quantity |
|---|---|
| `pxx`, `pxxy` | $D_xp_x,D_yp_x$ |
| `pyyx`, `pyy` | $D_xp_y,D_yp_y$ |
| `pspx` | $\psi^2q$ |
| `pxpx`, `pypy` | $D_x(\psi^2q),D_y(\psi^2q)$ |
| `tpbetax`, `tpbetay` | Components of $-\psi\mathbf p(1+\alpha|\mathbf p|^2)$ |
| `tnux`, `tnuy` | Components of $\texttt{anup}\,\psi\nabla\delta c$ |
| `tpkapx`, `tpkapy` | Components of $\texttt{akap}\nabla(\psi^2q)$ |
| `advpx`, `advpy` | Components of $-\texttt{vm}(\mathbf p\cdot\nabla)\mathbf p$ |

The program forms `addpx` and `addpy` by adding the relaxation, coupling divided by `bet`, and bending terms, subtracting `advpx*x` or `advpy*x`, and finally dividing by `x`. Thus the advection contribution after assembly is $+\texttt{vm}(\mathbf p\cdot\nabla)\mathbf p=-\widetilde v_{mt}(\mathbf p\cdot\nabla)\mathbf p$.

## Time integration and initial conditions

Each timestep first calculates `vmfield`, then completes three full-grid intermediate passes. The final pass assembles all right-hand sides from the old fields. The updates are

```fortran
x(i,j,inew)=x(i,j,iold)+(addpsitotal+addcompsi)*dt
cp(i,j,inew)=cp(i,j,iold)+addtotcp*dt
cm(i,j,inew)=cm(i,j,iold)+addtotcm*dt
plx(i,j,inew)=plx(i,j,iold)+addpx*dt
ply(i,j,inew)=ply(i,j,iold)+addpy*dt
cdif(i,j,inew)=cp(i,j,inew)-cm(i,j,inew)
```

After output, `iold` and `inew` are exchanged. The supplied settings are `k=340`, `dx=0.015`, `dt=0.000075`, and `ntime=200000000`. The periodic grid has $340\times340$ sites and side length `k*dx=5.1` in reduced coordinates.

Initial fields are generated using `ran2`:

| Field | Initial distribution |
|---|---|
| $\psi$ | `avgpsi=0.2` plus a uniform perturbation in $[-0.025,0.025]$ |
| $\hat c_+$ | `avgcp=0.15` plus a uniform perturbation in $[-0.005,0.005]$ |
| $\hat c_-$ | `avgcm=0.15` plus a uniform perturbation in $[-0.005,0.005]$ |
| $p_x,p_y$ | Each component has `avgp=0` plus a uniform perturbation in $[-0.025,0.025]$ |

At the start of each timestep, the program checks the old fields for positive microtubule fraction, nonnegative motor fractions, total fraction below one, and finite tested values. A failed check stops the run.

## Input and output

The program reads 900 filename entries from each of `fnamerhol`, `fnamerholp`, `fnamerholp5`, and `fnamerholp8`. The snapshot filenames are taken from `names`, loaded from `fnamerhol`.

Snapshots are saved at steps 10000, 20000, 40000, and subsequent doublings. The first snapshot uses `names(101)`. Each snapshot contains

| Column | Quantity |
|---|---|
| 1 | Grid index `i` |
| 2 | Grid index `j` |
| 3 | `da1`: $\psi$ |
| 4 | `da2`: $\hat c_+$ |
| 5 | `da3`: $\hat c_-$ |
| 6 | `da4`: $p_x$ |
| 7 | `da5`: $p_y$ |

`data_comp_polarity.dat` contains one parameter row in the order `vp, Dmt, Dm, anup, azetap, akap, ara`.
