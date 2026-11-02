# Simulation of Eqs. S10a–S10c

Documentation for `S10_commented_code.f`. 

## Fields and notation

The program evolves five fields on a two-dimensional periodic square grid.

| Code | Quantity |
|---|---|
| `x` | Microtubule area fraction, $\psi$ |
| `cp`, `cm` | Motor area fractions, $\hat c_+$ and $\hat c_-$ |
| `plx`, `ply` | Polarization components, $p_x$ and $p_y$ |
| `cdif` | Motor area-fraction difference, $\hat c_+-\hat c_-$ |
| `iold`, `inew` | Current and updated time buffers |
| `p(i)`, `m(i)` | Next and previous grid indices, with periodic wrapping |

Below, all derivatives are with respect to the dimensionless coordinates of S10. Define $u=1-\hat c_+-\hat c_--\psi$ and $g=\tilde\nu/\tilde v_0$.

## Parameters

Values marked “variable” are varied between simulation runs. The derived coefficient `azetap` is recalculated as `azeta*arho/lp**2`.

| Code | Definition or meaning | Value |
|---|---|---:|
| `Dr` | Rotational rate used for time scaling | 0.1 |
| `lp`, `Le` | Length scale and microtubule length $L$ | 5 |
| `d` | Radius $d$ | 0.02 |
| `lsq` | $2L$ | 10 |
| `v0` | Dimensional motor speed $v_0$ | 0.25 |
| `alpha`, `bet` | $\alpha$, $\beta$ | between 100 and 1000, 1 |
| `Dpsi`, `Dc` | Dimensional microtubule and motor diffusivities | 0.1, 1 |
| `arho` | $\rho_0=1/(\pi d^2)$ | 795.774715 |
| `anu` | Effective dimensional coupling $\nu$ | 0.007 |
| `azeta` | Dimensional motor gradient coefficient | variable  |
| `akap` | $\tilde\kappa/\psi_0$ | 0.00005 |
| `Dmt` | $\tilde D_{mt}=D_\psi/(D_rL^2)$ | 0.04 |
| `Dm` | $\tilde D=D_c/(D_rL^2)$ | 0.4 |
| `vp` | $\tilde v_0=v_0/(D_rL)$ | 0.5 |
| `anup` | $g=\nu\rho_0/L=\tilde\nu/\tilde v_0$ | 1.114085 |
| `azetap` | $2\tilde\zeta/(\hat c_0^++\hat c_0^-)$ | variable  |
| `ara` | $1/\mathfrak R=D_c\ln(L/d)/(\pi d v_0)$ | 351.507119 |

The SI definition $\tilde\nu=\nu\rho_0\tilde v_0/L$ gives $\tilde\nu\simeq0.55704$. The code uses **`anup` = $\tilde\nu/\tilde v_0$**, not $\tilde\nu$.

The array `vmfield` stores the local value of $-\tilde v_{mt}$ from the SI. At each grid point, `vm=vmfield(i,j)` retrieves this value:

$$
\texttt{vmfield}=\texttt{vm}
=\tilde v_0\frac{\hat c_+-\hat c_-}{\hat c_++\hat c_-}
=-\tilde v_{mt}.
$$

It is set to zero when the total motor density is zero.

## S10a: motor-density equations

| Code | Calculated quantity |
|---|---|
| `conpxp`, `conpyp` | Components of $\mathfrak R^{-1}(-\tilde v_0+\texttt{vm})\psi\hat c_+\mathbf p$ |
| `conmxm`, `conmym` | Components of $\mathfrak R^{-1}(\tilde v_0+\texttt{vm})\psi\hat c_-\mathbf p$ |
| `acpxx+acpyy`, `acmxx+acmyy` | $\Delta\hat c_+$ and $\Delta\hat c_-$ |
| `agdcpxx`, `agdcpyy` | Components of $(\hat c_+/u)\nabla\hat c_+$ |
| `agdcmxx`, `agdcmyy` | Components of $(\hat c_-/u)\nabla\hat c_-$ |
| `bpsx`, `bpsy` | Components of $\psi\mathbf p$ |
| `amucp`, `amucm` | $\mu_\pm^{\mathrm{grad}}=-\texttt{azetap} \Delta\hat c_\pm\pm g \pi d/(2L)\nabla\cdot(\psi\mathbf p)$ |
| `ajcpx`, `ajcpy` | Components of $\hat c_+\nabla\mu_+^{\mathrm{grad}}$ |
| `ajcmx`, `ajcmy` | Components of $\hat c_-\nabla\mu_-^{\mathrm{grad}}$ |
| `addtotcp`, `addtotcm` | Complete right-hand sides of S10a |

The program takes the divergence of each stored vector. Diffusion, self-gradient crowding, and gradient-coupling contributions are multiplied by `Dm`; active transport is added separately. These stored vectors appear with a **positive divergence** in the equations.

## S10b: microtubule-density equation

Define

$$
A=|\mathbf p|^2\left(1+\frac{\alpha}{2}|\mathbf p|^2\right),\qquad
B=\mathbf p\cdot\nabla(\hat c_+-\hat c_-),\qquad
H=\psi(\nabla\cdot\mathbf p)^2.
$$

| Code | Calculated quantity |
|---|---|
| `agrp2`, `dpcdf`, `aggp2` | $A$, $B$, $H$ |
| `conpx`, `conpy` | Components of $\texttt{vm}\psi\mathbf p$ |
| `addcompsi` | $\nabla\cdot(\texttt{vm}\psi\mathbf p)$ |
| `agdpx`, `agdpy` | Components of $(2L/\pi d)(\psi/u)\nabla\psi$ |
| `tlapsi` | $\Delta\psi+(2L/\pi d)\nabla\cdot[(\psi/u)\nabla\psi]$ |
| `divpsx`, `divpsy` | Components of $\psi\nabla A$ |
| `tbetpsi` | $(\beta/2)\nabla\cdot(\psi\nabla A)$ |
| `dpcdfpsi1`, `dpcdfpsi2` | Components of $\psi\nabla B$ |
| `tnupsict` | $-g\nabla\cdot(\psi\nabla B)$ |
| `aggpx2s`, `aggpy2s` | Components of $\psi\nabla H$ |
| `tkapsi` | $\beta \texttt{akap} \nabla\cdot(\psi\nabla H)$ |

`addpsitotal` is `Dmt` times the sum of `tlapsi`, `tbetpsi`, `tnupsict`, and `tkapsi`. The complete density right-hand side is `addpsitotal+addcompsi`.

## S10c: polarization equation

`polOP` is $|\mathbf p|^2$, and `pspx` is the scalar $\psi^2\nabla\cdot\mathbf p$.

| Code contribution, after division by `x` | Term in the polarization equation |
|---|---|
| `tpbetax`, `tpbetay` | $-\mathbf p(1+\alpha \vert \mathbf p \vert^2)$ |
| `tnux/bet`, `tnuy/bet` | $(g/\beta)\nabla(\hat c_+-\hat c_-)$ |
| `tpkapx`, `tpkapy` | $(\texttt{akap}/\psi)\nabla[\psi^2\nabla\cdot\mathbf p]$ |
| `-advpx*x`, `-advpy*x` | $\texttt{vm}(\mathbf p\cdot\nabla)\mathbf p=-\tilde v_{mt}(\mathbf p\cdot\nabla)\mathbf p$ |

The elastic contribution is the conservative-form discretization of

$$
\frac{\texttt{akap}}{\psi}\nabla[\psi^2\nabla\cdot\mathbf p]
=\tilde\kappa\left[
2(\nabla\cdot\mathbf p)\nabla\left(\frac{\psi}{\psi_0}\right)
+\frac{\psi}{\psi_0}\nabla(\nabla\cdot\mathbf p)\right],
\qquad \tilde\kappa=\texttt{akap} \psi_0.
$$

`addpx` and `addpy` are the complete polarization right-hand sides.

## Numerical steps

The number of grid points `k`, grid spacing `dx`, and timestep `dt` are selected for each run. First derivatives use centered differences; Laplacians use the nearest-neighbor five-point stencil.

1. **Initialize the fields.** Set the mean microtubule area fraction to 0.2, select the mean motor area fractions, and initialize each polarization component around zero. Add uniform random perturbations using `ran2`.

2. **Check the fields and calculate microtubule transport.** At each timestep, check the current fields and denominators. Calculate $-\tilde v_{mt}=\tilde v_0(\hat c_+-\hat c_-)/(\hat c_++\hat c_-)$ and store it in `vmfield`. Set this value to zero when the total motor area fraction is zero.

3. **First grid sweep: calculate intermediate quantities.** At every grid point, calculate $A$, $B$, and $H$ defined above, together with $\psi\mathbf p$ and $\psi^2\nabla\cdot\mathbf p$. Calculate the active transport vectors and the crowding vectors $(2L/\pi d)(\psi/u)\nabla\psi$ and $(\hat c_\pm/u)\nabla\hat c_\pm$.

4. **Second grid sweep: differentiate the intermediate quantities.** Calculate $\psi\nabla A$, $\psi\nabla B$, and $\nabla H$ from the completed arrays. Calculate the motor Laplacians and $\nabla\cdot(\psi\mathbf p)$, then construct the motor potentials `amucp` and `amucm` defined above.

5. **Third grid sweep: prepare vectors for the density equations.** Multiply $\nabla H$ by $\psi$ and store its components in `aggpx2s` and `aggpy2s`. This vector supplies the elastic contribution to the microtubule-density equation through $\tilde D_{mt}\beta\,\texttt{akap}\,\nabla\cdot(\psi\nabla H)$. Calculate the gradients of `amucp` and `amucm` and multiply them by the respective motor area fractions to form $\hat c_\pm\nabla\mu_\pm^{\mathrm{grad}}$.

6. **Calculate the time derivatives and advance the fields.** Take the required divergences and assemble the complete right-hand sides of S10a–S10c. Use forward Euler, $f^{n+1}=f^n+\texttt{dt}\,R_f^n$, to write the updated fields into `inew`.

7. **Complete the timestep.** Update `cdif`, write a snapshot when scheduled, and swap `iold` and `inew`.

Each grid sweep is completed over the whole grid before the next begins, so derivatives use completed neighboring intermediate values. All intermediate quantities and right-hand sides are evaluated from the fields at the current time level `iold`.

`ieee_is_finite` detects non-finite values. The program stops if the tested fields are non-finite, $\psi\leq0$, either motor density is negative, or $\psi+\hat c_++\hat c_-\geq1$. It does not clip the densities.

## Input and output

The program reads 900 filename entries from each of `fnamerhol`, `fnamerholp`, `fnamerholp5`, and `fnamerholp8`. Snapshot filenames come from `fnamerhol`, beginning with entry 101; names should fit the eight-character strings declared in the code.

Snapshots are saved at $steps*D_r=$ 1, 2, 4, and successive doublings within the run. Each row contains

```text
i  j  psi  c_plus  c_minus  p_x  p_y
```

`data_comp_polarity.dat` records `vp Dmt Dm anup azetap akap ara` at the first step. `data.dat` is opened but receives no data in this version.
