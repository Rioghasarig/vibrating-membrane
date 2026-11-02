# Damped circular membrane: nonuniform density

## 1. Damped wave equation

For a circular membrane of radius 1, nonuniform positive density $\rho(x,y)$, constant tension $T$, and constant damping coefficient $\gamma$, the wave equation is

$$
\boxed{\rho(x,y)z_{tt}+\gamma z_t=T\Delta z+\operatorname{Re}\{V(x,y)e^{i\omega t}\}.}
$$

The fixed boundary condition is

$$
z(x,y,t)=0,\qquad x^2+y^2=1.
$$

## 2. Eigenfunction expansion

Let $\phi_n(x,y)$ satisfy

$$
\boxed{-T\Delta\phi_n=\rho(x,y)\omega_n^2\phi_n}
$$

with $\phi_n=0$ on the boundary. Normalize using **weighted orthogonality**:

$$
\boxed{\int_D\rho(x,y)\phi_m\phi_n\,dx\,dy=\delta_{mn},}
$$

where $D=\{(x,y):x^2+y^2<1\}$. Expand

$$
\boxed{z(x,y,t)=\sum_nq_n(t)\phi_n(x,y).}
$$

## 3. Substitute into the wave equation

$$
z_t=\sum_nq_n'\phi_n
$$

$$
z_{tt}=\sum_nq_n''\phi_n
$$

$$
\begin{aligned}
\Delta z&=\sum_nq_n\Delta\phi_n\\
&=-\sum_n\frac{\rho(x,y)\omega_n^2}{T}q_n\phi_n.
\end{aligned}
$$

Substituting,

$$
\rho\sum_nq_n''\phi_n+\gamma\sum_nq_n'\phi_n
=-\rho\sum_n\omega_n^2q_n\phi_n+\operatorname{Re}\{Ve^{i\omega t}\}.
$$

Therefore,

$$
\boxed{\sum_n\left[\rho(x,y)(q_n''+\omega_n^2q_n)+\gamma q_n'\right]\phi_n
=\operatorname{Re}\{Ve^{i\omega t}\}.}
$$

## 4. Apply weighted orthogonality

Multiply by $\phi_m$ and integrate over $D$:

$$
\begin{aligned}
&\sum_n(q_n''+\omega_n^2q_n)\int_D\rho\phi_n\phi_m\,dx\,dy\\
&\quad+\gamma\sum_nq_n'\int_D\phi_n\phi_m\,dx\,dy\\
&=\operatorname{Re}\left\{e^{i\omega t}\int_DV\phi_m\,dx\,dy\right\}.
\end{aligned}
$$

Using $\int_D\rho\phi_n\phi_m\,dx\,dy=\delta_{nm}$,

$$
\boxed{q_m''+\omega_m^2q_m+\gamma\sum_nC_{mn}q_n'
=\operatorname{Re}\{V_me^{i\omega t}\}}
$$

where

$$
\boxed{C_{mn}=\int_D\phi_m(x,y)\phi_n(x,y)\,dx\,dy}
$$

and

$$
\boxed{V_m=\int_DV(x,y)\phi_m(x,y)\,dx\,dy.}
$$

Unlike the uniform-density case, constant damping generally couples different modes.

## 5. Solve the coupled oscillators

Assume steady-state particular solutions

$$
q_n^{(p)}(t)=\operatorname{Re}\{Q_ne^{i\omega t}\}.
$$

Then

$$
(q_n^{(p)})'=\operatorname{Re}\{i\omega Q_ne^{i\omega t}\}
$$

$$
(q_n^{(p)})''=\operatorname{Re}\{-\omega^2Q_ne^{i\omega t}\}.
$$

Substituting into the modal equation,

$$
\operatorname{Re}\left\{\left[(\omega_m^2-\omega^2)Q_m+i\gamma\omega\sum_nC_{mn}Q_n\right]e^{i\omega t}\right\}
=\operatorname{Re}\{V_me^{i\omega t}\}.
$$

Equating complex amplitudes,

$$
\boxed{(\omega_m^2-\omega^2)Q_m+i\gamma\omega\sum_nC_{mn}Q_n=V_m.}
$$

Define

$$
(\Omega^2)_{mn}=\omega_n^2\delta_{mn},\qquad
C_{mn}=\int_D\phi_m\phi_n\,dx\,dy.
$$

Then

$$
\left(\Omega^2-\omega^2I+i\gamma\omega C\right)\mathbf Q=\mathbf V.
$$

Therefore,

$$
\boxed{\mathbf Q=\left(\Omega^2-\omega^2I+i\gamma\omega C\right)^{-1}\mathbf V.}
$$

This inverse can be computed numerically after truncating the expansion to finitely many modes.

## 6. Reconstruct the displacement

$$
\begin{aligned}
z_{\mathrm{ss}}(x,y,t)
&=\sum_nq_n^{(p)}(t)\phi_n(x,y)\\
&=\sum_n\operatorname{Re}\{Q_ne^{i\omega t}\}\phi_n(x,y)\\
&=\operatorname{Re}\left\{e^{i\omega t}\sum_nQ_n\phi_n(x,y)\right\}.
\end{aligned}
$$

Thus,

$$
\boxed{z_{\mathrm{ss}}(x,y,t)=\operatorname{Re}\left\{e^{i\omega t}\sum_nQ_n\phi_n(x,y)\right\},}
$$

where the $Q_n$ are determined by the coupled linear system in Step 5.

This is the **steady-state solution**. The full solution also includes transient terms determined by the initial conditions.
