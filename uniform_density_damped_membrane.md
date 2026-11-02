# Damped circular membrane: uniform density

## 1. Damped wave equation

For a circular membrane of radius 1, uniform density $\rho$, constant tension $T$, and damping coefficient $\gamma$, the wave equation is

$$
\boxed{\rho z_{tt}+\gamma z_t=T\Delta z+\operatorname{Re}\{V(x,y)e^{i\omega t}\}.}
$$

The fixed boundary condition is

$$
z(x,y,t)=0,\qquad x^2+y^2=1.
$$

## 2. Eigenfunction expansion

Let $\phi_n(x,y)$ satisfy

$$
\boxed{-T\Delta\phi_n=\rho\omega_n^2\phi_n}
$$

with $\phi_n=0$ on the boundary. Normalize them so that

$$
\int_D\phi_m\phi_n\,dx\,dy=\delta_{mn},
$$

where $D=\{(x,y):x^2+y^2<1\}$. Expand

$$
\boxed{z(x,y,t)=\sum_n q_n(t)\phi_n(x,y).}
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
&=-\sum_n\frac{\rho\omega_n^2}{T}q_n\phi_n.
\end{aligned}
$$

Substituting,

$$
\rho\sum_nq_n''\phi_n+\gamma\sum_nq_n'\phi_n
=-\rho\sum_n\omega_n^2q_n\phi_n+\operatorname{Re}\{Ve^{i\omega t}\}.
$$

Therefore,

$$
\boxed{\sum_n\left(\rho q_n''+\gamma q_n'+\rho\omega_n^2q_n\right)\phi_n
=\operatorname{Re}\{Ve^{i\omega t}\}.}
$$

## 4. Apply orthogonality

Multiply by $\phi_m$ and integrate over $D$:

$$
\begin{aligned}
&\sum_n\left(\rho q_n''+\gamma q_n'+\rho\omega_n^2q_n\right)
\int_D\phi_n\phi_m\,dx\,dy\\
&\qquad=\operatorname{Re}\left\{e^{i\omega t}\int_DV\phi_m\,dx\,dy\right\}.
\end{aligned}
$$

Using $\int_D\phi_n\phi_m\,dx\,dy=\delta_{nm}$,

$$
\boxed{\rho q_m''+\gamma q_m'+\rho\omega_m^2q_m
=\operatorname{Re}\{V_me^{i\omega t}\}}
$$

where **the driving-force coefficient** is

$$
\boxed{V_m=\int_DV(x,y)\phi_m(x,y)\,dx\,dy.}
$$

## 5. Solve each oscillator

Assume a steady-state particular solution

$$
q_m^{(p)}(t)=\operatorname{Re}\{Q_me^{i\omega t}\}.
$$

Then

$$
(q_m^{(p)})'=\operatorname{Re}\{i\omega Q_me^{i\omega t}\}
$$

$$
(q_m^{(p)})''=\operatorname{Re}\{-\omega^2Q_me^{i\omega t}\}.
$$

Substituting into the modal equation,

$$
\operatorname{Re}\left\{(-\rho\omega^2+i\gamma\omega+\rho\omega_m^2)Q_me^{i\omega t}\right\}
=\operatorname{Re}\{V_me^{i\omega t}\}.
$$

Equating complex amplitudes,

$$
\left[\rho(\omega_m^2-\omega^2)+i\gamma\omega\right]Q_m=V_m.
$$

Therefore,

$$
\boxed{Q_m=\frac{V_m}{\rho(\omega_m^2-\omega^2)+i\gamma\omega}.}
$$

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
\boxed{z_{\mathrm{ss}}(x,y,t)=\operatorname{Re}\left\{e^{i\omega t}\sum_n\frac{V_n\phi_n(x,y)}{\rho(\omega_n^2-\omega^2)+i\gamma\omega}\right\}.}
$$

This is the **steady-state solution**. The full solution also includes transient terms determined by the initial conditions.
