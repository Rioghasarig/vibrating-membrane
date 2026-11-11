# Wave equations and forced membrane vibrations

Let \(\nabla^2=\partial_x^2+\partial_y^2\), let \(c>0\), and let \(\rho_S>0\). Complex constants and amplitudes are written in **bold**. The physical displacement is the real part of the complex solution when complex notation is used.

## 1. The two-dimensional wave equation

Consider

\[
\frac{\partial^2u}{\partial t^2}=c^2\nabla^2u.
\]

We seek separated solutions of the form

\[
u(x,y,t)=\boldsymbol{A}e^{\boldsymbol{\gamma}t}\Psi(x,y),
\]

where \(\boldsymbol{A}\) and \(\boldsymbol{\gamma}\) are complex constants. Substituting into the wave equation gives

\[
\boldsymbol{\gamma}^2\boldsymbol{A}e^{\boldsymbol{\gamma}t}\Psi
=c^2\boldsymbol{A}e^{\boldsymbol{\gamma}t}\nabla^2\Psi.
\]

For a nonzero solution, canceling the common factor and separating variables gives

\[
\frac{\nabla^2\Psi}{\Psi}=\frac{\boldsymbol{\gamma}^2}{c^2}=-k^2.
\]

Thus

\[
\boxed{\nabla^2\Psi=-k^2\Psi},\qquad
\boxed{\boldsymbol{\gamma}^2=-c^2k^2}.
\]

For any selected \(k\), we solve the spatial equation for \(\Psi\), with any applicable boundary conditions, and obtain

\[
\boxed{\boldsymbol{\gamma}=\pm ick=\pm i\omega},
\qquad \omega=ck.
\]

When \(k\) is real, these are oscillatory modes. Boundary conditions on a bounded membrane restrict which values of \(k\) are admissible.

## 2. The damped wave equation

Consider

\[
\frac{\partial^2u}{\partial t^2}
+2\beta\frac{\partial u}{\partial t}
-c^2\nabla^2u=0,
\qquad \beta\geq0.
\]

Again assume

\[
u(x,y,t)=\boldsymbol{A}e^{\boldsymbol{\gamma}t}\Psi(x,y).
\]

Substitution into the damped wave equation yields

\[
\boldsymbol{\gamma}^2\boldsymbol{A}e^{\boldsymbol{\gamma}t}\Psi
+2\beta\boldsymbol{\gamma}\boldsymbol{A}e^{\boldsymbol{\gamma}t}\Psi
-c^2\boldsymbol{A}e^{\boldsymbol{\gamma}t}\nabla^2\Psi=0.
\]

Canceling the common factor and separating variables gives

\[
\frac{\nabla^2\Psi}{\Psi}
=\frac{\boldsymbol{\gamma}^2+2\beta\boldsymbol{\gamma}}{c^2}
=-k^2.
\]

Consequently,

\[
\boxed{\nabla^2\Psi=-k^2\Psi},\qquad
\boxed{\boldsymbol{\gamma}^2+2\beta\boldsymbol{\gamma}=-c^2k^2}.
\]

With \(\omega=ck\), the temporal equation is

\[
\boldsymbol{\gamma}^2+2\beta\boldsymbol{\gamma}+\omega^2=0,
\]

so that

\[
\boxed{\boldsymbol{\gamma}=-\beta\pm\sqrt{\beta^2-\omega^2}}.
\]

In the underdamped case \(\omega^2>\beta^2\), define

\[
\omega_d=\sqrt{\omega^2-\beta^2}.
\]

Then

\[
\boxed{\boldsymbol{\gamma}=-\beta\pm i\omega_d},
\]

and the corresponding separated solutions are

\[
u(x,y,t)=\boldsymbol{A}e^{-\beta t}e^{\pm i\omega_dt}\Psi(x,y).
\]

## 3. Forced, damped circular membrane with constant pressure

Consider

\[
\frac{\partial^2u}{\partial t^2}
+2\beta\frac{\partial u}{\partial t}
-c^2\nabla^2u
=\frac{P}{\rho_S}e^{i\omega_ft},
\]

where \(P\) is constant, on the disk \(0\leq r<a\) with fixed boundary

\[
u(a,\theta,t)=0.
\]

We first find an **easy particular solution**, without imposing the boundary condition. Because the forcing is independent of position, try

\[
u_p(t)=\boldsymbol{A}_p e^{i\omega_ft}.
\]

Its Laplacian is zero. Substitution into the forced equation gives

\[
(-\omega_f^2+2i\beta\omega_f)\boldsymbol{A}_p e^{i\omega_ft}
=\frac{P}{\rho_S}e^{i\omega_ft}.
\]

Hence, for \(-\omega_f^2+2i\beta\omega_f\neq0\),

\[
\boxed{u_p(t)=\frac{P}{\rho_S(-\omega_f^2+2i\beta\omega_f)}e^{i\omega_ft}}.
\]

This particular solution does not satisfy \(u(a,\theta,t)=0\). We therefore seek

\[
\boxed{u=u_p+u_h},
\]

where \(u_h\) solves the *homogeneous* damped wave equation and cancels \(u_p\) on the boundary.

From Section 2, homogeneous separated solutions satisfy

\[
\nabla^2\Psi=-\boldsymbol{k}^2\Psi,
\qquad
\boldsymbol{\gamma}^2+2\beta\boldsymbol{\gamma}=-c^2\boldsymbol{k}^2.
\]

To cancel a boundary value proportional to \(e^{i\omega_ft}\), select \(\boldsymbol{\gamma}=i\omega_f\). The homogeneous separation relation then gives

\[
-\omega_f^2+2i\beta\omega_f=-c^2\boldsymbol{k}^2,
\qquad
\boxed{\boldsymbol{k}^2=\frac{\omega_f^2-2i\beta\omega_f}{c^2}}.
\]

For the rotationally symmetric solution, the spatial equation in polar coordinates is

\[
\Psi''(r)+\frac1r\Psi'(r)+\boldsymbol{k}^2\Psi(r)=0.
\]

The solution regular at \(r=0\) is \(J_0(\boldsymbol{k}r)\). Thus choose

\[
u_h(r,t)=\boldsymbol{B}J_0(\boldsymbol{k}r)e^{i\omega_ft}.
\]

The boundary condition gives

\[
\boldsymbol{A}_p+\boldsymbol{B}J_0(\boldsymbol{k}a)=0,
\qquad
\boldsymbol{B}=-\frac{\boldsymbol{A}_p}{J_0(\boldsymbol{k}a)}.
\]

Therefore the time-harmonic solution of the fixed-boundary problem is

\[
\boxed{
 u(r,t)=\frac{P e^{i\omega_ft}}{\rho_S(-\omega_f^2+2i\beta\omega_f)}
 \left[1-\frac{J_0(\boldsymbol{k}r)}{J_0(\boldsymbol{k}a)}\right],
 \qquad \boldsymbol{k}^2=\frac{\omega_f^2-2i\beta\omega_f}{c^2}.
}
\]

This formula assumes the displayed denominators are nonzero. For \(\beta>0\) and \(\omega_f>0\), they are nonzero. At zero forcing frequency, the static problem can instead be solved directly: \(u(r)=P(a^2-r^2)/(4\rho_Sc^2)\). Initial conditions, if prescribed, may require additional homogeneous transients.

## 4. Forced, damped circular membrane with spatially varying pressure

Now let \(P=P(r,\theta)\) vary across the disk:

\[
\frac{\partial^2u}{\partial t^2}
+2\beta\frac{\partial u}{\partial t}
-c^2\nabla^2u
=\frac{P(r,\theta)}{\rho_S}e^{i\omega_ft},
\qquad u(a,\theta,t)=0.
\]

We use the **same construction**, \(u=u_p+u_h\). A spatially constant particular solution no longer works. Instead, seek a particular solution with the forcing's time dependence,

\[
u_p(r,\theta,t)=\boldsymbol{U}_p(r,\theta)e^{i\omega_ft}.
\]

Substituting into the forced equation gives

\[
(-\omega_f^2+2i\beta\omega_f)\boldsymbol{U}_p-c^2\nabla^2\boldsymbol{U}_p
=\frac{P(r,\theta)}{\rho_S}.
\]

Using the same complex separation constant as in Section 3,

\[
\boldsymbol{k}^2=\frac{\omega_f^2-2i\beta\omega_f}{c^2},
\]

we have the inhomogeneous Helmholtz equation

\[
\boxed{(\nabla^2+\boldsymbol{k}^2)\boldsymbol{U}_p
=-\frac{P(r,\theta)}{\rho_Sc^2}}.
\]

One particular solution, constructed without enforcing the boundary condition, is the free-space Green-function integral

\[
\boxed{
\boldsymbol{U}_p(\mathbf r)
=\frac{i}{4\rho_Sc^2}
\int_{|\mathbf r'|<a}
P(\mathbf r')H_0^{(1)}\!\left(\boldsymbol{k}|\mathbf r-\mathbf r'|\right)\,dA'.
}
\]

Here \(H_0^{(1)}\) is the Hankel function of the first kind; this expression uses a nonzero \(\boldsymbol{k}\) and an appropriate outgoing-wave branch. (Other particular solutions are also possible.)

To enforce the fixed boundary, add a homogeneous solution with the same time dependence:

\[
u_h(r,\theta,t)=\boldsymbol{H}(r,\theta)e^{i\omega_ft},
\qquad
(\nabla^2+\boldsymbol{k}^2)\boldsymbol{H}=0.
\]

From the homogeneous separation equation, a general regular spatial correction may be expanded as

\[
\boldsymbol{H}(r,\theta)
=\sum_{m=-\infty}^{\infty}
\boldsymbol{B}_mJ_m(\boldsymbol{k}r)e^{im\theta}.
\]

Expand the boundary values of the already-chosen particular solution in a Fourier series:

\[
\boldsymbol{U}_p(a,\theta)
=\sum_{m=-\infty}^{\infty}\boldsymbol{p}_m e^{im\theta},
\qquad
\boldsymbol{p}_m=\frac1{2\pi}\int_0^{2\pi}
\boldsymbol{U}_p(a,\theta)e^{-im\theta}\,d\theta.
\]

Since \(u(a,\theta,t)=0\), each Fourier coefficient must satisfy

\[
\boldsymbol{p}_m+\boldsymbol{B}_mJ_m(\boldsymbol{k}a)=0,
\qquad
\boldsymbol{B}_m=-\frac{\boldsymbol{p}_m}{J_m(\boldsymbol{k}a)}.
\]

Thus the time-harmonic solution is

\[
\boxed{
 u(r,\theta,t)=e^{i\omega_ft}
 \left[
 \boldsymbol{U}_p(r,\theta)
 -\sum_{m=-\infty}^{\infty}
 \boldsymbol{p}_m\frac{J_m(\boldsymbol{k}r)}{J_m(\boldsymbol{k}a)}e^{im\theta}
 \right].
}
\]

The first term is a particular solution of the forced equation; the sum is a homogeneous solution chosen to cancel its boundary values. The formula assumes the denominators are nonzero and that the series and integral converge appropriately. For strictly positive damping and forcing frequency, \(\boldsymbol{k}^2\) is nonreal, so these denominators do not vanish. As before, an initial-value problem may additionally require homogeneous transient solutions.
