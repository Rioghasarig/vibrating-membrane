# Equations used by the simulator

This is every equation the code evaluates, in the order it runs, with links to the lines that do it. Sections 1, 2 and 4 follow [`membrane-fem.md`](../membrane-fem.md). Section 3 derives the damping and force terms by virtual displacements. The force term isn't in the note, and the damping term matches the note's step 4. Sections 5 and 6 give the force formulas and the time stepping, which the note doesn't contain.

## 1. The equation being solved

$$
M \ddot{\mathbf q} + C \dot{\mathbf q} + K \mathbf q = \mathbf F(t)
$$

$\mathbf q$ holds the displacements of the $N$ interior nodes. Boundary nodes are always zero.

Compared with the note:
- $C$ is the note's damping matrix from step 4 (sections 3.8, 3.9 and 4).
- $\mathbf F(t)$ is not in the note (sections 3.7 and 5).
- With no drives, $\mathbf F = 0$, and this is exactly the note's $M\ddot{\mathbf q} + C\dot{\mathbf q} + K\mathbf q = 0$.

## 2. Mass and stiffness matrices

[fem.py:56-79](membrane/fem.py#L56-L79)

For a triangle with vertices $(x_1,y_1), (x_2,y_2), (x_3,y_3)$, using cyclic indices $(i,j,k)$:

$$
b_i = y_j - y_k, \qquad c_i = x_k - x_j, \qquad 2A = b_1 c_2 - b_2 c_1
$$

The hat-function gradients are constant on the triangle:

$$
\nabla \varphi_i = \frac{1}{2A}(b_i, c_i)
$$

Element stiffness matrix:

$$
K^e_{ij} = T \int_{\Omega_e} \nabla\varphi_i \cdot \nabla\varphi_j \, dA = \frac{T}{4A}\left(b_i b_j + c_i c_j\right)
$$

Element mass matrix:

$$
M^e = \sigma \frac{A}{12}
\begin{pmatrix} 2 & 1 & 1 \\ 1 & 2 & 1 \\ 1 & 1 & 2 \end{pmatrix}
$$

The element matrices are summed into matrices over all nodes, and the boundary rows and columns are then dropped. This gives the $N \times N$ matrices $M$ and $K$, as the note describes.

The code also keeps the mass matrix with $\sigma = 1$ and every column, which sections 3.7 and 5 use:

$$
B_{ij} = \int \varphi_i \varphi_j \, dA, \qquad i \text{ interior}, \; j \text{ any node}
$$

## 3. Deriving the force and damping terms by virtual displacements

The note reaches its equations through Lagrange's equations, with damping added as a generalized force found by virtual work. This section uses the principle of virtual displacements for the whole equation, which makes external forces easy to add. This section derives the same equations that way, then uses it to get the driving force (section 5) and the damping matrix (section 4). Section 3.9 checks that the damping matrix agrees with the note's.

### 3.1 Equation of motion of a small piece of membrane

Take a small patch of membrane with area $dA$. Three vertical forces act on it.

- **Tension.** Around the patch's edge, the tension pulls with force $T$ per unit length. For small slopes, its vertical component along an edge with outward normal $\mathbf n$ is $T\,\nabla u \cdot \mathbf n$ per unit length. By the divergence theorem, the total over the edge is
  $$
  \oint T\,\nabla u \cdot \mathbf n \, ds = \int T\,\Delta u \, dA \;\approx\; T\,\Delta u \, dA
  $$
- **External force.** Any other force is written $f(x,y,t)$ per unit area, so it contributes $f\,dA$. Section 3.7 uses $f = p$ for pressure, and section 3.8 uses $f = -c\,u_t$ for damping.
- **Mass times acceleration.** The patch has mass $\sigma\,dA$ and acceleration $u_{tt}$.

Newton's second law, divided by $dA$, gives

$$
\sigma\, u_{tt} = T\,\Delta u + f
$$

This holds at every point of the membrane. With $f = 0$, it is the Euler–Lagrange equation of the note's Lagrangian.

### 3.2 The virtual displacement statement

A virtual displacement $\delta u(x,y)$ is any small shape the membrane could be imagined to move by, with time frozen. It is admissible if it respects the constraint, which here means $\delta u = 0$ on the rim.

Move every force from section 3.1 to one side, multiply by $\delta u$, and integrate over the membrane:

$$
\int_D \left( \sigma\, u_{tt} - T\,\Delta u - f \right) \delta u \, dA = 0 \qquad \text{for every admissible } \delta u
$$

The bracket is zero everywhere, so this certainly holds. The converse is also true: if the integral is zero for every admissible $\delta u$, the bracket must be zero everywhere. Otherwise, a $\delta u$ concentrated where the bracket is nonzero would make the integral nonzero. So this statement is equivalent to the equation of motion.

Each term is a virtual work: force per unit area times the virtual displacement, summed over the membrane. The inertia term $-\sigma u_{tt}$ is treated as a force, which is d'Alembert's principle.

### 3.3 Integrate the tension term by parts

The tension term contains $\Delta u$, a second derivative. The hat functions have kinks along the triangle edges, so their second derivatives don't exist as ordinary functions. The fix is to move one derivative onto $\delta u$.

Green's first identity states

$$
\int_D \delta u\, \Delta u \, dA = \oint_{\partial D} \delta u\, \nabla u \cdot \mathbf n \, ds - \int_D \nabla \delta u \cdot \nabla u \, dA
$$

The boundary integral is zero, because $\delta u = 0$ on the rim. Substituting the rest into section 3.2 and rearranging gives

$$
\int_D \sigma\, u_{tt}\, \delta u \, dA + \int_D T\, \nabla u \cdot \nabla \delta u \, dA = \int_D f\, \delta u \, dA
\qquad \text{for every admissible } \delta u \tag{$*$}
$$

The second term is the change in potential energy $\int \tfrac{T}{2}|\nabla u|^2\,dA$ caused by $\delta u$. This is the potential energy in the note's Lagrangian.

### 3.4 Restrict both $u$ and $\delta u$ to the hat functions

As in the note, approximate the displacement by the hat functions of the $N$ interior nodes:

$$
u(x,y,t) \approx \sum_{j=1}^{N} q_j(t)\, \varphi_j(x,y)
$$

Only the coefficients $q_j$ depend on time, and only the hat functions depend on position. So the derivatives needed in $(*)$ are

$$
u_{tt} = \sum_{j=1}^{N} \ddot q_j\, \varphi_j,
\qquad
\nabla u = \sum_{j=1}^{N} q_j\, \nabla \varphi_j
$$

Virtual displacements are chosen from the same family:

$$
\delta u = \sum_{i=1}^{N} \delta q_i\, \varphi_i
$$

The $\delta q_i$ are arbitrary numbers. Every interior hat function is zero on the rim, so every such $\delta u$ is automatically admissible.

From here on, the integrals are over the mesh, a polygon inscribed in the circle, rather than the exact disk.

### 3.5 Substitute into each term

Each term of $(*)$ is now a double sum. The $q$'s and $\delta q$'s don't depend on position, so they come outside the integrals.

**Inertia term:**

$$
\int \sigma \left( \sum_j \ddot q_j\, \varphi_j \right) \left( \sum_i \delta q_i\, \varphi_i \right) dA
= \sum_i \delta q_i \sum_j \left( \sigma \int \varphi_i \varphi_j \, dA \right) \ddot q_j
= \sum_i \delta q_i \sum_j M_{ij}\, \ddot q_j
$$

$M_{ij} = \sigma \int \varphi_i \varphi_j \, dA$ is the note's mass matrix. Splitting the integral into one integral per triangle gives the note's assembly $M = \sum_e \hat M^e$.

**Tension term:**

$$
\int T \left( \sum_j q_j\, \nabla\varphi_j \right) \cdot \left( \sum_i \delta q_i\, \nabla\varphi_i \right) dA
= \sum_i \delta q_i \sum_j \left( T \int \nabla\varphi_i \cdot \nabla\varphi_j \, dA \right) q_j
= \sum_i \delta q_i \sum_j K_{ij}\, q_j
$$

$K_{ij}$ is the note's stiffness matrix.

**External force term:**

$$
\int f \left( \sum_i \delta q_i\, \varphi_i \right) dA
= \sum_i \delta q_i \underbrace{\int f\, \varphi_i \, dA}_{Q_i}
$$

Putting the three terms back into $(*)$ and moving everything to one side:

$$
\sum_{i=1}^{N} \delta q_i \left[ \sum_j M_{ij}\, \ddot q_j + \sum_j K_{ij}\, q_j - Q_i \right] = 0
\qquad \text{for every choice of } \delta q_1, \dots, \delta q_N
$$

### 3.6 One equation per node

Choose $\delta q_k = 1$ and every other $\delta q_i = 0$. This means $\delta u = \varphi_k$: node $k$ moves up by one unit and all other nodes stay still. Every term of the sum vanishes except the $k$-th, so its bracket must be zero:

$$
\sum_j M_{kj}\, \ddot q_j + \sum_j K_{kj}\, q_j = Q_k
$$

Doing this for each $k = 1, \dots, N$ gives the full system:

$$
M \ddot{\mathbf q} + K \mathbf q = \mathbf Q,
\qquad
Q_i = \int f\, \varphi_i \, dA
$$

With $f = 0$, this is the note's equation. $Q_i$ is the generalized force on node $i$: the virtual work $f$ does when node $i$ alone moves by one unit. All that remains is to evaluate $Q_i$ for each kind of external force.

### 3.7 Pressure and point forces

**Pressure.** Set $f = p(x,y,t)$:

$$
F_i = \int p\, \varphi_i \, dA
$$

To compute this, the code interpolates $p$ with the hat functions of every node, including those on the rim:

$$
p(x,y,t) \approx \sum_{j \in \text{all nodes}} p_j(t)\, \varphi_j(x,y),
\qquad p_j(t) = p(x_j, y_j, t)
$$

Substituting:

$$
F_i = \sum_{j \in \text{all nodes}} \left( \int \varphi_i \varphi_j \, dA \right) p_j = \sum_j B_{ij}\, p_j
$$

This is why $B$ has rim columns but no rim rows. The row index $i$ comes from $\delta u$, which is zero on the rim. The column index $j$ comes from $p$, which is not.

**Point force.** A force $f(t)$ in newtons at $(x_0, y_0)$ is a force per unit area concentrated at one point. Write it with the Dirac delta $\delta_D$, which is unrelated to the virtual $\delta$:

$$
f = f(t)\, \delta_D(x - x_0)\, \delta_D(y - y_0)
$$

The integral picks out the value of $\varphi_i$ at that point:

$$
F_i = \int f(t)\, \delta_D(x - x_0)\, \delta_D(y - y_0)\, \varphi_i \, dA = f(t)\, \varphi_i(x_0, y_0)
$$

### 3.8 Distributed damping

Model damping as a viscous force per unit area that opposes the local velocity, with coefficient $c(x,y)$ in N·s/m³:

$$
f = -c\, u_t,
\qquad
u_t = \sum_{j=1}^{N} \dot q_j\, \varphi_j
$$

Its generalized force is found the same way as the inertia term in section 3.5:

$$
Q_i = -\int c \left( \sum_j \dot q_j\, \varphi_j \right) \varphi_i \, dA
= -\sum_j \left( \int c\, \varphi_i \varphi_j \, dA \right) \dot q_j
= -\sum_j C_{ij}\, \dot q_j
$$

With both pressure and damping acting, $f = p - c\,u_t$. The damping part moves to the left-hand side:

$$
M \ddot{\mathbf q} + C \dot{\mathbf q} + K \mathbf q = \mathbf F,
\qquad
C_{ij} = \int c\, \varphi_i \varphi_j \, dA
$$

$C$ is symmetric and has the same sparsity pattern as $M$. It is not diagonal. For uniform $c$, $C = (c/\sigma)\, M$ exactly. In that case each mode of the membrane obeys $\ddot q + (c/\sigma)\dot q + \omega^2 q = 0$, so its amplitude decays as $e^{-c t / 2\sigma}$.

### 3.9 Agreement with the note's step 4

Step 4 of the note derives damping with Lagrange's equations rather than for the whole equation of motion at once. Only the damping force is treated by virtual work. Its virtual work is

$$
\delta W = -\iint_D c\, u_t\, \delta u \, dA = -\sum_i \delta q_i \sum_j C_{ij}\, \dot q_j
$$

The coefficient of $\delta q_i$ is the generalized force $Q_i = -\sum_j C_{ij}\, \dot q_j$, which is the same as in section 3.8. The note adds it to the right-hand side of Lagrange's equations, whose left-hand side gives $M\ddot{\mathbf q} + K\mathbf q$, so both routes give

$$
M \ddot{\mathbf q} + C \dot{\mathbf q} + K \mathbf q = \mathbf F
$$

with $\mathbf F = 0$ in the note. Either way, $C$ is not diagonal: $C_{ij}$ is nonzero whenever nodes $i$ and $j$ share a triangle. The code uses $C$ exactly as derived, with no lumping onto the diagonal.

## 4. Damping matrix (note, step 4)

[fem.py:86-91](membrane/fem.py#L86-L91) and [fem.py:93-99](membrane/fem.py#L93-L99)

This is the matrix $C_{ij} = \int c\, \varphi_i \varphi_j \, dA$ from section 3.8 and the note's step 4. The user gives the damping coefficient per unit area, $c(x,y)$, in N·s/m³: a number, a function, or nothing for no damping.

The integral is split into triangles, in the same way as $M$. The code takes $c$ as constant on each triangle, equal to its value $c_e$ at the triangle's centroid. The element matrix is then the element mass matrix with $\sigma$ replaced by $c_e$:

$$
C^e = c_e \frac{A}{12}
\begin{pmatrix} 2 & 1 & 1 \\ 1 & 2 & 1 \\ 1 & 1 & 2 \end{pmatrix}
$$

The element matrices are summed into a matrix over all nodes, and the boundary rows and columns are dropped, exactly as for $M$ and $K$.

- **Uniform $c$:** the centroid value is exact, so $C = (c/\sigma)\, M$ with no approximation.
- **Non-uniform $c$:** treating $c$ as constant on each triangle is the only approximation in $C$. Its error shrinks as the mesh is refined. A sharp edge in $c$, like the pad in the `damped-patch` scenario, follows the triangle edges rather than the exact curve.

## 5. Driving force (not from the note)

[loads.py](membrane/loads.py)

These formulas are derived in section 3.7. The force comes from the work a pressure $p$ does on the membrane:

$$
F_i(t) = \int p(x,y,t)\, \varphi_i \, dA
$$

With $p$ interpolated onto all mesh nodes as $\mathbf p(t)$, this becomes $\mathbf F(t) = B\, \mathbf p(t)$. Each drive is one of these:

| Drive | Force | Code |
|---|---|---|
| Pressure field | $\mathbf F(t) = B\,\mathbf p(t)$, with $p_j(t) = p(x_j, y_j, t)$ | [loads.py:129-136](membrane/loads.py#L129-L136) |
| Plane wave | $p_j(t) = s\!\left(t - \dfrac{\sin\theta\,(x_j\cos\alpha + y_j\sin\alpha)}{c_\text{air}}\right)$, then as above | [loads.py:106-120](membrane/loads.py#L106-L120) |
| Point force | $F_i(t) = f(t)\,\varphi_i(x_0,y_0)$ | [loads.py:123-126](membrane/loads.py#L123-L126) |

For the plane wave:
- $s(t)$ is the pressure waveform.
- $\theta$ is the angle from the membrane normal.
- $\alpha$ is the direction the wave travels across the membrane.
- $c_\text{air} = 343$ m/s.

When there are several drives, their forces are added.

## 6. Time stepping (not from the note)

[solver.py:64-126](membrane/solver.py#L64-L126)

The note stops at the differential equation. To step through time, the code uses the Newmark average-acceleration method. Write $\mathbf v = \dot{\mathbf q}$ and $\mathbf a = \ddot{\mathbf q}$, and use subscript $n$ for time $t_n = n\,\Delta t$.

### Time step

[solver.py:64-67](membrane/solver.py#L64-L67)

The default step is half the time a wave takes to cross one element:

$$
\Delta t = \frac{h}{2c}, \qquad h = \frac{a}{n_\text{rings}}, \qquad c = \sqrt{T/\sigma}
$$

The step is then shrunk slightly so that a whole number of steps fills the duration.

### Starting acceleration

[solver.py:97-99](membrane/solver.py#L97-L99)

From the initial displacement $\mathbf q_0$ and velocity $\mathbf v_0$:

$$
M \mathbf a_0 = \mathbf F(0) - C \mathbf v_0 - K \mathbf q_0
$$

### Update rule

The method assumes the acceleration over a step is the average of its values at the two ends:

$$
\mathbf q_{n+1} = \mathbf q_n + \Delta t\, \mathbf v_n + \frac{\Delta t^2}{4}\left(\mathbf a_n + \mathbf a_{n+1}\right)
$$

$$
\mathbf v_{n+1} = \mathbf v_n + \frac{\Delta t}{2}\left(\mathbf a_n + \mathbf a_{n+1}\right)
$$

$\mathbf a_{n+1}$ is unknown, so the code splits each update into a known part (the prediction) and a part that depends on $\mathbf a_{n+1}$.

**Predict** ([solver.py:107-108](membrane/solver.py#L107-L108)):

$$
\mathbf q^* = \mathbf q_n + \Delta t\, \mathbf v_n + \frac{\Delta t^2}{4}\mathbf a_n,
\qquad
\mathbf v^* = \mathbf v_n + \frac{\Delta t}{2}\mathbf a_n
$$

The updates then read:

$$
\mathbf q_{n+1} = \mathbf q^* + \frac{\Delta t^2}{4}\mathbf a_{n+1},
\qquad
\mathbf v_{n+1} = \mathbf v^* + \frac{\Delta t}{2}\mathbf a_{n+1}
$$

**Solve** ([solver.py:109](membrane/solver.py#L109)). The equation of motion must hold at the next time:

$$
M \mathbf a_{n+1} + C \mathbf v_{n+1} + K \mathbf q_{n+1} = \mathbf F(t_{n+1})
$$

Substituting the updates and moving the known terms to the right gives one linear system for $\mathbf a_{n+1}$:

$$
\underbrace{\left(M + \frac{\Delta t}{2}C + \frac{\Delta t^2}{4}K\right)}_{S} \mathbf a_{n+1}
= \mathbf F(t_{n+1}) - C \mathbf v^* - K \mathbf q^*
$$

$S$ doesn't change between steps. The code factors it once before the loop, at [solver.py:84](membrane/solver.py#L84). `solve(...)` on line 109 uses that factorization to solve the system above. This is where $M$ enters each step.

**Correct** ([solver.py:110-111](membrane/solver.py#L110-L111)):

$$
\mathbf q_{n+1} = \mathbf q^* + \frac{\Delta t^2}{4}\mathbf a_{n+1},
\qquad
\mathbf v_{n+1} = \mathbf v^* + \frac{\Delta t}{2}\mathbf a_{n+1}
$$

## 7. Quantities recorded for the plots

These don't affect the simulation. They are computed from $\mathbf q$ and $\mathbf v$ after each step.

**Energy** ([solver.py:120-121](membrane/solver.py#L120-L121)):

$$
E_\text{kin} = \tfrac12 \mathbf v^\mathsf{T} M \mathbf v, \qquad
E_\text{pot} = \tfrac12 \mathbf q^\mathsf{T} K \mathbf q
$$

**Work in and energy dissipated** ([solver.py:115-117](membrane/solver.py#L115-L117)). These are running totals. Each step adds:

$$
\Delta W = \left(\mathbf q_{n+1} - \mathbf q_n\right)^\mathsf{T} \frac{\mathbf F_n + \mathbf F_{n+1}}{2},
\qquad
\Delta D = \Delta t\, \bar{\mathbf v}^\mathsf{T} C \bar{\mathbf v},
\qquad
\bar{\mathbf v} = \frac{\mathbf v_n + \mathbf v_{n+1}}{2}
$$

With these definitions, $E_{n+1} - E_n = \Delta W - \Delta D$ holds exactly for the update rule in section 6.

**Probe displacement** ([solver.py:122](membrane/solver.py#L122)). The displacement at a point $(x_0, y_0)$ is

$$
u(x_0, y_0) = \sum_i \varphi_i(x_0, y_0)\, q_i
$$

where the $\varphi_i(x_0, y_0)$ are the barycentric coordinates of the point in the triangle that contains it.

**Animation frames** ([solver.py:124-126](membrane/solver.py#L124-L126)). $\mathbf q$ is stored at evenly spaced steps. Nothing is interpolated between steps.

## 8. Initial conditions

[loads.py:22-75](membrane/loads.py#L22-L75)

These produce $\mathbf q_0$ or $\mathbf v_0$ by evaluating a function at the interior nodes. The one formula that needs explaining is the strike. With a Gaussian shape $g_i$ centred on the strike point, it sets the velocity

$$
v_{0,i} = \frac{J\, g_i}{\sigma \sum_j A_j g_j}
$$

so the membrane's total momentum $\sigma \sum_i A_i v_{0,i}$ equals the impulse $J$.
