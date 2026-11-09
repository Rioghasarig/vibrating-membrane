# Vibrating Circular Membrane

A circular membrane of radius $a$ with a fixed boundary, uniform surface mass density $\sigma$, and uniform tension $T$ (force per unit length). The membrane occupies the disk $D = \{(x, y) : x^2 + y^2 \le a^2\}$, and the small transverse displacement is $u(x,y,t)$.

## Lagrangian

The Lagrangian is kinetic minus potential energy:

$$
L[u,\dot u] = \iint_D \left[ \frac{\sigma}{2} u_t^2 - \frac{T}{2} |\nabla u|^2 \right] dx \, dy
$$

The gradient term is:

$$
|\nabla u|^2 = u_x^2 + u_y^2
$$

So the Lagrangian density (per unit area) is:

$$
\mathcal{L} = \frac{\sigma}{2} u_t^2 - \frac{T}{2} \left( u_x^2 + u_y^2 \right)
$$

The fixed boundary gives the condition:

$$
u(x,y,t) = 0 \quad \text{for } x^2 + y^2 = a^2
$$

### Where the terms come from

- **Kinetic energy:** each area element $dx\,dy$ has mass $\sigma\,dx\,dy$ and moves with speed $u_t$.
- **Potential energy:** this is the tension times the increase in area. The exact expression is $T\left(\sqrt{1 + |\nabla u|^2} - 1\right)$, which reduces to $\frac{T}{2}|\nabla u|^2$ for small slopes ($|\nabla u| \ll 1$). This approximation is what makes the theory linear.

## Finite element discretization

Write $\mathbf{x} = (x, y)$ for a point on the mesh.

### 1. Discretize $u$

Mesh the disk with linear triangles, placing the boundary nodes on the circle $x^2 + y^2 = a^2$. Let nodes $1, \dots, N$ be the interior nodes, and let $\varphi_i$ be the piecewise-linear hat function with $\varphi_i(\mathbf{x}_j) = \delta_{ij}$. Then approximate:

$$
u(\mathbf{x}, t) \approx \sum_{i=1}^{N} q_i(t)\, \varphi_i(\mathbf{x})
$$

The boundary nodes are left out of the sum, so $u = 0$ holds on the mesh boundary.

### 2. Substitute into the Lagrangian

The integral over the mesh splits into a sum of integrals over the triangles $\Omega_e$:

$$
L \approx \sum_{e} \left[ \frac{\sigma}{2} \sum_{i,j \in e} \dot q_i \dot q_j \int_{\Omega_e} \varphi_i \varphi_j \, dA
- \frac{T}{2} \sum_{i,j \in e} q_i q_j \int_{\Omega_e} \nabla\varphi_i \cdot \nabla\varphi_j \, dA \right]
$$

Here $i, j \in e$ means $i$ and $j$ run over the interior nodes that are vertices of triangle $e$. These are the only terms that survive, because $\varphi_i$ is zero on any triangle that doesn't have node $i$ as a vertex.

Written in matrix form, with $\mathbf{q} = (q_1, \dots, q_N)^\mathsf{T}$, this is:

$$
L = \frac{1}{2} \dot{\mathbf{q}}^\mathsf{T} M \dot{\mathbf{q}} - \frac{1}{2} \mathbf{q}^\mathsf{T} K \mathbf{q}
$$

where $M$ and $K$ are the $N \times N$ mass and stiffness matrices. Each is built from the matrices of the individual triangles.

Let $\varphi^e_1, \varphi^e_2, \varphi^e_3$ be the hat functions of the three vertices of triangle $e$, and let $n^e_i$ be the global node number of vertex $i$.

Each triangle $e$ has a $3 \times 3$ mass matrix $M^e$:

$$
M^e_{ij} = \sigma \int_{\Omega_e} \varphi^e_i \varphi^e_j \, dA, \qquad i, j \in \{1, 2, 3\}
$$

Each $M^e$ is placed in an $N \times N$ matrix $\hat M^e$, whose entries at the triangle's global node numbers are:

$$
\hat M^e_{n^e_i \, n^e_j} = M^e_{ij}
$$

All other entries of $\hat M^e$ are zero, and any entry of $M^e$ that involves a boundary vertex is dropped, because boundary nodes have no row or column in $M$. The full mass matrix is the sum of these $N \times N$ matrices:

$$
M = \sum_{e} \hat M^e \quad \text{(mass matrix)}
$$

Each triangle $e$ also has a $3 \times 3$ stiffness matrix $K^e$:

$$
K^e_{ij} = T \int_{\Omega_e} \nabla\varphi^e_i \cdot \nabla\varphi^e_j \, dA, \qquad i, j \in \{1, 2, 3\}
$$

Each $K^e$ is placed in an $N \times N$ matrix $\hat K^e$, whose entries at the triangle's global node numbers are:

$$
\hat K^e_{n^e_i \, n^e_j} = K^e_{ij}
$$

All other entries of $\hat K^e$ are zero, and any entry of $K^e$ that involves a boundary vertex is dropped, because boundary nodes have no row or column in $K$. The full stiffness matrix is the sum of these $N \times N$ matrices:

$$
K = \sum_{e} \hat K^e \quad \text{(stiffness matrix)}
$$

The resulting matrices $M$ and $K$ are sparse, symmetric, and positive definite.

### 3. Apply Lagrange's equations to the coordinates $q_i$

Setting

$$
\frac{d}{dt}\frac{\partial L}{\partial \dot q_i} - \frac{\partial L}{\partial q_i} = 0
$$

for each $i$ gives $N$ coupled linear equations:

$$
M \ddot{\mathbf{q}} + K \mathbf{q} = 0
$$

### 4. Add damping

Suppose each area element feels a damping force proportional to its velocity, with coefficient $c(x, y) \ge 0$ per unit area:

$$
f_d = -c(x, y)\, u_t
$$

If only part of the membrane is damped, set $c = 0$ elsewhere.

This force can't be written as part of $L$. Instead, it enters Lagrange's equations as a generalized force $Q_i$ on each coordinate:

$$
\frac{d}{dt}\frac{\partial L}{\partial \dot q_i} - \frac{\partial L}{\partial q_i} = Q_i, \qquad i = 1, \dots, N
$$

The generalized forces are found with the principle of virtual work.

**Virtual displacement.** Give the membrane a virtual displacement $\delta u$: a small imagined change of shape at a fixed instant, consistent with the fixed boundary. Using the same hat functions as for $u$:

$$
\delta u(\mathbf{x}) = \sum_{i=1}^{N} \delta q_i\, \varphi_i(\mathbf{x})
$$

The $\delta q_i$ are arbitrary and independent. Every interior hat function is zero on the boundary, so $\delta u = 0$ there, as the constraint requires.

**Virtual work of the damping force.** The force on area element $dA$ is $f_d\, dA$. It moves through $\delta u$, so the total virtual work is:

$$
\delta W = \iint_D f_d\, \delta u \, dA = -\iint_D c\, u_t\, \delta u \, dA
$$

**Substitute the hat-function expansions.** With $u_t \approx \sum_j \dot q_j \varphi_j$:

$$
\delta W = -\iint_D c \left( \sum_{j=1}^{N} \dot q_j \varphi_j \right) \left( \sum_{i=1}^{N} \delta q_i \varphi_i \right) dA
= -\sum_{i=1}^{N} \delta q_i \sum_{j=1}^{N} \left( \iint_D c\, \varphi_i \varphi_j \, dA \right) \dot q_j
$$

The $\dot q_j$ and $\delta q_i$ don't depend on position, so they come outside the integral.

**Read off the generalized forces.** The generalized force $Q_i$ is defined as the coefficient of $\delta q_i$ in the virtual work:

$$
\delta W = \sum_{i=1}^{N} Q_i\, \delta q_i
$$

Because the $\delta q_i$ are independent, matching coefficients gives:

$$
Q_i = -\sum_{j=1}^{N} C_{ij}\, \dot q_j, \qquad C_{ij} = \iint_D c\, \varphi_i \varphi_j \, dA
$$

$C$ is the $N \times N$ damping matrix. As in step 2, the integral splits into a sum over triangles, so $C$ is built from the triangles in the same way as $M$ and $K$. Each triangle $e$ has a $3 \times 3$ damping matrix $C^e$:

$$
C^e_{ij} = \int_{\Omega_e} c\, \varphi^e_i \varphi^e_j \, dA, \qquad i, j \in \{1, 2, 3\}
$$

Each $C^e$ is placed in an $N \times N$ matrix $\hat C^e$, whose entries at the triangle's global node numbers are:

$$
\hat C^e_{n^e_i \, n^e_j} = C^e_{ij}
$$

All other entries of $\hat C^e$ are zero, and any entry of $C^e$ that involves a boundary vertex is dropped. The full damping matrix is the sum of these $N \times N$ matrices:

$$
C = \sum_{e} \hat C^e \quad \text{(damping matrix)}
$$

$C$ is sparse and symmetric, with the same pattern of nonzero entries as $M$. It is not diagonal: moving one node moves the membrane over every triangle around that node, so the damping force acts on the neighbouring nodes too. If $c$ is uniform, $C = (c/\sigma)\, M$.

The left-hand side of Lagrange's equations is the same as in step 3:

$$
\frac{d}{dt}\frac{\partial L}{\partial \dot q_i} = \sum_{j=1}^{N} M_{ij} \ddot q_j,
\qquad
-\frac{\partial L}{\partial q_i} = \sum_{j=1}^{N} K_{ij} q_j
$$

Setting this equal to $Q_i$ and moving $Q_i$ to the left gives:

$$
\sum_{j=1}^{N} M_{ij} \ddot q_j + \sum_{j=1}^{N} C_{ij} \dot q_j + \sum_{j=1}^{N} K_{ij} q_j = 0, \qquad i = 1, \dots, N
$$

In matrix form, this is:

$$
M \ddot{\mathbf{q}} + C \dot{\mathbf{q}} + K \mathbf{q} = 0
$$
