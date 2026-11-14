# Godunov Marching Method for Two-Dimensional Supersonic Flows of an Ideal Gas


Many problems of supersonic aerodynamics are described by the steady Euler equations in a domain where the flow is supersonic everywhere. Examples are flow past pointed bodies, flow past blunt bodies downstream of the sonic line, flow in the supersonic part of a nozzle, underexpanded and overexpanded jets, and air intakes. Such flows have an important property: disturbances propagate only downstream, inside Mach wedges. If the direction $x$ is chosen so that the velocity component along it is supersonic everywhere, the solution in the section $x = \mathrm{const}$ is completely determined by what happens upstream. The steady problem can then be solved by "marching": stepping from one section to the next, just as unsteady equations are integrated in time.

The Godunov marching method carries over to this situation the idea of S. K. Godunov's classical scheme [1] for unsteady flows. The solution on each section is represented as piecewise constant, and the fluxes across cell boundaries are obtained from the exact solution of the problem of the breakup of an arbitrary discontinuity. The only difference is that the role of time is now played by the marching coordinate $x$, and the role of the one-dimensional discontinuity breakup is played by the steady interaction of two uniform supersonic streams. The method was developed in the early 1970s by M. Ya. Ivanov, A. N. Kraiko and co-authors [3] and is described in detail in the monograph [2].

The method is conservative and monotone, captures shocks and tangential discontinuities, and at the same time requires orders of magnitude less computation than time-marching to a steady state: the two-dimensional steady problem is reduced to a sequence of one-dimensional ones.

## 1. Equations and Applicability Condition

We consider steady planar flows of an inviscid, non-heat-conducting perfect gas. Here $x$ is the marching coordinate (usually along the body or channel), $y$ is the transverse coordinate, $u, v$ are the velocity components along $x$ and $y$, and $p$ and $\rho$ are the pressure and density.

In divergence form the conservation laws read

$$
\frac{\partial \mathbf F}{\partial x} + \frac{\partial \mathbf G}{\partial y} = 0,
$$

$$
\mathbf F = \begin{pmatrix} \rho u \\ \rho u^2  + p \\ \rho u v \\ \rho u H \end{pmatrix},\qquad
\mathbf G = \begin{pmatrix} \rho v \\ \rho u v \\ \rho v^2 + p \\ \rho v H \end{pmatrix},
$$

where the total enthalpy is ($\gamma$ is the ratio of specific heats)

$$
H = \frac{\gamma}{\gamma - 1}\,\frac{p}{\rho} + \frac{q^2}{2}.
$$

We will also need
$$c = \sqrt{\gamma p/\rho} \quad \text{— speed of sound}$$
$$q = \sqrt{u^2 + v^2} \quad \text{— flow speed}$$
$$M = q/c \quad \text{— Mach number}$$
$$\theta = \arctan(v/u) \quad \text{— flow angle}$$
$$\mu = \arcsin(1/M) \quad \text{— Mach angle}$$

### 1.1. Characteristics and the Marching Condition

**Characteristics.** In a supersonic flow a small disturbance, for example from a small bump on a wall, does not spread in all directions but only downstream and only along certain lines. These lines are called characteristics. Three characteristics pass through every point of the flow:

- the **streamline** with slope $dy/dx = \tan\theta$. Gas particles move along it. A particle does not exchange heat with its neighbours and experiences no friction, so along a streamline (away from shocks) its entropy and total enthalpy $H$ are conserved;
- the **two Mach lines** $C^+$ and $C^-$. They are inclined upward and downward from the streamline by the Mach angle $\mu$, and weak acoustic disturbances of pressure and flow direction travel along them. Their slopes are

$$
\lambda^\pm = \tan(\theta \pm \mu) = \frac{uv \pm c\sqrt{q^2 - c^2}}{u^2 - c^2}.
$$

```
  y
  ^                                         C+ :  dy/dx = tan(θ + μ)
  |                                      /
  +--> x                               /
                                     /
       UPSTREAM                    /  μ             DOWNSTREAM
                                 /
   flow   ------------------->  *--------------------  streamline :  dy/dx = tan θ
                                 \
                                   \  μ
                                     \
                                       \
                                         C- :  dy/dx = tan(θ − μ)
```

Mach lines exist only for $M > 1$: only then is the angle $\mu = \arcsin(1/M)$ defined. A disturbance from a point acts only inside the wedge between $C^+$ and $C^-$ downstream (the Mach wedge), while the flow upstream "does not know" about it. A system of equations with this property is called hyperbolic. For $M < 1$ there are no Mach lines, a disturbance from any point is felt throughout the domain, including upstream, and marching is impossible.

**Marching condition.** The marching method computes the flow step by step along the $x$ axis. First the gas parameters are prescribed on the vertical line $x = x_0$ (called a section). From them the parameters on the neighbouring line $x_1 = x_0 + \Delta x$ are found, then from $x_1$ those on $x_2$, and so on, as in time integration: each step uses only what has already been computed. The condition $M > 1$ is not enough for this. All characteristics must go **forward** in $x$. If some Mach line turned backward, the current line would be influenced by the flow to the right of it, which has not been computed yet.

Hence no Mach line may reach the vertical: the angles $\theta \pm \mu$ must lie in the interval $(-\pi/2,\ \pi/2)$. This is equivalent to $|\theta| + \mu < \pi/2$, or $\cos\theta > \sin\mu = 1/M$. Multiplying both sides by $q$ and using $q\cos\theta = u$ and $q/M = c$, we obtain the **marching condition**:

$$
u > c.
$$

That is, the velocity component along the $x$ axis must be supersonic. If the flow is strongly inclined to the $x$ axis, this condition may be violated even for $M > 1$. In that case a curvilinear marching coordinate is chosen, roughly along the flow.

When the condition holds, no characteristic is tangent to the line $x = \mathrm{const}$. Therefore the solution for $x > x_0$ is determined only by the data on the section $x = x_0$ and the boundary conditions, and the coordinate $x$ plays the role of time.

**What happens at the limit of applicability.** As $u \to c$, the denominator $u^2 - c^2$ in the formula for $\lambda^\pm$ tends to zero. One of the Mach lines becomes vertical (parallel to the $y$ axis), and the marching stops. In a computation this shows up as a sharp decrease of the step size.

## 2. Computational Grid

The computational domain is bounded below by the line $y = y_w(x)$ (body surface, channel wall or symmetry plane) and above by the line $y = y_e(x)$ (bow shock, jet boundary, wall or simply a line in the undisturbed flow). The marching coordinate is divided into layers with a variable step $\Delta x_n = x_{n+1} - x_n$, chosen during the computation.

On each layer the segment $[y_w, y_e]$ is divided into $N$ cells. Cells are numbered with integer indices $j = 1, \dots, N$ from bottom to top, and the boundaries between them (nodes) with half-integer indices: cell $j$ lies between nodes $y_{j-1/2}$ and $y_{j+1/2}$. The lowest node $y_{1/2}$ lies on the lower boundary of the domain, the highest $y_{N+1/2}$ on the upper one.

The position of a node is given by the fraction $\eta_{j+1/2}$ of the domain height at which it lies:

$$
y^n_{j+1/2} = y_w(x_n) + \eta_{j+1/2}\,\bigl[y_e(x_n) - y_w(x_n)\bigr],\qquad j = 0, 1, \dots, N.
$$

The fractions increase with the node number from $\eta_{1/2} = 0$ (lower boundary) to $\eta_{N+1/2} = 1$ (upper boundary). The set of fractions is chosen once and is the same on all layers: as $y_w(x)$ and $y_e(x)$ change, the grid stretches with the domain, but in every section it divides it in the same proportions. For a uniform grid $\eta_{j+1/2} = j/N$. If clustering toward a wall or a shock is needed, non-uniform fractions are used, e.g. for $N = 4$: $0;\ 0.1;\ 0.25;\ 0.5;\ 1$ (clustering toward the wall).

Between layers $x_n$ and $x_{n+1}$, nodes with the same index are connected by straight segments. This produces quadrilateral cells: two sides lie on the layers, the other two are the inclined faces $j \pm 1/2$ with slopes

$$
s_{j+1/2} = \frac{y^{n+1}_{j+1/2} - y^n_{j+1/2}}{\Delta x_n}.
$$

```
              x_n                      x_{n+1}
               |                          |
 y^n_{j+1/2}   *--------------------------*   y^{n+1}_{j+1/2}    face j+1/2, slope s_{j+1/2}
               |                          |
     F^n_j     |         cell j           |   F^{n+1}_j
               |                          |
 y^n_{j-1/2}   *--------------------------*   y^{n+1}_{j-1/2}    face j-1/2, slope s_{j-1/2}
               |                          |
   (faces are straight segments, in general inclined and not parallel to each other)
```

On layer $x_n$ the gas parameters $W^n_j = (p, \rho, u, v)^n_j$ are assumed constant in each cell. This is the analogue of the piecewise-constant spatial distribution in Godunov's scheme. More precisely, $\mathbf F^n_j$ is the average of the mass, momentum and energy fluxes across the segment $[y^n_{j-1/2},\ y^n_{j+1/2}]$ of length

$$
h^n_j = y^n_{j+1/2} - y^n_{j-1/2}.
$$

For simplicity the number of cells $N$ is the same on all layers here. The scheme does not require this: it relies only on the conservation law for each cell, so cells can be split or merged during marching, for example when the domain strongly expands or contracts. Such regridding is not considered below.

## 3. Difference Scheme

### 3.1. Integral Conservation Law for a Cell

Integrate the equations over the cell $\Omega_j$ enclosed between layers $x_n$ and $x_{n+1}$. By Green's theorem (counter-clockwise traversal)

$$
\oint_{\partial\Omega_j} \bigl(\mathbf F\,dy - \mathbf G\,dx\bigr) = 0.
$$

The side on layer $x_n$ contributes $-h^n_j \mathbf F^n_j$, the side on layer $x_{n+1}$ contributes $h^{n+1}_j \mathbf F^{n+1}_j$. On an inclined face $dy = s\,dx$, and its contribution is $\pm\int(\mathbf G - s\mathbf F)\,dx$. The quantity

$$
\boldsymbol\Phi = \mathbf G - s\,\mathbf F
$$

is the flux through a moving (inclined) face per unit length in $x$. If $\boldsymbol\Phi$ is constant along the face, its contribution is exactly $\pm\Delta x_n\,\boldsymbol\Phi$.

### 3.2. Scheme

We obtain the scheme

$$
\boxed{\;
h^{n+1}_j \mathbf F^{n+1}_j = h^n_j \mathbf F^n_j - \Delta x_n\bigl(\boldsymbol\Phi_{j+1/2} - \boldsymbol\Phi_{j-1/2}\bigr)
\;}
$$

The flux through a face is computed from the state $W^*$ on that face:

$$
\boldsymbol\Phi = \begin{pmatrix} m \\ m u - s p \\ m v + p \\ m H \end{pmatrix},\qquad m = \rho\,(v - s u).
$$

Here $m$ is the mass flow rate through the face, and the terms with $p$ are the pressure forces acting on the face with normal $(-s, 1)$.

Godunov's key idea lies in the choice of $W^*$. At the point $(x_n, y^n_{j+1/2})$ two uniform streams meet: $W_j$ from below and $W_{j+1}$ from above. Their interaction is a steady discontinuity breakup problem (Section 5). Its solution is self-similar, i.e. it depends only on the slope of the ray $\xi = (y - y^n_{j+1/2})/(x - x_n)$, and face $j+1/2$ is exactly a ray from this point with slope $\xi = s_{j+1/2}$. Therefore the state is constant along the whole face and equals the breakup solution on the ray $\xi = s_{j+1/2}$. The flux through the face is computed exactly as long as waves from the neighbouring faces have not reached it; this condition limits the step.

### 3.3. Geometric Conservation Law

The grid in the marching method is moving: from layer to layer the nodes shift and the cell length changes ($h^{n+1}_j \ne h^n_j$). In the "$x$ is time" analogy this is a computation on a moving grid. The scheme must be consistent with this motion: the grid deformation itself must not generate fluxes of mass, momentum and energy. The minimal requirement is that a uniform flow remains uniform under any motion of the nodes. This is called the **geometric conservation law** (GCL).

**Geometric identity.** Substituting the definitions of $h$ and $s$ from Section 2, we get

$$
\Delta x_n\,(s_{j+1/2} - s_{j-1/2}) = \bigl(y^{n+1}_{j+1/2} - y^n_{j+1/2}\bigr) - \bigl(y^{n+1}_{j-1/2} - y^n_{j-1/2}\bigr) = h^{n+1}_j - h^n_j. \tag{3.3}
$$

The change in cell length equals the difference between the displacements of its upper and lower faces over the step. (Green's theorem for $\oint dy = 0$ around the cell contour gives the same result.) The identity is exact for any node placement, because the slopes $s$ are taken as the slopes of chords between the same nodes from which the lengths $h$ are computed. Therefore in the code the face slope must be computed precisely from the chord, not as the derivative $y'(x)$.

**Uniform flow.** Suppose all cells on layer $x_n$ have the same state $W$. For $W_L = W_R$ the breakup problem produces no waves, and on any face, for any slope, the state equals $W$. Hence $\boldsymbol\Phi_{j\pm1/2} = \mathbf G - s_{j\pm1/2}\mathbf F$, and $\mathbf G$ cancels in the flux difference:

$$
\boldsymbol\Phi_{j+1/2} - \boldsymbol\Phi_{j-1/2} = -\bigl(s_{j+1/2} - s_{j-1/2}\bigr)\,\mathbf F.
$$

Substituting this into the scheme of Section 3.2 and using (3.3), we obtain

$$
h^{n+1}_j \mathbf F^{n+1}_j = h^n_j \mathbf F + \bigl(h^{n+1}_j - h^n_j\bigr)\mathbf F = h^{n+1}_j\,\mathbf F,
$$

that is, $\mathbf F^{n+1}_j = \mathbf F$, and by Section 4 also $W^{n+1}_j = W$.

This shows the role of the term $-s\,\mathbf F$ in the flux $\boldsymbol\Phi$: it accounts for the motion of the face, which "sweeps" a strip of gas. If $\mathbf G$ were used instead of $\boldsymbol\Phi$, as for a horizontal face, we would get $\mathbf F^{n+1}_j = \mathbf F\, h^n_j / h^{n+1}_j$: in an expanding domain the flux would weaken for no physical reason.

Preservation of a uniform flow is the first thing worth checking in a program.

## 4. From Fluxes to Gas Parameters

After a step the vector $\mathbf F^{n+1}_j = (f_1, f_2, f_3, f_4)$ is known, and $p, \rho, u, v$ must be recovered from it. Immediately we find

$$
v = \frac{f_3}{f_1},\qquad H = \frac{f_4}{f_1},\qquad P \equiv \frac{f_2}{f_1} = u + \frac{p}{\rho u}.
$$

Substituting $p = f_1(P - u)$ and $\rho = f_1/u$ into the expression for $H$, we obtain a quadratic equation for $u$:

$$
(\gamma + 1)\,u^2 - 2\gamma P\,u + 2(\gamma - 1)\Bigl(H - \frac{v^2}{2}\Bigr) = 0.
$$

Its roots are

$$
u_{1,2} = \frac{\gamma P \pm \sqrt{D}}{\gamma + 1},\qquad D = \gamma^2 P^2 - 2(\gamma^2 - 1)\Bigl(H - \frac{v^2}{2}\Bigr),
$$

and their product $u_1 u_2 = a_*^2 = \frac{2(\gamma-1)}{\gamma+1}\bigl(H - \frac{v^2}{2}\bigr)$ is Prandtl's relation for a normal shock. The two roots are two states with identical fluxes $\mathbf F$, i.e. the states ahead of and behind a shock perpendicular to the $x$ axis. The marching (supersonic in $x$) solution corresponds to the larger root:

$$
u = \frac{\gamma P + \sqrt{D}}{\gamma + 1},\qquad p = f_1\,(P - u),\qquad \rho = \frac{f_1}{u}.
$$

The larger root satisfies $u \ge a_*$, which is equivalent to $u \ge c$. Therefore the scheme cannot produce a cell with $u < c$: instead $D < 0$ arises, meaning that no state with such fluxes exists. This is the analogue of channel choking, when one tries to pass a larger mass flow through a section than the sonic regime allows. A negative or near-zero discriminant means that the computation has gone beyond the applicability limits of the marching method (Section 9) or that the step is too large.

## 5. Steady Discontinuity Breakup Problem

The scheme of Section 3.2 requires the flux $\boldsymbol\Phi$ on every face, and for this the gas state $W^*$ on it is needed. But on layer $x_n$ the solution is piecewise constant: at face $j+1/2$ the parameters jump from $W_j$ to $W_{j+1}$. Which of the two states to substitute into $\boldsymbol\Phi$ (or some third one) is not clear in advance.

Simple recipes do not work. The average $(W_j + W_{j+1})/2$ gives a central scheme: it is unstable or produces oscillations at shocks. Taking the "upwind" state is not possible either: in a supersonic flow disturbances travel across the flow in both directions, upward along $C^+$ and downward along $C^-$, so there is no upwind side.

Godunov proposed not to guess but to solve exactly an auxiliary physical problem: what happens downstream if two different uniform supersonic streams meet at one point? This is the **steady discontinuity breakup problem**. The discontinuity itself generally cannot persist: the pressures and velocity directions on its two sides differ. Downstream it breaks up into two waves (shocks or rarefaction fans) and a slip line between them. The state that this solution gives on the face is what is substituted into $\boldsymbol\Phi$.

Everything required of the scheme follows automatically from the exact solution:

- information reaches the face only from the side from which the characteristics actually carry it;
- the scheme is monotone and produces no oscillations at shocks;
- captured shocks satisfy the Rankine–Hugoniot relations.

Besides the fluxes, the solution gives the wave slopes used to choose the step $\Delta x_n$ (Section 7). The same problem, with minor changes, is solved at the domain boundaries: at a wall, at a free jet boundary and at a fitted bow shock (Section 6). Hence the solver of this problem is the central part of the program.

The problem is fully analogous to the one-dimensional unsteady Riemann problem of the classical Godunov scheme. The correspondence is:

| | One-dimensional unsteady problem | Steady supersonic problem |
|---|---|---|
| Role of time | $t$ | $x$ |
| Role of space | $x$ | $y$ |
| Initial data | two gas states left and right of a point | two supersonic streams below and above a point |
| Outer waves | shock wave or rarefaction wave | oblique shock or Prandtl–Meyer fan |
| Middle discontinuity | contact discontinuity | slip line |
| Continuous across it | pressure $p^*$ and velocity $u^*$ | pressure $p^*$ and flow angle $\theta^*$ |
| Self-similar variable | $x/t$ | $\xi = y/x$ |

The chapter is organized as follows. Section 5.1 describes the structure of the solution and shows that the problem reduces to a single equation for the pressure $p^*$. Sections 5.2 and 5.3 give formulas for the two types of waves. In Section 5.4 the equation for $p^*$ is written out and solved. Sections 5.5 and 5.6 show how to find the state on a face from $p^*$. Section 5.7 discusses cases where no solution exists.

### 5.1. Statement and Structure of the Solution

**Statement.** On the line $x = x_n$ at the point $y = y_{j+1/2}$ two uniform streams meet: the lower one $W_L = W_j$ (for $y < y_{j+1/2}$) and the upper one $W_R = W_{j+1}$ (for $y > y_{j+1/2}$). In both $u > c$. The flow for $x > x_n$ is to be found.

It is convenient to describe each state not by the velocity components but by its magnitude and angle: $(p, \rho, q, \theta)$. The angle $\theta$ is measured from the $x$ axis counter-clockwise, so "the flow turned upward" means that $\theta$ increased.

**Self-similarity.** The problem contains no characteristic length: there is only a point and two uniform streams. If the flow pattern is stretched by any factor about the discontinuity point, it remains a solution of the same problem. Hence the solution depends not on $x$ and $y$ separately but only on the slope of the ray drawn from the discontinuity point:

$$
\xi = \frac{y - y_{j+1/2}}{x - x_n}.
$$

All parameters are constant along each such ray, and the boundaries between solution regions are also rays.

**Three elements of the solution.** Disturbances propagate only along characteristics, of which there are three families (Section 1.1). Therefore three elements emanate from the discontinuity point:

- a wave of the $C^-$ family going downward into the lower stream;
- a slip line (tangential discontinuity) — the streamline leaving the discontinuity point. It separates gas that came from the lower stream from gas that came from the upper one;
- a wave of the $C^+$ family going upward into the upper stream.

Each of the two waves is either an oblique shock or a Prandtl–Meyer rarefaction fan.

```
  x = x_n
     |                                     xi_R : shock or fan C+
     |    W_R                        _.-´
     |                         _.-´
     |                   _.-´            W_R*
     |             _.-´
     |       _.-´
     *-----------------------------------------  xi_0 : slip line
     |       `-._
     |             `-._              W_L*
     |                   `-._
     |    W_L                  `-._
     |                               `-._  xi_L : shock or fan C-
```

The waves divide the domain into four zones with constant parameters: the undisturbed streams $W_L$ and $W_R$ outside, and two new states $W_{L*}$ and $W_{R*}$ between the waves and the slip line. It is $W_{L*}$ and $W_{R*}$ that must be found.

**Conditions on the slip line.** Gas does not cross the slip line, so on both sides of it the velocity is directed along it: the flow angles are equal. The pressures are also equal, otherwise the line would not be in equilibrium. Thus

$$
p_{L*} = p_{R*} = p^*,\qquad \theta_{L*} = \theta_{R*} = \theta^*.
$$

The other quantities (density, flow speed, entropy, total enthalpy) generally differ on the two sides of the slip line: these are two different gases, each with its own history. The total enthalpy $H_k$ ($k = L, R$) does not change in either a shock or a fan, so $H_{L*} = H_L$ and $H_{R*} = H_R$.

**An intuitive picture: a fluid wall.** For each stream the slip line is the same as a solid wall inclined at angle $\theta^*$. The upper stream flows over this wall. If the wall is turned into it ($\theta^* > \theta_R$), we get flow over a wedge with a compression shock. If the wall turns away from it ($\theta^* < \theta_R$), we get flow around a convex corner with a rarefaction fan. For the lower stream everything is mirrored. Both problems are familiar from a gas dynamics course. The only new thing here is that the slope of the "wall" is not known in advance. It must be chosen so that the pressures on its two sides coincide.

**Each wave is specified by a single number.** Let the incoming stream $W_k$ and the pressure $p$ behind the wave be known. Then everything is known:

- the wave type: for $p > p_k$ it is a shock, for $p < p_k$ a rarefaction fan. There are no other options: a rarefaction shock is forbidden by the second law of thermodynamics (entropy would decrease in it), and a compression wave emanating from a single point is impossible because its characteristics converge rather than diverge;
- all parameters behind the wave, including the angle by which it turns the flow (formulas of Sections 5.2 and 5.3).

Denote by $f_k(p)$ the angle by which the wave in stream $k$ turns the flow while changing its pressure from $p_k$ to $p$. The sign is chosen as follows: in compression ($p > p_k$) $f_k > 0$, in rarefaction ($p < p_k$) $f_k < 0$, and for $p = p_k$ there is no wave and $f_k = 0$.

The direction of turning depends on the wave family:

| | Compression (shock), $p^* > p_k$ | Rarefaction (fan), $p^* < p_k$ |
|---|---|---|
| Upper stream, wave $C^+$ | turns upward, $\theta$ increases | turns downward, $\theta$ decreases |
| Lower stream, wave $C^-$ | turns downward, $\theta$ decreases | turns upward, $\theta$ increases |

This is easy to remember: a shock turns the flow toward the side to which the shock itself is inclined, and a fan turns it to the opposite side. A check example is a wedge: the flow over a wall turned upward is compressed by a $C^+$ shock and turns upward. In formulas:

$$
\theta_{R*} = \theta_R + f_R(p^*),\qquad \theta_{L*} = \theta_L - f_L(p^*).
$$

**What the problem reduces to.** Equating these two angles gives one equation in one unknown $p^*$ (Section 5.4). All other parameters are expressed explicitly in terms of $p^*$ (Section 5.5). The one-dimensional Riemann problem is solved the same way, except that there $p^*$ determines not the angle but the velocity behind each wave.

The function $f_k$ is the same for both streams: it depends only on $p$, $p_k$ and $M_k$. The difference between the upper and lower streams is only in the sign with which it enters the formula for the angle.

### 5.2. Oblique Shock

This branch of the function $f_k$ applies for $p \ge p_k$.

```
                                           /  shock (wave C+)
                                         /
                                       /        W_k*
                                     /
      W_k                          /
                                 /                     _.-´  slip line
   ----------------------->    / )               _.-´
                             /    β        _.-´      )
                           /       ) _.-´            δ
                         /     _.-´ )                )
                       *----------------------------------------  direction of W_k
                       O

   arc β: from the direction of W_k to the shock
   arc δ: from the direction of W_k to the slip line
```

The figure shows a $C^+$ shock in the upper stream. The incoming stream $W_k$ passes through the shock, is compressed from $p_k$ to $p$ and turns upward by the angle $\delta_k$. The shock itself is inclined to the direction of $W_k$ at a larger angle $\beta_k$, so the streamlines behind it are less steep than the shock. For a $C^-$ shock in the lower stream the picture is mirrored: both the shock and the flow behind it deflect downward.

**Where the formulas come from.** An oblique shock is a normal shock with a motion along the front superimposed on it. Decompose the velocity of the incoming stream into the component normal to the front $q_k\sin\beta_k$ and the tangential component $q_k\cos\beta_k$. The tangential component does not change across the shock, while the normal component undergoes an ordinary normal shock with Mach number $M_k\sin\beta_k$. Therefore the pressure ratio $r = p/p_k$ is given by the normal shock formula with $M$ replaced by $M_k\sin\beta_k$:

$$
r = 1 + \frac{2\gamma}{\gamma + 1}\bigl(M_k^2\sin^2\beta_k - 1\bigr).
$$

In our problem the pressure, not the shock angle, is given. Therefore this formula must be inverted:

$$
\sin\beta_k = \frac{1}{M_k}\sqrt{1 + \frac{\gamma + 1}{2\gamma}(r - 1)}.
$$

**Flow deflection angle.** It follows from the fact that the tangential velocity component is preserved while the normal one decreases. The derivation takes four steps.

*Step 1. Velocities before and after the shock.* Denote the velocity components relative to the front: normal $q_n$ and tangential $q_\tau$. Ahead of the shock

$$
q_n = q_k\sin\beta_k,\qquad q_\tau = q_k\cos\beta_k.
$$

Behind the shock the tangential component is the same, $q_\tau' = q_\tau$. The normal one is found from mass conservation: the mass flux per unit area of the front $\rho_k q_n = \rho\, q_n'$ is the same on both sides, so

$$
q_n' = \varepsilon\, q_n,\qquad \varepsilon = \frac{\rho_k}{\rho} < 1.
$$

The gas has been compressed, and the normal velocity has decreased by the same factor by which the density increased.

*Step 2. Angle between the flow and the front.* Ahead of the shock the velocity vector makes an angle $\beta_k$ with the front: $\tan\beta_k = q_n/q_\tau$. Behind the shock the flow has turned by $\delta_k$ toward the front, and the angle between them has become $\beta_k - \delta_k$:

$$
\tan(\beta_k - \delta_k) = \frac{q_n'}{q_\tau} = \varepsilon\tan\beta_k. \tag{1}
$$

Since $\varepsilon < 1$, the velocity vector has been pressed toward the front — this is the flow deflection. By the formula for the tangent of the difference $\beta_k - \delta_k$:
$$
\tan(\beta_k - \delta_k) = \frac{\tan\beta_k - \tan\delta_k}{1 + \tan\beta_k \tan\delta_k}. \tag{2}
$$
Substituting (2) into (1) and solving for $\delta_k$:
$$
\tan\delta_k = \frac{\tan\beta_k - \varepsilon\tan\beta_k}{1 + \varepsilon\tan^2\beta_k}
= \frac{(1 - \varepsilon)\sin\beta_k\cos\beta_k}{1 - (1 - \varepsilon)\sin^2\beta_k}.
$$

The second equality is obtained by multiplying the numerator and denominator by $\cos^2\beta_k$ and replacing $\cos^2\beta_k = 1 - \sin^2\beta_k$.

*Step 3. Eliminating the density.* Write the momentum conservation normal to the front (using mass conservation $\rho q_n = \rho' q_n'$): the pressure rise equals the loss of momentum flux,

$$
p - p_k = \rho_k q_n\,(q_n - q_n') = \rho_k q_n^2\,(1 - \varepsilon).
$$

Substitute $q_n = q_k\sin\beta_k$ and note that
$$M_k^2 = \frac{q_k^2}{c_k^2} = \frac{q_k^2 \rho_k}{\gamma p_k} \Rightarrow \rho_k q_k^2 = \gamma p_k M_k^2$$

$$
(1 - \varepsilon)\sin^2\beta_k = \frac{r - 1}{\gamma M_k^2}.
$$

This is exactly the combination appearing in the formula of Step 2. After substitution

$$
\tan\delta_k = \frac{r - 1}{\gamma M_k^2 - (r - 1)}\,\cot\beta_k.
$$

*Step 4. Eliminating the shock angle.* From the formula for $\sin\beta_k$ obtained above,

$$
M_k^2\sin^2\beta_k = \frac{(\gamma + 1)\,r + (\gamma - 1)}{2\gamma},
$$

hence

$$
\cot^2\beta_k = \frac{M_k^2 - M_k^2\sin^2\beta_k}{M_k^2\sin^2\beta_k}
= \frac{2\gamma M_k^2 - (\gamma - 1) - (\gamma + 1)\,r}{(\gamma + 1)\,r + (\gamma - 1)}.
$$

Substituting the square root of this expression into the result of Step 3, we obtain the deflection angle as a function of $r$ and $M_k$ only:

$$
\tan\delta_k(p) = \frac{r - 1}{\gamma M_k^2 - (r - 1)}\,
\sqrt{\frac{2\gamma M_k^2 - (\gamma - 1) - (\gamma + 1)\,r}{(\gamma + 1)\,r + (\gamma - 1)}}.
$$

The density behind the shock is given by the Hugoniot adiabat:

$$
\frac{\rho}{\rho_k} = \frac{(\gamma + 1)\,r + (\gamma - 1)}{(\gamma - 1)\,r + (\gamma + 1)}.
$$

**What is needed for what.** Only the deflection angle enters the equation for $p^*$: on the shock branch $f_k(p) = \delta_k(p)$. The angle $\beta_k$ does not take part in the equation but determines where the shock lies: taking into account the inclination of the incoming stream, the shock ray has slope $\tan(\theta_k + \beta_k)$ for a $C^+$ wave and $\tan(\theta_k - \beta_k)$ for a $C^-$ wave (Section 5.5).

**Limiting case check.** For $r = 1$ the formulas give $\delta_k = 0$ and $\sin\beta_k = 1/M_k$, i.e. $\beta_k = \mu_k$: an infinitely weak shock does not turn the flow and coincides with a Mach line.

**Weak and strong branches.** The pressure behind the shock can vary from $p_k$ ($r = 1$) to the pressure behind a normal shock, $r_n = 1 + \frac{2\gamma}{\gamma+1}(M_k^2 - 1)$, at which $\beta_k = \pi/2$. The dependence of $\delta_k$ on $r$ on this interval is non-monotone: the angle grows from zero to the maximum $\delta_{\max}(M_k)$ and then decreases to zero at the normal shock. The increasing part is called the weak branch, the decreasing part the strong branch. Only the weak branch is used. There are two reasons: it continuously turns into a Mach line as $r \to 1$, and behind a weak shock the flow usually remains supersonic, whereas behind a strong one it is subsonic and marching is impossible.

The right end of the weak branch is the point of maximum $\delta_k$. It is inconvenient to express it in terms of $r$, but in terms of the shock angle it is found explicitly:

$$
\sin^2\beta_{\max} = \frac{(\gamma + 1)M_k^2 - 4 + \sqrt{(\gamma + 1)\bigl[(\gamma + 1)M_k^4 + 8(\gamma - 1)M_k^2 + 16\bigr]}}{4\gamma M_k^2}.
$$

Hence the largest admissible pressure ratio is

$$
r_{\max} = 1 + \frac{2\gamma}{\gamma + 1}\bigl(M_k^2\sin^2\beta_{\max} - 1\bigr).
$$

The angle $\beta_{\max}$ itself is needed only as an intermediate quantity for $r_{\max}$. Thus the shock branch is defined for $1 \le r \le r_{\max}(M_k)$.

### 5.3. Prandtl–Meyer Rarefaction Wave

This branch of the function $f_k$ applies for $p < p_k$.

```
                                         /  fan head (characteristic C+)
                                      /         /              /  fan tail
                                   /       /           /
      W_k                       /   ) /        /
                             /   /  μ  /          W_k*
   ----------------->     / /  /     )
                       *---------------------------------------------  direction of W_k
                       O   `-._              )
                                 `-._        Δθ
                                       `-._ )
                                             `-._  slip line

   arc μ:  from the direction of W_k to the fan head, μ_k = arcsin(1/M_k)
   arc Δθ: from the direction of W_k to the slip line, Δθ = ω(M) − ω(M_k) = −f_k(p)
```

The figure shows a $C^+$ rarefaction wave in the upper stream. This is the same flow as around a convex corner. A fan of $C^+$ characteristics emanates from point O:

- the fan head is the characteristic of the incoming stream $W_k$. It is inclined to its direction at the Mach angle $\mu_k$;
- inside the fan the flow accelerates continuously, the pressure drops from $p_k$ to $p$, and the velocity turns downward, away from the wave, by the angle $\Delta\theta$;
- the fan tail is the characteristic of the stream $W_{k*}$. It is inclined to the slip line at the angle $\mu_{k*} < \mu_k$ (the flow has accelerated, the Mach angle has decreased).

Unlike the shock, there is no discontinuity here: the parameters change continuously from ray to ray and are constant on each ray of the fan. For a $C^-$ wave in the lower stream the picture is mirrored: the fan spreads downward and the flow turns upward.

**Parameters behind the fan.** Since there are no discontinuities, the entropy of each particle is conserved and the flow is isentropic: $p/\rho^\gamma = \mathrm{const}$. Hence the density $\rho = \rho_k\,(p/p_k)^{1/\gamma}$. The stagnation pressure is also conserved, so the Mach number after expansion to pressure $p$ follows from the usual isentropic formula:

$$
M_k^2(p) = \frac{2}{\gamma - 1}\left[\Bigl(1 + \frac{\gamma - 1}{2}M_k^2\Bigr)\Bigl(\frac{p_k}{p}\Bigr)^{\frac{\gamma - 1}{\gamma}} - 1\right].
$$

**Turning angle.** It is expressed through the Prandtl–Meyer function

$$
\omega(M) = \sqrt{\frac{\gamma + 1}{\gamma - 1}}\,\arctan\sqrt{\frac{\gamma - 1}{\gamma + 1}(M^2 - 1)} - \arctan\sqrt{M^2 - 1}.
$$

Its meaning: $\omega(M)$ is the angle by which a sonic flow ($M = 1$) must be turned in a rarefaction fan to accelerate it to Mach number $M$. The function increases with $M$. To accelerate the flow from $M_k$ to $M$, it must be turned by $\omega(M) - \omega(M_k)$.

The same in terms of invariants: in a $C^+$ fan the quantity $\theta + \omega(M)$ is conserved at all points, in a $C^-$ fan the quantity $\theta - \omega(M)$. These invariants will be needed in Section 5.6.

According to the adopted sign convention, the turning in rarefaction is negative:

$$
f_k(p) = \omega(M_k) - \omega\bigl(M_k(p)\bigr) < 0,\qquad p < p_k.
$$

### 5.4. Equation for the Pressure

**Turning function.** Combine both branches:

$$
f_k(p) = \begin{cases}
\delta_k(p), & p \ge p_k \quad\text{(shock)},\\[4pt]
\omega(M_k) - \omega\bigl(M_k(p)\bigr), & p < p_k \quad\text{(rarefaction wave)}.
\end{cases}
$$

The function $f_k$ is continuous, monotonically increasing (on the shock branch up to $p = p_k\,r_{\max}$), and vanishes at $p = p_k$. The branches join smoothly: at $p = p_k$ they share the derivative

$$
f_k'(p_k) = Z_k \equiv \frac{\sqrt{M_k^2 - 1}}{\rho_k q_k^2} = \frac{\sqrt{M_k^2 - 1}}{\gamma p_k M_k^2}.
$$

The coefficient $Z_k$ is known from the linear theory of supersonic flows: a weak wave of any type that changes the pressure by $\Delta p$ turns the flow by the angle $Z_k\,\Delta p$.

**Equation.** The condition $\theta_{L*} = \theta_{R*}$, i.e. $\theta_L - f_L(p^*) = \theta_R + f_R(p^*)$, gives

$$
\Theta(p) \equiv f_L(p) + f_R(p) - (\theta_L - \theta_R) = 0.
$$

This is the exact analogue of the equation $f_L(p) + f_R(p) - (u_L - u_R) = 0$ of the one-dimensional problem. The function $\Theta(p)$ increases monotonically, so the root is unique if it exists.

**Graphical meaning.** Each stream corresponds to a curve in the $(\theta, p)$ plane — all states into which it can be brought by a single wave. The upper-stream curve $\theta = \theta_R + f_R(p)$ passes through the point $(\theta_R, p_R)$ and goes up and to the right: the higher the pressure, the larger the angle. The lower-stream curve $\theta = \theta_L - f_L(p)$ passes through $(\theta_L, p_L)$ and goes up and to the left. On each curve the shock branch lies above the initial point and the rarefaction branch below it. The intersection of the curves is the sought $(\theta^*, p^*)$.

```
   p
   ^
   |            \           /
   |              \       /
   |                \   /
 p*|- - - - - - - - - X            intersection point (θ*, p*)
   |                /   \
   |              /       \
   |      W_R → o           o ← W_L
   |          /               \
   +----------------------------------------> θ
                θ_R   θ*    θ_L
```

In the figure the streams converge ($\theta_L > \theta_R$) at equal pressures. The curves intersect above both initial points, so both waves are shocks. The relative position of the points $W_L$ and $W_R$ immediately shows which waves will arise:

- streams converge ($\theta_L > \theta_R$), pressures close — two shocks;
- streams diverge ($\theta_L < \theta_R$), pressures close — two fans;
- pressures differ strongly — the intersection lies between $p_L$ and $p_R$: a shock in the lower-pressure stream and a fan in the higher-pressure one.

Formally the wave types are determined by the sign of $\Theta$ at $p_{\min} = \min(p_L, p_R)$ and $p_{\max} = \max(p_L, p_R)$: if $\Theta(p_{\min}) > 0$ the root lies below both pressures (two fans), if $\Theta(p_{\max}) < 0$ above both (two shocks), otherwise between them (a shock and a fan). Knowing the wave types in advance is not necessary for the solution: the function $f_k$ itself selects the proper branch by the sign of $p - p_k$.

**Initial guess.** Replace both functions by their linear approximations $f_k(p) \approx Z_k\,(p - p_k)$. The equation becomes linear, and its root is

$$
p^{(0)} = \frac{\theta_L - \theta_R + Z_L p_L + Z_R p_R}{Z_L + Z_R}.
$$

This is the solution of the problem in the acoustic approximation. If $p^{(0)} \le 0$ (strong rarefaction), one takes, for example, $p^{(0)} = 0.1\min(p_L, p_R)$.

**Iterations.** The root is refined by Newton's method:

$$
p^{(i+1)} = p^{(i)} - \frac{\Theta\bigl(p^{(i)}\bigr)}{\Theta'\bigl(p^{(i)}\bigr)},\qquad \Theta' = f_L' + f_R'.
$$

The derivatives are computed as follows:

- on the rarefaction branch — analytically: $f_k'(p) = \dfrac{\sqrt{M^2 - 1}}{\gamma\,p\,M^2}$, where $M = M_k(p)$. This is the same formula as for $Z_k$, but with the current $p$ and $M$;
- on the shock branch it is simpler to compute it numerically, e.g. $f_k'(p) \approx \dfrac{f_k(p + \varepsilon p) - f_k(p)}{\varepsilon p}$ with $\varepsilon \sim 10^{-6}$, or to use the secant method instead of Newton's.

Iterations stop when $|p^{(i+1)} - p^{(i)}|/p^{(i)} < 10^{-8}$; usually 2–4 iterations suffice. Since $\Theta$ is monotone, the iterations are easy to safeguard: if the next approximation falls outside the interval $(0,\ p_{\max}^{\text{sh}}]$, where $p_{\max}^{\text{sh}}$ is the weak-branch limit from Section 5.7, the step is replaced by bisection.

In regions of smooth flow the discontinuities between neighbouring cells are small, and the linear solution $p^* \approx p^{(0)}$ itself is often sufficient. The exact solution is needed where strong shocks reach the face.

### 5.5. Parameters Behind the Waves and Wave Positions

The pressure $p^*$ has been found. The rest is computed from explicit formulas.

**Slip line angle:**

$$
\theta^* = \theta_L - f_L(p^*) = \theta_R + f_R(p^*).
$$

Both formulas must give the same number; their agreement is a convenient check of iteration convergence.

**States $W_{L*}$ and $W_{R*}$.** For each $k = L, R$ compute $r_k = p^*/p_k$ and then:

- density behind a shock ($r_k > 1$): $\rho_{k*} = \rho_k \dfrac{(\gamma + 1) r_k + \gamma - 1}{(\gamma - 1) r_k + \gamma + 1}$;

- density behind a fan ($r_k < 1$): $\rho_{k*} = \rho_k\, r_k^{1/\gamma}$;

- the flow speed in both cases follows from total enthalpy conservation: $q_{k*}^2 = 2H_k - \dfrac{2\gamma}{\gamma - 1}\dfrac{p^*}{\rho_{k*}}$;

- velocity components: $u_{k*} = q_{k*}\cos\theta^*$, $\ v_{k*} = q_{k*}\sin\theta^*$;

- Mach number and Mach angle (needed for fan slopes): $M_{k*}^2 = \dfrac{\rho_{k*}\, q_{k*}^2}{\gamma p^*}$, $\ \mu_{k*} = \arcsin\dfrac{1}{M_{k*}}$.

**Ray slopes.** The solution regions are separated by rays with the following slopes:

| Element | If shock | If rarefaction fan |
|---|---|---|
| wave $C^-$ (lower stream) | $\xi_L = \tan(\theta_L - \beta_L)$ | from $\xi_L^{\text{h}} = \tan(\theta_L - \mu_L)$ to $\xi_L^{\text{t}} = \tan(\theta^* - \mu_{L*})$ |
| slip line | $\xi_0 = \tan\theta^*$ | $\xi_0 = \tan\theta^*$ |
| wave $C^+$ (upper stream) | $\xi_R = \tan(\theta_R + \beta_R)$ | from $\xi_R^{\text{t}} = \tan(\theta^* + \mu_{R*})$ to $\xi_R^{\text{h}} = \tan(\theta_R + \mu_R)$ |

The superscripts "h" and "t" denote the head and tail of the fan. The fan head is the Mach line of the undisturbed stream $W_k$, the tail is the Mach line of the state $W_{k*}$. The angle $\beta_k$ is taken from the formula of Section 5.2 with $r = r_k$.

As long as $u > c$ in all four states, all angles in the table lie in the interval $(-\pi/2,\ \pi/2)$, the tangents are finite, and the rays are ordered from bottom to top: wave $C^-$, slip line, wave $C^+$.

The extreme slopes are needed for the step choice (Section 7): $\xi^-$ is $\xi_L$ or $\xi_L^{\text{h}}$, and $\xi^+$ is $\xi_R$ or $\xi_R^{\text{h}}$.

### 5.6. State on a Cell Face

A cell face is a ray with slope $s$ emanating from the same point. The solution on a ray is constant, so it only remains to determine which zone this ray falls into.

```
              W_L           W_L*            W_R*           W_R
   ξ  ----------------|----------------|----------------|---------------->
                     ξ_L              ξ_0              ξ_R
                   wave C-          slip line         wave C+
```

If a wave is a shock, its boundary is a single point on the $\xi$ axis. If it is a fan, the point is replaced by an interval between head and tail, inside which the state varies smoothly.

State selection:

- $s$ below the $C^-$ wave (less than $\xi_L$ or $\xi_L^{\text{h}}$) — $W_L$;
- $s$ between the $C^-$ wave and the slip line — $W_{L*}$;
- $s$ between the slip line and the $C^+$ wave — $W_{R*}$;
- $s$ above the $C^+$ wave (greater than $\xi_R$ or $\xi_R^{\text{h}}$) — $W_R$;
- $s$ inside a fan — the state in the fan on that ray.

**State inside a fan.** We need the Mach number $M$ and flow angle $\theta$ on the ray with slope $s$. There are two equations.

For the lower fan ($C^-$ family): a fan ray is a $C^-$ characteristic, i.e. $\theta - \mu(M) = \arctan s$, and throughout the fan the invariant $\theta - \omega(M) = \theta_L - \omega(M_L)$ is conserved (Section 5.3). Eliminating $\theta$, we obtain a single equation for $M$:

$$
\omega(M) - \mu(M) = \arctan s - \theta_L + \omega(M_L),\qquad \theta = \arctan s + \mu(M).
$$

For the upper fan ($C^+$ family) similarly: $\theta + \mu(M) = \arctan s$ and $\theta + \omega(M) = \theta_R + \omega(M_R)$, whence

$$
\omega(M) - \mu(M) = \theta_R + \omega(M_R) - \arctan s,\qquad \theta = \arctan s - \mu(M).
$$

The left-hand side $\omega(M) - \mu(M)$ increases monotonically with $M$, so the root is unique. It lies between $M_k$ (at the fan head) and $M_{k*}$ (at the tail) and is easily found by Newton's method with the derivatives

$$
\frac{d\omega}{dM} = \frac{\sqrt{M^2 - 1}}{M\bigl(1 + \frac{\gamma - 1}{2}M^2\bigr)},\qquad
\frac{d\mu}{dM} = -\frac{1}{M\sqrt{M^2 - 1}}.
$$

From the found $M$ the other parameters are obtained along the isentrope from the state $W_k$ ($k = L$ for the lower fan, $k = R$ for the upper one):

$$
p = p_k\left[\frac{1 + \frac{\gamma - 1}{2}M_k^2}{1 + \frac{\gamma - 1}{2}M^2}\right]^{\frac{\gamma}{\gamma - 1}},\qquad
\rho = \rho_k\Bigl(\frac{p}{p_k}\Bigr)^{1/\gamma},\qquad
q = M\sqrt{\frac{\gamma p}{\rho}},
$$

and $u = q\cos\theta$, $v = q\sin\theta$.

**Remarks.**

- If the face coincides with the slip line ($s = \tan\theta^*$), the mass flow through it is zero, and the flux $\boldsymbol\Phi = (0,\ -s p^*,\ p^*,\ 0)^T$ is the same on both sides of the discontinuity. Therefore it does not matter which of the states $W_{L*}$, $W_{R*}$ is chosen in this case.
- Interior grid lines usually run almost along streamlines, so most often the face falls into the $W_{L*}$ or $W_{R*}$ zone. But the case of a fan crossing the face (for example, near a contour corner) must be handled.

### 5.7. When There Is No Solution

The equation $\Theta(p) = 0$ may have no solution in the class of supersonic flows. Three situations are possible.

1. **Shock detachment.** The weak shock branch is bounded from above, so the pressure cannot exceed $p_{\max}^{\text{sh}} = \min\bigl(p_L\, r_{\max}(M_L),\ p_R\, r_{\max}(M_R)\bigr)$. If $\Theta(p_{\max}^{\text{sh}}) < 0$, the streams converge at too large an angle: the required deflection exceeds $\delta_{\max}$. This is the same phenomenon as shock detachment from a wedge with too large an angle. In a real flow a Mach reflection or a detached shock with a subsonic region arises, and the marching method is not applicable.

2. **Vacuum formation.** For strongly diverging streams the pressure drops to zero before the angles equalize: $\Theta(0) > 0$. A flow cannot be turned by a rarefaction fan by an arbitrarily large angle: on expansion to zero pressure $M \to \infty$, while $\omega$ tends to a finite limit $\omega_{\max} = \frac{\pi}{2}\bigl(\sqrt{(\gamma+1)/(\gamma-1)} - 1\bigr)$ (for $\gamma = 1.4$ this is $130.45^\circ$). The vacuum condition is $\theta_R - \theta_L > 2\omega_{\max} - \omega(M_L) - \omega(M_R)$. An empty zone forms between the streams. In practice this happens only at the boundary of a jet exhausting into vacuum.

3. **Subsonic state behind a strong shock.** A solution exists, but the condition $u > c$ is violated in $W_{L*}$ or $W_{R*}$. For a single face this is acceptable, but if such a zone grows, a cell will get $D < 0$ after the update (Section 4).

In the program the first two cases are detected before the iterations from the signs of $\Theta(p_{\max}^{\text{sh}})$ and $\Theta(0)$, the third by checking $u > c$ in the obtained states. In all three cases the computation should be stopped with a diagnostic message: this is not a solver error but a departure from the method's applicability limits (Section 9).

### 5.8. Computational Procedure

**Computational procedure.** Input: states $W_L$, $W_R$ and face slope $s$. Output: the state on the face and the slopes of the outer waves $\xi^-$, $\xi^+$.

1. For $k = L, R$ compute $q_k$, $\theta_k$, $M_k$, $H_k$, $Z_k$, $\omega(M_k)$, $r_{\max}(M_k)$.
2. Check the existence of a solution: $\Theta(p_{\max}^{\text{sh}}) \ge 0$ and $\Theta(0) < 0$ (Section 5.7).
3. Find the initial guess $p^{(0)}$ and refine it by iterations to $p^*$ (Section 5.4).
4. Compute $\theta^*$, then $\rho_{k*}$, $q_{k*}$, $M_{k*}$ for both streams (Section 5.5).
5. Compute the slopes of the waves and the slip line (table of Section 5.5).
6. From the face slope $s$ select the zone and return the state in it (Section 5.6).

Steps 1–5 do not depend on the face slope. Therefore in the overall algorithm (Section 8) the problem is first solved on all faces, the step $\Delta x_n$ is chosen from the wave slopes, and only then, when the face slopes are known, step 6 is performed.


## 6. Initial and Boundary Conditions

### 6.1. Wall and Symmetry Plane

Let the lower boundary be an impermeable wall $y = y_w(x)$ with slope $s_w = \bigl(y_w(x_{n+1}) - y_w(x_n)\bigr)/\Delta x_n$ over the step and angle $\theta_w = \arctan s_w$. The gas of cell $1$ lies above the wall, so a single $C^+$ wave arises, turning the flow along the wall:

$$
f_1(p_w) = \theta_w - \theta_1.
$$

If the wall is turned into the flow ($\theta_w > \theta_1$), this is a shock, otherwise a rarefaction wave. Since $v - s_w u = 0$ at the wall, the flux through it contains only pressure:

$$
\boldsymbol\Phi_{1/2} = (0,\ -s_w p_w,\ p_w,\ 0)^T.
$$

The same result is obtained by introducing a ghost cell with the mirrored state $(p_1, \rho_1, q_1, 2\theta_w - \theta_1)$ and solving the ordinary breakup problem: by symmetry the slip line coincides with the wall. For an upper wall (channel) the wave belongs to the $C^-$ family, and the equation becomes $f_N(p_w) = \theta_N - \theta_w$.

A symmetry plane is a wall with $s_w = 0$.

It is convenient to place contour corner points exactly on layers $x_n$, reducing the step for this purpose: then the wall is straight on every step, and a corner generates the exact solution of the flow-turning problem.

### 6.2. Outer Boundary

Two approaches are possible.

**Boundary in the undisturbed flow (bow shock capturing).** The upper boundary is placed well above all disturbances, for example as a ray with slope greater than that of the expected bow shock. The undisturbed-flow flux is computed through it: $\boldsymbol\Phi = \mathbf G(W_\infty) - s\,\mathbf F(W_\infty)$. This approach is simple, but some cells are wasted on uniform flow, and the bow shock is smeared over several cells.

**Bow shock fitting.** The upper boundary is the shock itself. On each step the breakup problem between cell $N$ (below) and the undisturbed flow $W_\infty$ (above) is solved. Its $C^+$ wave is the bow shock: its slope is $s_e = \tan(\theta_\infty + \beta_\infty)$, and the new boundary position is $y^{n+1}_e = y^n_e + \Delta x_n\, s_e$. Since the face moves with the shock, by the Rankine–Hugoniot relations the flux through it is the same on both sides of the shock, and it is convenient to compute it from $W_\infty$. If $p^* \le p_\infty$, the shock has degenerated, and the boundary follows the characteristic with slope $\tan(\theta_\infty + \mu_\infty)$.

Bow shock fitting reduces the number of cells and gives accurate parameters behind the shock. Interior shocks are still captured.

### 6.3. Free Jet Boundary

At the boundary of a jet exhausting into a medium with pressure $p_a$, the pressure is prescribed, and the boundary itself is a streamline. The gas lies below the boundary, so the wave belongs to the $C^-$ family:

$$
\theta_b = \theta_N - f_N(p_a),\qquad s_e = \tan\theta_b,\qquad \boldsymbol\Phi_{N+1/2} = (0,\ -s_e p_a,\ p_a,\ 0)^T.
$$

The boundary moves according to $y^{n+1}_e = y^n_e + \Delta x_n\, s_e$.

### 6.4. Initial Data

On the initial layer $x_0$ a distribution of parameters supersonic in $x$ ($u > c$) is needed. Its sources are:

- for a pointed body — a uniform free stream or the exact wedge solution near the apex;
- for a blunt body — the solution of the nose-region problem (usually by time-marching to a steady state) on a line downstream of the sonic line, where $u > c$ everywhere;
- for a nozzle — the solution in the transonic region on a line where $u > c$.

## 7. Choice of the Marching Step

The step $\Delta x_n$ is limited by the same requirement as the time step in Godunov's scheme: the waves generated by the discontinuity breakups at faces $j - 1/2$ and $j + 1/2$ must not reach the opposite face of cell $j$ within the step. Then the state on each face throughout the step is determined only by "its own" breakup problem, and the flux $\boldsymbol\Phi$ is computed exactly. For cell $j$ this gives

$$
\Delta x_n\bigl(\xi^{+}_{j-1/2} - s_{j+1/2}\bigr) \le h^n_j,\qquad
\Delta x_n\bigl(s_{j-1/2} - \xi^{-}_{j+1/2}\bigr) \le h^n_j,
$$

where $\xi^{+}_{j-1/2}$ is the slope of the upper wave (shock or fan head) from the problem at face $j - 1/2$, and $\xi^{-}_{j+1/2}$ is the slope of the lower wave from the problem at face $j + 1/2$.

The breakup solutions at the points of layer $x_n$ do not depend on $\Delta x_n$, so they can be found before choosing the step, and the exact wave slopes can be substituted into these conditions. The face slopes $s$ are then estimated from the slopes of the domain boundaries on layer $x_n$. A simpler and usually sufficient option is an a priori estimate from the Mach line slopes in the cells:

$$
\Delta x_n = \kappa\,\min_j \frac{h^n_j}{\lambda^+_j - \lambda^-_j},\qquad
\lambda^+ - \lambda^- = \frac{2c\sqrt{q^2 - c^2}}{u^2 - c^2},\qquad 0 < \kappa \le 1.
$$

Faces are usually close to streamlines ($s \approx \tan\theta$), and waves depart from them on both sides by roughly half the width of the characteristic fan. Therefore for $\kappa \le 1$ there remains roughly a twofold margin, which also covers the fact that shocks are steeper than the characteristics ahead of them. Typically $\kappa = 0.5$–$0.9$. If grid lines deviate noticeably from streamlines, it is safer to check the exact condition.

As $u \to c$, the difference $\lambda^+ - \lambda^-$ grows without bound and the step tends to zero. A sharp decrease of the step is a sign that the computation is approaching the method's applicability limit.

## 8. Algorithm of a Single Step

Given: nodes $y^n_{j+1/2}$ and parameters $W^n_j = (p, \rho, u, v)^n_j$, $j = 1, \dots, N$.

1. From $W^n_j$ compute $\mathbf F^n_j$, $H_j$, $M_j$, $\theta_j$, $\lambda^\pm_j$.
2. Solve the breakup problems at interior faces (pairs $W_j$, $W_{j+1}$) and at the outer boundary if a shock is fitted there or a free boundary is prescribed. Store $p^*$, $\theta^*$, the states $W_{L*}$, $W_{R*}$ and the wave slopes.
3. Choose the step $\Delta x_n$ (Section 7). If necessary, reduce it to land on a contour corner point or the final section.
4. Find the boundary positions on layer $x_{n+1}$: $y_w(x_{n+1})$ from the geometry, $y_e^{n+1}$ from the slope of the shock or free boundary, or from a prescribed line. Build the nodes $y^{n+1}_{j+1/2}$ and compute $s_{j+1/2}$ and $h^{n+1}_j$.
5. For each interior face select the state on the ray $\xi = s_{j+1/2}$ (Section 5.6) and compute $\boldsymbol\Phi_{j+1/2}$.
6. Compute the fluxes through the boundary faces (Section 6).
7. Find $\mathbf F^{n+1}_j$ by the scheme of Section 3.2.
8. Recover $W^{n+1}_j$ (Section 4) and check that $D \ge 0$.
9. Set $x_{n+1} = x_n + \Delta x_n$ and proceed to the next step.

Memory and work per step are proportional to $N$. The number of steps is of order $L/\Delta x \sim N L / Y$, where $L$ is the domain length and $Y$ its transverse size. The whole computation costs about as much as a few steps of a time-marching method on the same grid.

## 9. Limitations

- $u > c$ must hold everywhere. Subsonic zones (behind detached shocks, in Mach reflection, at the nose of a blunt body) and reversed flows cannot be computed by the marching method. Such problems are solved in a combined way: the subsonic and transonic part by time-marching, the supersonic part by marching in $x$.
- The model is inviscid: boundary layers, separation and mixing are not described.
- For a perfect gas the formulas of Sections 4 and 5 are written explicitly. For a gas with an arbitrary equation of state only the conversion from $\mathbf F$ to gas parameters and the breakup solution change; the structure of the scheme remains the same.

## Implementation

Code implementing the described method for computing an underexpanded jet exhausting from a nozzle is available here: <https://github.com/Semen1996/godunov-marshing-method>.

## References

1. Godunov S. K. A difference method for numerical calculation of discontinuous solutions of the equations of hydrodynamics // Matematicheskii Sbornik. 1959. Vol. 47 (89), No. 3. P. 271–306 (in Russian).
2. Godunov S. K., Zabrodin A. V., Ivanov M. Ya., Kraiko A. N., Prokopov G. P. Numerical Solution of Multidimensional Problems of Gas Dynamics. Moscow: Nauka, 1976 (in Russian).
3. Ivanov M. Ya., Kraiko A. N., Mikhailov N. V. A shock-capturing method for two-dimensional and three-dimensional supersonic flows. I // USSR Computational Mathematics and Mathematical Physics. 1972. Vol. 12, No. 2.
4. Kolgan V. P. Application of the principle of minimal values of the derivative to the construction of finite-difference schemes for calculating discontinuous solutions of gas dynamics // Uchenye Zapiski TsAGI. 1972. Vol. 3, No. 6. P. 68–77 (in Russian).
5. Rodionov A. V. A monotone second-order scheme for shock-capturing computation of non-equilibrium flows // USSR Computational Mathematics and Mathematical Physics. 1987. Vol. 27, No. 4. P. 585–593.
6. Loh C. Y., Hui W. H. A new Lagrangian method for steady supersonic flow computation // Journal of Computational Physics. 1990. Vol. 89.
7. Anderson J. D. Modern Compressible Flow: With Historical Perspective. 3rd ed. McGraw-Hill, 2003.
8. Toro E. F. Riemann Solvers and Numerical Methods for Fluid Dynamics. 3rd ed. Springer, 2009.
