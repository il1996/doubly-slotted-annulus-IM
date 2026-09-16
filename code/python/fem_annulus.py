# -*- coding: utf-8 -*-
"""
Scalar-potential finite-element solution on the bare gap annulus (formulation B).

Boundary conditions on the two circles are set arc by arc:
  * tooth-face arcs : Dirichlet, phi = U_j (stator) or U_r (rotor)
  * opening arcs    : either homogeneous Neumann (the manuscript's Phi_O = 0,
                      no radial flux) or Dirichlet as well (a smooth comparator
                      where the opening arcs carry a prescribed potential)
The fundamental of B_r at mid-gap is obtained from the exact nodal fluxes
(reactions) on the stator Dirichlet nodes: the flux through each tooth face
is the sum of the reactions on its nodes, and the p-th harmonic of the
tooth-flux staircase is formed.  For the smooth annulus the same functional is
evaluated analytically, so that the ratio is consistent.
"""
import numpy as np, gmsh, time
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from dsop import Machine, MU0
from fem_slots import assemble_P1


def build_annulus_mesh(M, phi_rot, lc, stator_open=True, rotor_open=True, verbose=False):
    """Annulus R_r < r < R_s with the face and opening arcs as separate
    boundary curves.  Returns nodes, triangles and, for each boundary node,
    the arc id (stator: 0..2Ns-1 alternating face/opening, rotor likewise)."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
    gmsh.model.add("annulus")
    geo = gmsh.model.geo
    Rs = M.Rs * 1e3; Rr = M.Rr * 1e3
    centre = geo.addPoint(0, 0, 0)

    def circle_arcs(R, N, b0, phi, tag0):
        """points at the face/opening junctions, arcs in between; returns arc tags
        and a list telling whether each arc is a face (True) or an opening."""
        tau = 2 * np.pi / N; a_o = b0 / R; a_f = tau - a_o
        angs = []
        for k in range(N):
            c = tau * k + phi
            angs.append(c - a_f / 2); angs.append(c + a_f / 2)
        pts = [geo.addPoint(R * np.cos(a), R * np.sin(a), 0) for a in angs]
        arcs = []; isface = []
        for i in range(2 * N):
            p1 = pts[i]; p2 = pts[(i + 1) % (2 * N)]
            # arc from p1 to p2 (counter-clockwise); split long arcs are fine (< pi)
            arcs.append(geo.addCircleArc(p1, centre, p2))
            isface.append(i % 2 == 0)
        return arcs, isface

    bs0 = M.bs0 * 1e3 if stator_open else 0.0
    br0 = M.br0 * 1e3 if rotor_open else 0.0
    arcs_s, face_s = circle_arcs(Rs, M.Ns, bs0, 0.0, 0)
    arcs_r, face_r = circle_arcs(Rr, M.Nr, br0, phi_rot, 0)
    lo = geo.addCurveLoop(arcs_s)
    li = geo.addCurveLoop(arcs_r)
    surf = geo.addPlaneSurface([lo, li])
    geo.synchronize()
    for i, a in enumerate(arcs_s):
        gmsh.model.addPhysicalGroup(1, [a], tag=1000 + i)
    for i, a in enumerate(arcs_r):
        gmsh.model.addPhysicalGroup(1, [a], tag=2000 + i)
    gmsh.model.addPhysicalGroup(2, [surf], tag=1)
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    ntags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)[:, :2] * 1e-3
    idx = {t: i for i, t in enumerate(ntags)}
    etypes, etags, enodes = gmsh.model.mesh.getElements(2)
    tri = np.vstack([np.array([idx[t] for t in en]).reshape(-1, 3) for et, en in zip(etypes, enodes) if et == 2])
    # boundary nodes per arc
    arc_nodes = {}
    for i, a in enumerate(arcs_s):
        nt, _, _ = gmsh.model.mesh.getNodes(1, a, includeBoundary=True)
        arc_nodes[('s', i)] = np.array([idx[t] for t in nt])
    for i, a in enumerate(arcs_r):
        nt, _, _ = gmsh.model.mesh.getNodes(1, a, includeBoundary=True)
        arc_nodes[('r', i)] = np.array([idx[t] for t in nt])
    gmsh.finalize()
    # drop orphan nodes (the centre point is meshed as a node but belongs to no triangle)
    used = np.unique(tri.ravel())
    remap = -np.ones(len(xyz), dtype=int); remap[used] = np.arange(len(used))
    xyz = xyz[used]; tri = remap[tri]
    arc_nodes = {k: remap[v][remap[v] >= 0] for k, v in arc_nodes.items()}
    return xyz, tri, arc_nodes, face_s, face_r


def solve_scalar(M, xyz, tri, arc_nodes, face_s, face_r, Us, Ur=0.0, opening='neumann',
                 Uopen_s=None, Uopen_r=None):
    """Solve -div(mu0 grad phi) = 0 with Dirichlet data on the faces (and on
    the openings when opening='dirichlet').  Returns phi and the nodal
    reactions R = K phi (flux leaving the domain at Dirichlet nodes)."""
    K, area = assemble_P1(xyz, tri, np.full(len(tri), MU0))
    n = len(xyz)
    dir_val = np.full(n, np.nan)
    for i in range(2 * M.Ns):
        nodes = arc_nodes[('s', i)]
        if face_s[i]:
            dir_val[nodes] = Us[i // 2]
        elif opening == 'dirichlet':
            dir_val[nodes] = Uopen_s[i // 2] if Uopen_s is not None else 0.0
    for i in range(2 * M.Nr):
        nodes = arc_nodes[('r', i)]
        if face_r[i]:
            dir_val[nodes] = Ur
        elif opening == 'dirichlet':
            dir_val[nodes] = Uopen_r[i // 2] if Uopen_r is not None else Ur
    # junction nodes are shared between a face arc and an opening arc; the face value wins
    for i in range(2 * M.Ns):
        if face_s[i]:
            dir_val[arc_nodes[('s', i)]] = Us[i // 2]
    for i in range(2 * M.Nr):
        if face_r[i]:
            dir_val[arc_nodes[('r', i)]] = Ur
    isdir = ~np.isnan(dir_val)
    free = np.where(~isdir)[0]; fixed = np.where(isdir)[0]
    phi = np.zeros(n); phi[fixed] = dir_val[fixed]
    Kff = K[free][:, free].tocsc(); Kfd = K[free][:, fixed]
    phi[free] = spla.spsolve(Kff, -Kfd @ phi[fixed])
    R = K @ phi                       # reactions: nonzero only at Dirichlet nodes
    return phi, R


def tooth_flux_fundamental(M, R, arc_nodes, face_s, per_length=True):
    """Flux entering each stator tooth face (from the reactions), and the
    p-th harmonic of the tooth-flux density staircase B_j = Phi_j/(w_j R_s)."""
    Phi = np.zeros(M.Ns)
    for i in range(2 * M.Ns):
        if face_s[i]:
            Phi[i // 2] += R[arc_nodes[('s', i)]].sum()
    # junction nodes belong to two face arcs? no: a junction joins a face and an opening; fine.
    tau = 2 * np.pi / M.Ns
    th = tau * np.arange(M.Ns)
    # tooth flux per unit length -> equivalent mean B over the pitch
    Bj = Phi / (tau * M.Rs)
    c = (2 / M.Ns) * np.sum(Bj * np.exp(-1j * M.p * th))
    return Phi, abs(c)


def smooth_tooth_fundamental_analytic(M, Us_teeth):
    """Same functional for the smooth annulus with the stator potential
    piecewise constant per tooth pitch (staircase, jumps at the slot centres):
    tooth flux = integral over the pitch of B_r(R_s).  Exact per harmonic."""
    tau = 2 * np.pi / M.Ns; th = tau * np.arange(M.Ns)
    # Fourier of the staircase u(theta) = Us_j on pitch j
    Nh = 4096
    n = np.arange(1, Nh + 1)
    # coefficients (cos, sin) of the staircase
    ec = np.exp(-1j * np.outer(n, th))
    sinc = np.sin(n * tau / 2) / (n * tau / 2)
    coef = (tau / np.pi) * sinc * (ec @ Us_teeth)          # a_n - i b_n
    # B_r at R_s: -mu0 dphi/dr = -(mu0/R_s) n [u_s cosh(nX) - u_r] / sinh(nX), u_r = 0
    Bn = -(MU0 / M.Rs) * n * np.cosh(n * M.X) / np.sinh(n * M.X) * coef
    # tooth flux density (mean over pitch j) = sum_n Bn * sinc_n * e^{i n th_j}  (real part)
    Bj = np.real((Bn * sinc) @ np.exp(1j * np.outer(n, th)))
    Phi = Bj * tau * M.Rs
    c = (2 / M.Ns) * np.sum(Bj * np.exp(-1j * M.p * th))
    return Phi, abs(c)


if __name__ == '__main__':
    import sys
    M = Machine()
    lc = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
    i3 = np.sqrt(2) * np.array([1.0, -0.5, -0.5]); Us = M.tooth_mmf(i3)
    t = time.time()
    # (1) check: smooth annulus, openings as Dirichlet at the tooth potential -> staircase
    xyz, tri, arcs, fs, fr = build_annulus_mesh(M, 0.0, lc)
    # opening arc i (odd) lies between tooth i//2 and tooth i//2+1: put half the jump? For the check,
    # use a genuinely smooth comparator: no openings at all (bs0 = 0) is impossible geometrically here,
    # so give the opening the potential of the tooth on its left (staircase shifted by half an opening).
    phi, R = solve_scalar(M, xyz, tri, arcs, fs, fr, Us, 0.0, opening='dirichlet',
                          Uopen_s=Us, Uopen_r=np.zeros(M.Nr))
    Phi, B1 = tooth_flux_fundamental(M, R, arcs, fs)
    Phi_a, B1_a = smooth_tooth_fundamental_analytic(M, Us)
    print(f"smooth check (Dirichlet openings): nodes {len(xyz)}  B1 FEM {B1:.6f}  analytic {B1_a:.6f}  ratio {B1/B1_a:.5f}  ({time.time()-t:.0f}s)")
    # (2) Neumann openings
    phi, R = solve_scalar(M, xyz, tri, arcs, fs, fr, Us, 0.0, opening='neumann')
    Phi, B1n = tooth_flux_fundamental(M, R, arcs, fs)
    print(f"Neumann openings: B1 {B1n:.6f}  kC(vs analytic smooth) = {B1_a/B1n:.5f}")


# ---------------------------------------------------------------- mid-gap
def midgap_harmonics_scalar(M, xyz, tri, phi, nsamp=8192, dv_frac=1.0 / 6.0):
    """|B_{r,n}| at mid-gap from the scalar potential sampled on two circles
    v_m +- dv (exact per harmonic for the annulus)."""
    from matplotlib.tri import Triangulation, LinearTriInterpolator
    T = Triangulation(xyz[:, 0], xyz[:, 1], tri); it = LinearTriInterpolator(T, phi)
    th = 2 * np.pi * np.arange(nsamp) / nsamp
    vm = np.log(M.Rm / M.Rr); dv = dv_frac * M.X
    out = []
    for v in (vm - dv, vm + dv):
        r = M.Rr * np.exp(v)
        a = np.asarray(it(r * np.cos(th), r * np.sin(th)).filled(np.nan))
        if np.isnan(a).any():
            raise RuntimeError('sampling outside mesh')
        out.append(np.fft.rfft(a) * 2.0 / nsamp)        # a_n - i b_n
    n = np.arange(len(out[0]))
    beta = (out[1] - out[0]) / (2 * np.sinh(np.maximum(n, 1) * dv))
    Bn = (MU0 * n / M.Rm) * np.abs(beta)
    return Bn


def piecewise_constant_fourier(edges, values, n):
    """(a_n - i b_n) for a function equal to values[k] on [edges[k], edges[k+1])."""
    e0 = np.asarray(edges); e1 = np.roll(e0, -1); e1 = np.where(e1 <= e0, e1 + 2 * np.pi, e1)
    n = np.asarray(n, float)[:, None]
    integ = (np.exp(-1j * n * e0[None, :]) - np.exp(-1j * n * e1[None, :])) / (1j * n)
    return (integ @ np.asarray(values)) / np.pi


def annulus_Bn_from_surface_coeffs(M, cs, cr, n, r=None):
    """|B_{r,n}| at radius r for stator/rotor surface Fourier coefficients cs, cr."""
    if r is None:
        r = M.Rm
    v = np.log(r / M.Rr); X = M.X
    a = n * np.cosh(n * v) / np.sinh(n * X); b = n * np.cosh(n * (X - v)) / np.sinh(n * X)
    return (MU0 / r) * np.abs(a * cs - b * cr)
