# -*- coding: utf-8 -*-
"""
Slot-cavity admittance element.

The annulus operator of the manuscript closes the slot openings with the
condition Phi_O = 0 (no flux crosses the bore over an opening).  Physically
the opening leads into a cavity bounded by iron at the potentials of the two
adjacent teeth, and the flux that crosses the opening terminates on those
walls.  This module computes, by a small two-dimensional finite-element
solution in the local frame of one slot, the admittance matrix of that
cavity:

    Phi = Q U,   U = [U_mouth(1..nO), U_left, U_right]

where the mouth potential is the piecewise-linear interpolation of the
opening-column values (the hat basis of the operator) between the two tooth
corners, whose potentials are those of the adjacent teeth.  Q is symmetric,
positive semi-definite, with null vector 1.

Stator cavity: isthmus + wedge, closed at the top of the conductor region by
the linear MMF ramp U_left -> U_right (uniform current density below; this is
the classical slot-leakage boundary and is exact for parallel walls).
Rotor cavity (no bar current): isthmus + pear, all walls at the rotor potential.
"""
import numpy as np, gmsh
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from fem_slots import assemble_P1, stator_slot_polygon, rotor_slot_polygon, BS0, HS0, BS1, HS1, BR0
from dsop import MU0


def _mesh_polygon(poly, lc_mouth, lc_far, mouth_y=0.0):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cav")
    geo = gmsh.model.geo
    pids = [geo.addPoint(x, y, 0) for (x, y) in poly]
    lids = [geo.addLine(pids[i], pids[(i + 1) % len(pids)]) for i in range(len(pids))]
    wl = geo.addCurveLoop(lids); surf = geo.addPlaneSurface([wl])
    geo.synchronize()
    f = gmsh.model.mesh.field
    fid = f.add("MathEval")
    f.setString(fid, "F", f"{lc_mouth} + {(lc_far - lc_mouth) / 2.0}*Abs(y - {mouth_y})")
    f.setAsBackgroundMesh(fid)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(2)
    ntags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)[:, :2]
    idx = {t: i for i, t in enumerate(ntags)}
    etypes, etags, enodes = gmsh.model.mesh.getElements(2)
    tri = np.vstack([np.array([idx[t] for t in en]).reshape(-1, 3) for et, en in zip(etypes, enodes) if et == 2])
    gmsh.finalize()
    used = np.unique(tri.ravel()); remap = -np.ones(len(xyz), int); remap[used] = np.arange(len(used))
    return xyz[used], remap[tri]


def cavity_element(kind, nO, L, lc_mouth=0.02, lc_far=0.6, basis='p1'):
    """Admittance matrix Q ((nO+2) x (nO+2)) of one slot cavity, in H per
    unit... i.e. flux (Wb) per ampere of potential, for the active length L.
    Column order: mouth columns 1..nO (from the left corner to the right),
    then U_left (tooth on the small-angle side), U_right."""
    if kind == 'stator':
        b0 = BS0
        # cavity polygon: mouth at y=0 (x in [-b0/2, b0/2]), isthmus, wedge, bottom at y = HS0+HS1
        poly = [(-b0 / 2, 0.0), (b0 / 2, 0.0), (b0 / 2, HS0), (BS1 / 2, HS0 + HS1),
                (-BS1 / 2, HS0 + HS1), (-b0 / 2, HS0)]
        ybot = HS0 + HS1
    else:
        b0 = BR0
        poly = rotor_slot_polygon(depth_extra=0.0)
        ybot = None
    xyz, tri = _mesh_polygon(poly, lc_mouth, lc_far)
    x, y = xyz[:, 0], xyz[:, 1]
    tol = 1e-7
    on_mouth = np.abs(y) < tol
    if kind == 'stator':
        on_bottom = np.abs(y - ybot) < tol
        on_left = (~on_mouth) & (~on_bottom) & (x < 0)
        on_right = (~on_mouth) & (~on_bottom) & (x > 0)
        # boundary nodes: those on any polygon edge. detect by distance to the edges
    # generic boundary detection: nodes on the polygon edges
    P = np.array(poly); nb = len(P)
    onb = np.zeros(len(xyz), bool)
    for i in range(nb):
        a = P[i]; b = P[(i + 1) % nb]; ab = b - a; l2 = ab @ ab
        t = np.clip(((xyz - a) @ ab) / l2, 0, 1)
        d = np.linalg.norm(xyz - (a + t[:, None] * ab), axis=1)
        onb |= d < 1e-6
    K, area = assemble_P1(xyz, tri, np.full(len(tri), MU0 * L * 1e-3))   # coordinates in mm -> K in H (scale-free in 2D, L in m: mu0*L)
    # NB: for the 2-D Laplacian the stiffness matrix is dimensionless in the
    # coordinates, so the mm units cancel; the factor mu0*L gives Wb/A.
    K = assemble_P1(xyz, tri, np.full(len(tri), MU0 * L))[0]
    n = len(xyz)
    # DOF groups on the boundary
    # mouth nodes sorted by x; mouth column centres
    mouth = np.where(on_mouth & onb)[0]
    xm = x[mouth]
    xc = -b0 / 2 + (np.arange(nO) + 0.5) * b0 / nO            # column centres
    xk = np.concatenate([[-b0 / 2], xc, [b0 / 2]])               # interpolation knots: corner, centres, corner
    ndof = nO + 2
    # boundary potential for unit value at dof j:
    def bc_vector(j):
        val = np.full(n, np.nan)
        # mouth: piecewise-linear interpolation between knots; knots values: corner_left = U_left (dof nO),
        # centres = U_k (dof k), corner_right = U_right (dof nO+1)
        kv = np.zeros(nO + 2)
        if j < nO:
            kv[j + 1] = 1.0
        elif j == nO:
            kv[0] = 1.0
        else:
            kv[-1] = 1.0
        if basis == 'p1':
            val[mouth] = np.interp(xm, xk, kv)
        else:   # p0: piecewise constant columns, corners belong to the walls
            col = np.clip(np.floor((xm + b0 / 2) / (b0 / nO)).astype(int), 0, nO - 1)
            v = np.zeros(len(mouth))
            if j < nO:
                v = (col == j).astype(float)
            val[mouth] = v
            # corner nodes exactly at the wall carry the wall potential
            cl = np.abs(xm + b0 / 2) < 1e-7; cr = np.abs(xm - b0 / 2) < 1e-7
            val[mouth[cl]] = 1.0 if j == nO else 0.0
            val[mouth[cr]] = 1.0 if j == nO + 1 else 0.0
        walls = onb & (~on_mouth)
        if kind == 'stator':
            left = walls & (x < 0) & (~on_bottom); right = walls & (x > 0) & (~on_bottom)
            val[left] = 1.0 if j == nO else 0.0
            val[right] = 1.0 if j == nO + 1 else 0.0
            bot = walls & on_bottom
            ramp = (x[bot] + BS1 / 2) / BS1              # 0 at the left wall, 1 at the right wall
            if j == nO:
                val[bot] = 1.0 - ramp
            elif j == nO + 1:
                val[bot] = ramp
            else:
                val[bot] = 0.0
        else:
            val[walls] = 1.0 if j == nO else 0.0       # rotor: all walls = "left" dof (single potential)
        return val
    fixed = np.where(onb)[0]; free = np.where(~onb)[0]
    Kff = K[free][:, free].tocsc(); Kfd = K[free][:, fixed]
    lu = spla.splu(Kff)
    Q = np.zeros((ndof, ndof))
    sols = []
    for j in range(ndof):
        val = bc_vector(j)
        u = np.zeros(n); u[fixed] = val[fixed]
        u[free] = lu.solve(-Kfd @ u[fixed])
        sols.append(u)
    # Q_ab = energy bilinear form = u_a^T K u_b
    for a in range(ndof):
        for b in range(ndof):
            Q[a, b] = sols[a] @ (K @ sols[b])
    Q = 0.5 * (Q + Q.T)
    if kind == 'rotor':
        # the two wall dofs are one and the same potential: merge them (right = left)
        Qm = np.zeros((nO + 2, nO + 2))
        Qm[:nO + 1, :nO + 1] = Q[:nO + 1, :nO + 1]
        # move the 'right' dof energy onto 'left' (it is the same wall)
        Qm[:nO + 1, nO] += Q[:nO + 1, nO + 1]; Qm[nO, :nO + 1] += Q[nO + 1, :nO + 1]; Qm[nO, nO] += Q[nO + 1, nO + 1]
        Q = Qm
    return Q, dict(nnodes=n, xyz=xyz, tri=tri)


if __name__ == '__main__':
    L = 164.782448048495e-3
    for kind in ['stator', 'rotor']:
        for nO in [2, 4, 8]:
            Q, info = cavity_element(kind, nO, L)
            print(kind, nO, "nodes", info['nnodes'], "null-vector residual", np.abs(Q @ np.ones(nO + 2)).max() / np.abs(Q).max(),
                  "eig min", np.linalg.eigvalsh(Q)[0])
            # total admittance mouth->walls: raise all mouth columns to 1, walls 0
            u = np.zeros(nO + 2); u[:nO] = 1.0
            print("   flux for unit mouth-wall potential difference: %.4e Wb/A  (parallel-plate estimate mu0 L b0/h: %.4e)" %
                  (u @ Q @ u, MU0 * L * (BS0 if kind == 'stator' else BR0) / (HS0 if kind == 'stator' else 1.0)))
