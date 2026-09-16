# -*- coding: utf-8 -*-
"""
Independent two-dimensional finite-element solution of the doubly slotted
air gap of the 18.5 kW machine, infinitely permeable iron, linear.

Two formulations are provided.

(A) Vector potential A_z on the real air domain (gap + stator slots + rotor
    slots).  Iron of infinite permeability is the boundary of the domain,
    on which H_t = 0, i.e. the homogeneous Neumann condition dA/dn = 0.
    The stator slot currents are the source.  This is the physical problem
    that Carter's factor approximates.

(B) Scalar potential phi on the bare gap annulus with the manuscript's
    opening condition: phi = U_j on the tooth faces (Dirichlet) and no
    radial flux over the openings (homogeneous Neumann).  This is the
    continuous problem the condensed operator of the manuscript
    discretises, and it gives the tiling limit that the operator sweeps
    approach.

Both return the fundamental (order p) of the radial flux density at
mid-gap, from which the slotting ratio k_C = B_g1(smooth)/B_g1(slotted) is
formed with the analytical smooth-annulus fundamental at the same MMF.

Geometry: Table 5 of the manuscript and the variables of the finite-element
project (IM_18kW_690V.aedt).  Units inside gmsh: millimetres.
"""
import numpy as np, gmsh, sys, time
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from matplotlib.tri import Triangulation, LinearTriInterpolator
from dsop import Machine, MU0

# ------------------------------------------------------------ slot shapes (mm)
# stator: semi-closed trapezoidal slot, RMxprt variables of the project
BS0, HS0 = 2.0, 0.5          # opening width, opening height
BS1, HS1 = 5.236, 2.5        # width at the end of the wedge, wedge height
BS2, HS2 = 8.5, 24.724       # width at the bottom, body height
# rotor: 2 mm x 1 mm isthmus, then two circles (5.804 and 2.038 mm) joined by tangents
BR0, HR0 = 2.0, 1.0
DR1, DR2, HR2 = 5.80425613768, 2.038213646116, 16.78981901905


def stator_slot_polygon(depth_extra=0.05):
    """Polygon of one stator slot in local coordinates (x tangential, y
    radial, y = 0 at the bore, y > 0 into the stator).  Counter-clockwise."""
    y0 = -depth_extra
    pts = [(-BS0 / 2, y0), (BS0 / 2, y0), (BS0 / 2, HS0), (BS1 / 2, HS0 + HS1),
           (BS2 / 2, HS0 + HS1 + HS2), (-BS2 / 2, HS0 + HS1 + HS2),
           (-BS1 / 2, HS0 + HS1), (-BS0 / 2, HS0)]
    return pts


def rotor_slot_polygon(depth_extra=0.05, narc=24):
    """Polygon of one rotor slot in local coordinates (y = 0 at the rotor
    surface, y > 0 into the rotor).  Isthmus + pear-shaped bar."""
    r1, r2 = DR1 / 2, DR2 / 2
    c1 = HR0 + r1                # centre of the upper circle
    c2 = c1 + HR2                # centre of the lower circle
    # tangent lines between the two circles: angle of the tangent
    d = c2 - c1
    alpha = np.arcsin((r1 - r2) / d)     # half-angle of the cone
    # right side, going down: from the isthmus to circle 1, then tangent, then circle 2
    pts = [(-BR0 / 2, -depth_extra), (BR0 / 2, -depth_extra), (BR0 / 2, HR0)]
    # circle 1 right part: from angle where x = BR0/2 (top) down to the tangent point
    th_top = np.arcsin((BR0 / 2) / r1)          # measured from the +y (upward) axis... use polar with y down
    # parametrise circle 1 by angle t from the upward vertical: point = (r1 sin t, c1 - r1 cos t)
    t_start = th_top
    t_tan = np.pi / 2 + alpha
    for t in np.linspace(t_start, t_tan, narc):
        pts.append((r1 * np.sin(t), c1 - r1 * np.cos(t)))
    # circle 2 from the tangent point round the bottom to the symmetric tangent point
    for t in np.linspace(t_tan, 2 * np.pi - t_tan, 2 * narc):
        pts.append((r2 * np.sin(t), c2 - r2 * np.cos(t)))
    # back up circle 1 on the left
    for t in np.linspace(2 * np.pi - t_tan, 2 * np.pi - t_start, narc):
        pts.append((r1 * np.sin(t), c1 - r1 * np.cos(t)))
    pts.append((-BR0 / 2, HR0))
    return pts


def place(pts, R, theta, inward):
    """Map local slot coordinates to the machine frame.  For a stator slot the
    slot goes outward (radius R + y); for a rotor slot inward (R - y)."""
    out = []
    er = (np.cos(theta), np.sin(theta)); et = (-np.sin(theta), np.cos(theta))
    for (x, y) in pts:
        r = R - y if inward else R + y
        out.append((r * er[0] + x * et[0], r * er[1] + x * et[1]))   # local Cartesian slot frame
    return out


def polygon_area(P):
    x = np.array([p[0] for p in P]); y = np.array([p[1] for p in P])
    return 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


# ------------------------------------------------------------ meshing (A)
def build_air_mesh(M, phi_rot, lc_gap=0.06, lc_far=2.0, stator_slots=True, rotor_slots=True,
                   verbose=False):
    """Mesh of the air domain (gap + slots).  Returns nodes (m), triangles,
    and per-element slot-body id (-1 if not in a stator slot body)."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
    gmsh.model.add("air")
    occ = gmsh.model.occ
    Rs = M.Rs * 1e3; Rr = M.Rr * 1e3
    disk_s = occ.addDisk(0, 0, 0, Rs, Rs)
    disk_r = occ.addDisk(0, 0, 0, Rr, Rr)
    tools_s = []
    if stator_slots:
        poly = stator_slot_polygon()
        for k in range(M.Ns):
            th = 2 * np.pi * (k + 0.5) / M.Ns
            P = place(poly, Rs, th, inward=False)
            pids = [occ.addPoint(x, y, 0) for (x, y) in P]
            lids = [occ.addLine(pids[i], pids[(i + 1) % len(pids)]) for i in range(len(pids))]
            wl = occ.addCurveLoop(lids); tools_s.append((2, occ.addPlaneSurface([wl])))
    tools_r = []
    if rotor_slots:
        poly = rotor_slot_polygon()
        for k in range(M.Nr):
            th = 2 * np.pi * (k + 0.5) / M.Nr + phi_rot
            P = place(poly, Rr, th, inward=True)
            pids = [occ.addPoint(x, y, 0) for (x, y) in P]
            lids = [occ.addLine(pids[i], pids[(i + 1) % len(pids)]) for i in range(len(pids))]
            wl = occ.addCurveLoop(lids); tools_r.append((2, occ.addPlaneSurface([wl])))
    # outer air = disk_s U stator slots ; rotor iron = disk_r \ rotor slots
    outer = [(2, disk_s)]
    if tools_s:
        outer, _ = occ.fuse([(2, disk_s)], tools_s, removeObject=True, removeTool=True)
    rotor = [(2, disk_r)]
    if tools_r:
        rotor, _ = occ.cut([(2, disk_r)], tools_r, removeObject=True, removeTool=True)
    air, _ = occ.cut(outer, rotor, removeObject=True, removeTool=True)
    occ.synchronize()
    # size field: fine in the gap, coarse away from it
    Rm = 0.5 * (Rs + Rr); g = Rs - Rr
    f = gmsh.model.mesh.field
    fid = f.add("MathEval")
    f.setString(fid, "F", f"{lc_gap} + {(lc_far - lc_gap) / 3.0}*Max(0, Abs(Sqrt(x*x+y*y) - {Rm}) - {0.5 * g})")
    f.setAsBackgroundMesh(fid)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)      # frontal-delaunay
    gmsh.model.mesh.generate(2)
    ntags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)[:, :2] * 1e-3           # to metres
    idx = {t: i for i, t in enumerate(ntags)}
    etypes, etags, enodes = gmsh.model.mesh.getElements(2)
    tri = np.vstack([np.array([idx[t] for t in en]).reshape(-1, 3) for et, en in zip(etypes, enodes) if et == 2])
    gmsh.finalize()
    used = np.unique(tri.ravel())
    remap = -np.ones(len(xyz), dtype=int); remap[used] = np.arange(len(used))
    xyz = xyz[used]; tri = remap[tri]
    return xyz, tri


def slot_body_ids(M, xyz, tri):
    """Stator-slot body index of each element (-1 outside the bodies), from
    the element centroid in the local slot frame."""
    c = xyz[tri].mean(axis=1)
    r = np.hypot(c[:, 0], c[:, 1]); ang = np.arctan2(c[:, 1], c[:, 0])
    y = r - M.Rs                           # depth into the stator (m)
    ids = np.full(len(tri), -1)
    inside = y > (HS0 + HS1) * 1e-3
    k = np.round((ang / (2 * np.pi)) * M.Ns - 0.5).astype(int) % M.Ns
    ids[inside] = k[inside]
    return ids


def assemble_P1(xyz, tri, coef):
    """Stiffness matrix of -div(coef grad u) with P1 triangles (coef per element)."""
    x = xyz[tri, 0]; y = xyz[tri, 1]
    b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1)
    c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1)
    area2 = (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])       # 2*area (signed)
    area = 0.5 * np.abs(area2)
    Ke = (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) / (4 * area[:, None, None]) * coef[:, None, None]
    I = np.repeat(tri, 3, axis=1).ravel(); J = np.tile(tri, (1, 3)).ravel()
    K = sp.coo_matrix((Ke.ravel(), (I, J)), shape=(len(xyz), len(xyz))).tocsr()
    return K, area


def solve_Az(M, xyz, tri, J_elem):
    """(1/mu0) K A = f with f_i = sum_e J_e area_e/3 ; gauge: A = 0 at node 0."""
    K, area = assemble_P1(xyz, tri, np.full(len(tri), 1.0 / MU0))
    f = np.zeros(len(xyz))
    np.add.at(f, tri.ravel(), np.repeat(J_elem * area / 3.0, 3))
    n = len(xyz); free = np.arange(1, n)
    Kf = K[free][:, free].tocsc()
    A = np.zeros(n)
    A[free] = spla.spsolve(Kf, f[free])
    return A, area


def midgap_fundamental_A(M, xyz, tri, A, nsamp=8192, r=None):
    """Fundamental (order p) of B_r on the circle r from the sampled A_z:
    B_r = (1/r) dA/dtheta  ->  |B_{r,p}| = (p/r) |A_p|  with A_p the complex
    Fourier amplitude 2/N sum A e^{-ip theta}.  Also returns the sampled A."""
    if r is None:
        r = M.Rm
    th = 2 * np.pi * np.arange(nsamp) / nsamp
    T = Triangulation(xyz[:, 0], xyz[:, 1], tri)
    interp = LinearTriInterpolator(T, A)
    a = interp(r * np.cos(th), r * np.sin(th))
    a = np.asarray(a.filled(np.nan))
    if np.isnan(a).any():
        raise RuntimeError("sampling circle leaves the mesh: %d NaN" % np.isnan(a).sum())
    Ap = (2.0 / nsamp) * np.sum(a * np.exp(-1j * M.p * th))
    coeffs = np.fft.rfft(a) * 2.0 / nsamp
    Bn = np.abs(coeffs) * np.arange(len(coeffs)) / r        # |B_{r,n}| for n>=0
    return (M.p / r) * abs(Ap), Bn, a, th


def smooth_Bg1_analytic(M, Us_staircase_fn, r=None):
    """Fundamental of B_r at radius r for the smooth annulus with the stator
    surface at the staircase potential F(theta) (jumps at the slot centres)
    and the rotor at 0.  F_p is the complex p-th Fourier amplitude of F."""
    if r is None:
        r = M.Rm
    # exact Fourier coefficient of a staircase with jumps dF_k at theta_k:
    #   F(theta) = sum_k dF_k H(theta - theta_k)  (mean irrelevant), p-th coeff:
    #   F_p = (1/(i p pi)) sum_k dF_k e^{-i p theta_k}  (for the series with cos/sin amplitude)
    Fp = Us_staircase_fn(M)
    v = np.log(r / M.Rr); X = M.X; n = M.p
    return MU0 * (n / r) * abs(Fp) * np.cosh(n * v) / np.sinh(n * X)


def staircase_Fp(M, I0=1.0):
    """p-th Fourier amplitude (cos/sin convention) of the MMF staircase built
    from the slot ampere-turns at the slot centres."""
    i3 = np.sqrt(2) * I0 * np.array([1.0, -0.5, -0.5])
    AT = M.slot_ampere_turns(i3)                 # jump at slot k, centred at 2pi(k+1/2)/Ns
    th = 2 * np.pi * (np.arange(M.Ns) + 0.5) / M.Ns
    # f(theta) = sum_k AT_k H(theta - th_k): coefficient a_n - i b_n = (1/pi) int f e^{-in theta}
    # int_{th_k}^{2pi} e^{-in theta} = (e^{-in th_k} - 1)/(i n)  -> the constant part drops for n>=1
    n = M.p
    c = (1 / np.pi) * np.sum(AT * (np.exp(-1j * n * th) - 1.0) / (1j * n))
    return c


def run_case(M, phi_rot, lc_gap=0.06, I0=1.0, stator_slots=True, rotor_slots=True, verbose=False):
    t0 = time.time()
    xyz, tri = build_air_mesh(M, phi_rot, lc_gap=lc_gap, stator_slots=stator_slots, rotor_slots=rotor_slots, verbose=verbose)
    ids = slot_body_ids(M, xyz, tri)
    i3 = np.sqrt(2) * I0 * np.array([1.0, -0.5, -0.5])
    AT = M.slot_ampere_turns(i3)
    K, area = assemble_P1(xyz, tri, np.ones(len(tri)))
    body_area = np.array([area[ids == k].sum() for k in range(M.Ns)])
    J = np.zeros(len(tri))
    m = ids >= 0
    J[m] = AT[ids[m]] / body_area[ids[m]]
    A, area = solve_Az(M, xyz, tri, J)
    bg1, Bn, a, th = midgap_fundamental_A(M, xyz, tri, A)
    Fp = staircase_Fp(M, I0)
    bg1_smooth = MU0 * (M.p / M.Rm) * abs(Fp) * np.cosh(M.p * np.log(M.Rm / M.Rr)) / np.sinh(M.p * M.X)
    return dict(kC=bg1_smooth / bg1, Bg1=bg1, Bg1_smooth=bg1_smooth, nnodes=len(xyz), ntri=len(tri),
                Bn=Bn, t=time.time() - t0, body_area=body_area, xyz=xyz, tri=tri, A=A)


if __name__ == '__main__':
    M = Machine()
    lc = float(sys.argv[1]) if len(sys.argv) > 1 else 0.08
    r = run_case(M, 0.0, lc_gap=lc, verbose=False)
    print(f"lc_gap={lc}: nodes {r['nnodes']}, tri {r['ntri']}, Bg1 slotted {r['Bg1']:.6f} T, smooth {r['Bg1_smooth']:.6f} T, kC = {r['kC']:.5f}   ({r['t']:.1f} s)")
    print("slot body areas (mm2):", r['body_area'][:3] * 1e6)
