# -*- coding: utf-8 -*-
"""
Slot-cavity admittance element with ARBITRARY mouth-column centres (graded
tiling), and the rotor-cavity flux source of the bar current.

cavity_element_knots(kind, xc, L)
    Generalises cavity.cavity_element, whose interpolation knots are the
    centres of nO uniform columns: here the knots are passed explicitly
    (xc, in mm from the mouth centre) so that the cavity map matches a graded
    opening tiling.  With uniform xc it reproduces cavity.cavity_element.

cavity_source_rotor(xc, L)
    Flux source vector of the rotor cavity per ampere of bar current (T5 of
    the brief); see its docstring.
"""
import numpy as np, sys
sys.dont_write_bytecode = True
import scipy.sparse.linalg as spla
from fem_slots import assemble_P1, rotor_slot_polygon, BS0, HS0, BS1, HS1, BR0, HR0
from dsop import MU0
from cavity import _mesh_polygon


def _boundary_mask(poly, xyz):
    P = np.array(poly); nb = len(P)
    onb = np.zeros(len(xyz), bool)
    for i in range(nb):
        a = P[i]; b = P[(i + 1) % nb]; ab = b - a; l2 = ab @ ab
        t = np.clip(((xyz - a) @ ab) / l2, 0, 1)
        d = np.linalg.norm(xyz - (a + t[:, None] * ab), axis=1)
        onb |= d < 1e-6
    return onb


def cavity_element_knots(kind, xc, L, lc_mouth=0.02, lc_far=0.6):
    """Admittance matrix Q ((nO+2) x (nO+2)) of one slot cavity for mouth
    columns centred at xc (mm, sorted, inside (-b0/2, b0/2)).  Column order:
    mouth columns, U_left, U_right (rotor: single wall dof at index nO, index
    nO+1 zero).  Piecewise-linear mouth potential between the knots
    [-b0/2, xc..., b0/2] (corner values = adjacent tooth potentials)."""
    xc = np.asarray(xc, float); nO = len(xc)
    if kind == 'stator':
        b0 = BS0
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
    onb = _boundary_mask(poly, xyz)
    K = assemble_P1(xyz, tri, np.full(len(tri), MU0 * L))[0]
    n = len(xyz)
    mouth = np.where(on_mouth & onb)[0]
    xm = x[mouth]
    xk = np.concatenate([[-b0 / 2], xc, [b0 / 2]])
    ndof = nO + 2
    walls = onb & (~on_mouth)
    if kind == 'stator':
        on_bottom = np.abs(y - ybot) < tol
        left = walls & (x < 0) & (~on_bottom); right = walls & (x > 0) & (~on_bottom); bot = walls & on_bottom
        ramp = (x[bot] + BS1 / 2) / BS1

    def bc_vector(j):
        val = np.full(n, np.nan)
        kv = np.zeros(nO + 2)
        if j < nO:
            kv[j + 1] = 1.0
        elif j == nO:
            kv[0] = 1.0
        else:
            kv[-1] = 1.0
        val[mouth] = np.interp(xm, xk, kv)
        if kind == 'stator':
            val[left] = 1.0 if j == nO else 0.0
            val[right] = 1.0 if j == nO + 1 else 0.0
            if j == nO:
                val[bot] = 1.0 - ramp
            elif j == nO + 1:
                val[bot] = ramp
            else:
                val[bot] = 0.0
        else:
            val[walls] = 1.0 if j == nO else 0.0
        return val
    fixed = np.where(onb)[0]; free = np.where(~onb)[0]
    Kff = K[free][:, free].tocsc(); Kfd = K[free][:, fixed]
    lu = spla.splu(Kff)
    sols = []
    for j in range(ndof):
        val = bc_vector(j)
        u = np.zeros(n); u[fixed] = val[fixed]
        u[free] = lu.solve(-Kfd @ u[fixed])
        sols.append(u)
    Q = np.zeros((ndof, ndof))
    for a in range(ndof):
        for b in range(ndof):
            Q[a, b] = sols[a] @ (K @ sols[b])
    Q = 0.5 * (Q + Q.T)
    if kind == 'rotor':
        Qm = np.zeros((nO + 2, nO + 2))
        Qm[:nO + 1, :nO + 1] = Q[:nO + 1, :nO + 1]
        Qm[:nO + 1, nO] += Q[:nO + 1, nO + 1]; Qm[nO, :nO + 1] += Q[nO + 1, :nO + 1]; Qm[nO, nO] += Q[nO + 1, nO + 1]
        Q = Qm
    return Q, dict(nnodes=n, xyz=xyz, tri=tri, mouth=mouth)


def cavity_source_rotor(xc, L, lc_mouth=0.02, lc_far=0.6):
    """Flux source vector s (nO+1,) of the rotor cavity per ampere of bar
    current, on the nO mouth columns and the wall dof.

    Decomposition (reduced formulation).  With a bar current I (along +z) the
    tooth on the right of the slot sits I BELOW the tooth on the left: on a
    contour that runs along the mouth in +x and closes through infinitely
    permeable iron, Ampere gives int_mouth H_x dx = +I.  The cavity field is
    the sum of
      (a) the Laplace field of cavity_element_knots('rotor', ...), driven by
          the mouth-column potentials measured from the (single) wall
          potential, and
      (b) a particular field A_J of the uniform current density J = I/A_bar in
          the bar region, with H_t = 0 on the iron walls and H_x = I/b0
          uniform across the mouth, i.e. the mouth potential is the linear
          ramp U_wall -> U_wall - I from the left corner to the right corner:
          this field carries the potential jump between the corners.
    For network column potentials U_k the flux delivered to the column and
    wall dofs is   Phi = Q [U; U_wall] + s I   with
          s = Phi_J - Q [ramp; 0],   ramp_k = -(xc_k + b0/2)/b0 ,
    Phi_J being the hat-weighted flux of (b) entering the cavity through the
    mouth columns and, on the wall dof, minus the total mouth flux plus the
    corner shares.  In a rectangular slot A_J depends on y only and Phi_J = 0;
    in the pear-shaped slot it is small but nonzero.  The approximation is the
    one made everywhere in the cavity model: the iron walls are infinitely
    permeable, so the reluctance drop of the rotor iron around the slot is not
    carried by the element (the network carries it).

    Units: coordinates in mm; the P1 stiffness of the Laplacian is invariant
    to the length unit, so with coef = 1/mu0, J in A/mm^2 and areas in mm^2
    the load is in A and A_z comes out in Wb/m.  Local frame: x = counter-
    clockwise tangent, y into the rotor; (x, y, z) is right-handed with z the
    machine axis.  B = curl(A z) = (dA/dy, -dA/dx), H_x = (1/mu0) dA/dy, and the
    outward normal on the mouth is -y, so (1/mu0) dA/dn_out = -H_x = -I/b0.
    Weak form: K A = f_J - f_b with f_b,i = (I/b0) int_mouth phi_i ds, and the
    pure Neumann problem is compatible because sum f_J = I = sum f_b."""
    xc = np.asarray(xc, float); nO = len(xc); b0 = BR0
    poly = rotor_slot_polygon(depth_extra=0.0)
    xyz, tri = _mesh_polygon(poly, lc_mouth, lc_far)
    x, y = xyz[:, 0], xyz[:, 1]
    onb = _boundary_mask(poly, xyz)
    on_mouth = (np.abs(y) < 1e-7) & onb
    K1, area = assemble_P1(xyz, tri, np.full(len(tri), 1.0 / MU0))
    cen = xyz[tri].mean(axis=1)
    inbar = cen[:, 1] > HR0 - 1e-9                                     # bar region = below the isthmus
    Abar = area[inbar].sum()                                            # mm^2
    I = 1.0
    J = np.zeros(len(tri)); J[inbar] = I / Abar                         # A/mm^2
    f = np.zeros(len(xyz))
    np.add.at(f, tri.ravel(), np.repeat(J * area / 3.0, 3))            # A
    mouth = np.where(on_mouth)[0]
    xm = x[mouth]; o = np.argsort(xm); mouth = mouth[o]; xm = xm[o]
    seg = np.diff(xm)                                                   # mm
    g = I / b0                                                          # A/mm
    fb = np.zeros(len(xyz))
    fb[mouth[:-1]] += g * seg / 2; fb[mouth[1:]] += g * seg / 2
    rhs = f - fb
    tot_src = f.sum(); tot_bnd = fb.sum(); sumrhs = rhs.sum()
    if abs(sumrhs) > 1e-9 * abs(tot_src):
        raise RuntimeError("Neumann problem not compatible: sum f_J = %.6e, sum f_b = %.6e" % (tot_src, tot_bnd))
    rhs = rhs - rhs.mean()
    n = len(xyz); fix = mouth[0]; free = np.array([i for i in range(n) if i != fix])
    Kff = K1[free][:, free].tocsc()
    A = np.zeros(n); A[free] = spla.spsolve(Kff, rhs[free])            # Wb/m, gauge A = 0 at the left corner
    # flux entering the cavity through each hat-weighted mouth column:
    #   Phi_k = L int chi_k B_y dx = L int chi_k (-dA/dx) dx = -L sum_seg chi_mid dA   (the mm cancel)
    xk = np.concatenate([[-b0 / 2], xc, [b0 / 2]])
    Am = A[mouth]; dA = np.diff(Am)
    Phi = np.zeros(nO + 1)
    for k in range(nO):
        kv = np.zeros(nO + 2); kv[k + 1] = 1.0
        chi = np.interp(xm, xk, kv)
        Phi[k] = -L * np.sum(0.5 * (chi[:-1] + chi[1:]) * dA)         # Wb per ampere
    kvL = np.zeros(nO + 2); kvL[0] = 1.0; chiL = np.interp(xm, xk, kvL)
    kvR = np.zeros(nO + 2); kvR[-1] = 1.0; chiR = np.interp(xm, xk, kvR)
    Phi_corners = -L * np.sum(0.5 * (chiL[:-1] + chiL[1:]) * dA) - L * np.sum(0.5 * (chiR[:-1] + chiR[1:]) * dA)
    Phi_mouth_total = -L * np.sum(dA)                                   # = -L (A_right - A_left)
    # wall dof: corner shares of the mouth flux, minus everything that entered through the mouth
    Phi[nO] = Phi_corners - Phi_mouth_total
    # tangential field check on the mouth: (1/mu0) dA/dy = H_x should be I/b0 = 0.5 A/mm
    return Phi, dict(Abar_mm2=Abar, tot_src=tot_src, tot_bnd=tot_bnd, sumrhs=sumrhs,
                     nnodes=n, A=A, xyz=xyz, tri=tri, mouth=mouth, Phi_mouth_total=Phi_mouth_total,
                     Phi_corners=Phi_corners)


def rotor_source_vector(nO, L, q=1.0, lc_mouth=None):
    """s = Phi_J - Q_r [ramp; 0] for uniform (q = 1) or graded mouth columns.
    Returns s (nO+1,), Q_r, Phi_J, ramp."""
    b0 = BR0
    if q == 1.0:
        xc = -b0 / 2 + (np.arange(nO) + 0.5) * b0 / nO
    else:
        from t4_graded import opening_centres_mm
        xc = opening_centres_mm(nO, q)
    if lc_mouth is None:
        lc_mouth = 0.02 if nO <= 8 else 0.01
    Qr, _ = cavity_element_knots('rotor', xc, L, lc_mouth=lc_mouth)
    PhiJ, info = cavity_source_rotor(xc, L, lc_mouth=lc_mouth)
    ramp = np.zeros(nO + 1); ramp[:nO] = -(xc + b0 / 2) / b0
    s = PhiJ - Qr[:nO + 1, :nO + 1] @ ramp
    return s, Qr, PhiJ, ramp, xc, info
