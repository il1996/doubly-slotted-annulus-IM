# -*- coding: utf-8 -*-
"""
Independent re-implementation (Python/NumPy) of the doubly slotted annulus
operator of the manuscript, Section 2, in the linear (infinitely permeable
iron) limit.

Nothing here is taken from the MATLAB chain: the operator is rebuilt from
equations (1)-(5) of the manuscript and from the geometry of Table 5.

Conventions
-----------
* annulus  R_r <= r <= R_s ,  X_g = ln(R_s/R_r),  v = ln(r/R_r)
* bore quantity  f(theta) = f_0 + sum_{n>=1} f_n^c cos(n theta) + f_n^s sin(n theta)
* magnetic energy (per harmonic, both quadratures)
      W_n = mu0 pi n L / 2 [ coth(n X)(|u_s|^2+|u_r|^2) - 2 <u_s,u_r>/sinh(n X) ]
      W_0 = pi mu0 L / X  (u_s0 - u_r0)^2         (homopolar, Lambda_0 = 2 pi mu0 L / X)
* surface potential  u(theta) = sum_k U_k chi_k(theta)  with basis chi_k
      C[n,k] = (1/pi) int chi_k cos(n theta),  S[n,k] = (1/pi) int chi_k sin(n theta)
      w[k]   = (1/2pi) int chi_k
* admittance operator  A = d^2 W / dU dU  (symmetric, positive semi-definite,
  null vector [1;1]).  Phi = A U is the flux each node must supply.
  (The manuscript writes Y = -A; the sign is immaterial to every ratio.)

Bases
-----
'p0'   : indicator of the column (piecewise constant)
'p1'   : symmetric hat of half-width equal to the column's own width
         (the manuscript's hat basis; partition of unity only on a uniform grid)
'p1a'  : asymmetric hat whose two half-widths are the distances to the
         neighbouring nodes (partition of unity on any grid; the remedy the
         manuscript proposes in Section 3.1)
"""
import numpy as np

MU0 = 4e-7 * np.pi


# ----------------------------------------------------------------- machine
class Machine:
    """18.5 kW, 48/44 cage induction machine, Table 5 of the manuscript."""
    def __init__(self):
        self.Ds = 163.777468232983e-3      # stator bore diameter, from the FE project
        self.Dr = 163.2915e-3              # rotor outer diameter (Table 5)
        self.Rs = self.Ds / 2
        self.Rr = self.Dr / 2
        self.g = self.Rs - self.Rr
        self.L = 164.782448048495e-3
        self.Ns, self.Nr, self.p = 48, 44, 2
        self.bs0, self.br0 = 2.0e-3, 2.0e-3   # slot openings
        self.Nph, self.kw1 = 144, 0.9250
        self.Ntc, self.a_par, self.layers, self.pitch = 18, 2, 2, 10
        self.f = 50.0
        self.w = 2 * np.pi * self.f
        self.X = np.log(self.Rs / self.Rr)
        self.Rm = 0.5 * (self.Rs + self.Rr)

    # winding MMF at the teeth --------------------------------------------
    def slot_ampere_turns(self, i3):
        """Ampere-turns of each of the 48 stator slots for phase currents
        i3 = (iA, iB, iC).  Double layer, q = 4, coil pitch 10 (of 12),
        18 conductors per coil side, 2 parallel paths."""
        Ns, q = self.Ns, 4
        # top layer: phase belt sequence A+ C- B+ A- C+ B- (60 deg belts)
        belt = [(0, +1), (2, -1), (1, +1), (0, -1), (2, +1), (1, -1)]
        top = np.zeros(Ns)
        for j in range(Ns):
            ph, sg = belt[(j // q) % 6]
            top[j] = sg * i3[ph]
        # bottom layer of slot j+pitch is the return side of the coil in top layer j
        bot = np.zeros(Ns)
        for j in range(Ns):
            bot[(j + self.pitch) % Ns] = -top[j]
        per_cond = self.Ntc / self.a_par          # conductors x (current per path)
        return per_cond * (top + bot)

    def tooth_mmf(self, i3):
        """Scalar potential of the 48 stator teeth (A) for currents i3.
        Tooth k lies between slot k-1 and slot k (slot k centred at
        theta = 2 pi (k + 1/2)/Ns, tooth k at 2 pi k/Ns).  Mean removed."""
        AT = self.slot_ampere_turns(i3)
        F = np.cumsum(AT)                # potential of tooth k+1 = F[k]
        F = np.roll(F, 1)                # tooth k
        return F - F.mean()

    def fundamental_check(self, I0=1.0):
        """Return the fundamental (order p) of the tooth MMF and the
        theoretical value (3/2)(4/pi)(Nph kw1/(2p)) sqrt2 I0."""
        i3 = np.sqrt(2) * I0 * np.array([1.0, -0.5, -0.5])
        F = self.tooth_mmf(i3)
        th = 2 * np.pi * np.arange(self.Ns) / self.Ns
        c = (2 / self.Ns) * np.sum(F * np.exp(-1j * self.p * th))
        theo = 1.5 * (4 / np.pi) * (self.Nph * self.kw1 / (2 * self.p)) * np.sqrt(2) * I0
        return abs(c), theo


# ------------------------------------------------------------------ tiling
def tile_surface(N, R, b0, nT, nO, phi=0.0):
    """Column centres, widths and neighbour half-widths of one slotted bore.

    Returns dict with
      th     : column centre angles (rad), ordered along the bore
      d      : column widths (rad)
      isface : True for tooth-face columns
      tooth  : tooth index of each face column (-1 for openings)
      hl, hr : distances to the previous / next column centre (for 'p1a')
    Tooth k is centred at 2 pi k / N + phi; the opening after tooth k is
    centred at 2 pi (k + 1/2)/N + phi.
    """
    tau = 2 * np.pi / N
    a_o = b0 / R                    # opening angle
    a_f = tau - a_o                 # face angle
    th, d, isface, tooth = [], [], [], []
    for k in range(N):
        c = tau * k + phi
        # face columns
        for i in range(nT):
            th.append(c - a_f / 2 + (i + 0.5) * a_f / nT); d.append(a_f / nT)
            isface.append(True); tooth.append(k)
        # opening columns
        for i in range(nO):
            th.append(c + a_f / 2 + (i + 0.5) * a_o / nO); d.append(a_o / nO)
            isface.append(False); tooth.append(-1)
    th = np.array(th); d = np.array(d); isface = np.array(isface); tooth = np.array(tooth)
    # neighbour distances (periodic)
    thp = np.roll(th, 1); thn = np.roll(th, -1)
    hl = (th - thp) % (2 * np.pi); hr = (thn - th) % (2 * np.pi)
    return dict(th=th, d=d, isface=isface, tooth=tooth, hl=hl, hr=hr, N=N, R=R)


def smooth_surface(N, phi=0.0):
    """One column per tooth pitch, paving the bore exactly (comparator (a))."""
    tau = 2 * np.pi / N
    th = tau * np.arange(N) + phi
    d = np.full(N, tau)
    return dict(th=th, d=d, isface=np.ones(N, bool), tooth=np.arange(N),
                hl=np.full(N, tau), hr=np.full(N, tau), N=N, R=None)


# ---------------------------------------------------------- projections
def projections(surf, Nh, basis):
    """C[n,k], S[n,k] for n = 1..Nh and the homopolar weights w[k]."""
    th, d, hl, hr = surf['th'], surf['d'], surf['hl'], surf['hr']
    n = np.arange(1, Nh + 1, dtype=float)[:, None]
    if basis == 'p0':
        amp = (2.0 / (n * np.pi)) * np.sin(n * d[None, :] / 2)          # real, even
        T = amp * np.exp(-1j * n * th[None, :])
        w = d / (2 * np.pi)
    elif basis == 'p1':
        h = d[None, :]
        amp = (4.0 / (np.pi * n ** 2 * h)) * np.sin(n * h / 2) ** 2
        T = amp * np.exp(-1j * n * th[None, :])
        w = d / (2 * np.pi)                     # integral of a hat of half-width d is d
    elif basis == 'p1a':
        HL = hl[None, :]; HR = hr[None, :]
        Tk = (1 - np.exp(-1j * n * HR)) / (n ** 2 * HR) + (1 - np.exp(1j * n * HL)) / (n ** 2 * HL)
        T = Tk * np.exp(-1j * n * th[None, :]) / np.pi
        w = 0.5 * (hl + hr) / (2 * np.pi)
    else:
        raise ValueError(basis)
    C = np.ascontiguousarray(T.real)
    S = np.ascontiguousarray(-T.imag)
    return C, S, w


# ------------------------------------------------------------- operator
def assemble(M, surf_s, surf_r, Nh, basis):
    """Fine-grid admittance operator A (positive semi-definite) of the
    annulus with the two tiled surfaces.  Returns A and the projections."""
    Cs, Ss, ws = projections(surf_s, Nh, basis)
    Cr, Sr, wr = projections(surf_r, Nh, basis)
    n = np.arange(1, Nh + 1, dtype=float)
    x = n * M.X
    coth = 1.0 / np.tanh(x)
    csch = 1.0 / np.sinh(x)
    Kd = MU0 * np.pi * M.L * n * coth          # diagonal kernel
    Ko = MU0 * np.pi * M.L * n * csch          # coupling kernel
    Ass = (Cs.T * Kd) @ Cs + (Ss.T * Kd) @ Ss
    Arr = (Cr.T * Kd) @ Cr + (Sr.T * Kd) @ Sr
    Asr = -((Cs.T * Ko) @ Cr + (Ss.T * Ko) @ Sr)
    Lam0 = 2 * np.pi * MU0 * M.L / M.X
    ms, mr = len(ws), len(wr)
    A = np.zeros((ms + mr, ms + mr))
    A[:ms, :ms] = Ass + Lam0 * np.outer(ws, ws)
    A[ms:, ms:] = Arr + Lam0 * np.outer(wr, wr)
    A[:ms, ms:] = Asr - Lam0 * np.outer(ws, wr)
    A[ms:, :ms] = A[:ms, ms:].T
    proj = dict(Cs=Cs, Ss=Ss, Cr=Cr, Sr=Sr, ws=ws, wr=wr, Kd=Kd, Ko=Ko, n=n)
    return A, proj


def condense(A, surf_s, surf_r, opening='neumann', slot_adm=None):
    """Condense the fine operator onto the tooth nodes.
    opening='neumann' : Phi_O = 0 on the opening columns (the manuscript, eq. 5)
    slot_adm          : optional matrix added to the opening block before the
                        Schur complement (cavity admittance of the slots)."""
    ms = len(surf_s['th']); mr = len(surf_r['th'])
    isf = np.concatenate([surf_s['isface'], surf_r['isface']])
    tooth = np.concatenate([surf_s['tooth'], np.where(surf_r['tooth'] >= 0, surf_r['tooth'] + surf_s['N'], -1)])
    T = np.where(isf)[0]; O = np.where(~isf)[0]
    ATT = A[np.ix_(T, T)]; ATO = A[np.ix_(T, O)]; AOO = A[np.ix_(O, O)].copy()
    if slot_adm is not None:
        AOO += slot_adm
    if len(O):
        X = np.linalg.solve(AOO, ATO.T)
        Sc = ATT - ATO @ X
    else:
        Sc = ATT; X = None
    nt = surf_s['N'] + surf_r['N']
    P = np.zeros((len(T), nt))
    P[np.arange(len(T)), tooth[T]] = 1.0
    Yc = P.T @ Sc @ P
    Yc = 0.5 * (Yc + Yc.T)
    return dict(Yc=Yc, T=T, O=O, P=P, Xo=X, tooth=tooth, ms=ms, mr=mr)


def solve_linear(M, cond, Us):
    """Infinite iron: stator teeth at potentials Us (48), rotor teeth at a
    common floating potential fixed by zero net rotor flux.
    Returns the tooth potentials U (92) and the tooth fluxes Phi (92)."""
    Yc = cond['Yc']; Ns = M.Ns
    Yss = Yc[:Ns, :Ns]; Ysr = Yc[:Ns, Ns:]; Yrr = Yc[Ns:, Ns:]
    one = np.ones(M.Nr)
    Ur = -(one @ (Ysr.T @ Us)) / (one @ Yrr @ one)
    U = np.concatenate([Us, Ur * one])
    Phi = Yc @ U
    return U, Phi


def fine_potentials(cond, U):
    """Recover the fine-grid potentials from the tooth potentials."""
    T, O, P, X = cond['T'], cond['O'], cond['P'], cond['Xo']
    nfine = cond['ms'] + cond['mr']
    u = np.zeros(nfine)
    UT = P @ U
    u[T] = UT
    if len(O):
        u[O] = -X @ UT
    return u


def gap_field_harmonics(M, proj, u, cond, r=None):
    """Fourier coefficients (cos, sin) of B_r at radius r (default mid-gap)
    from the fine-grid surface potentials u.  B_r = -mu0 dphi/dr."""
    if r is None:
        r = M.Rm
    ms = cond['ms']
    us, ur = u[:ms], u[ms:]
    usc = proj['Cs'] @ us; uss = proj['Ss'] @ us
    urc = proj['Cr'] @ ur; urs = proj['Sr'] @ ur
    n = proj['n']; v = np.log(r / M.Rr); X = M.X
    # dphi/dv = n [ u_s cosh(n v) - u_r cosh(n (X - v)) ] / sinh(n X)
    a = n * np.cosh(n * v) / np.sinh(n * X)
    b = n * np.cosh(n * (X - v)) / np.sinh(n * X)
    Brc = -(MU0 / r) * (a * usc - b * urc)
    Brs = -(MU0 / r) * (a * uss - b * urs)
    return Brc, Brs


def Bg1(M, proj, u, cond, r=None):
    Brc, Brs = gap_field_harmonics(M, proj, u, cond, r)
    k = M.p - 1
    return np.hypot(Brc[k], Brs[k])


def Xm_from_Bg1(M, bg1, I0):
    """X_m = w kw1 Nph ((2/pi) Bg1 tau_p L) / (sqrt2 I0), tau_p = pi Ds/(2p)."""
    tp = np.pi * M.Ds / (2 * M.p)
    return M.w * M.kw1 * M.Nph * ((2 / np.pi) * bg1 * tp * M.L) / np.sqrt(2) / I0


# ------------------------------------------------------------ invariants
def invariants(A, proj, M):
    """I-1 (uniform potential drives no flux) and I-4 (homopolar permeance)."""
    nfine = A.shape[0]
    one = np.ones(nfine)
    I1 = np.max(np.abs(A @ one)) / np.max(np.abs(A))
    ms = len(proj['ws'])
    # raise the stator uniformly by 1, rotor at 0 -> total flux = Lambda_0
    u = np.zeros(nfine); u[:ms] = 1.0
    Phi = A @ u
    Lam_meas = np.sum(Phi[:ms])
    Lam0 = 2 * np.pi * MU0 * M.L / M.X
    I4 = abs(Lam_meas - Lam0) / Lam0
    sym = np.linalg.norm(A - A.T) / np.linalg.norm(A)
    return dict(I1=I1, I4=I4, Lam_meas=Lam_meas, Lam0=Lam0, sym=sym)


# ----------------------------------------------------------- one k_C run
def slotting_ratio(M, nT, nO, Nh, basis='p1', phi=0.0, I0=1.0, opening='neumann',
                   slot_adm=None, return_all=False):
    """k_C = X_m^smooth / X_m^slotted at fixed MMF, infinite iron.
    The smooth comparator is definition (a) of the manuscript: one node per
    tooth, arcs paving the bore, same basis, same Nh."""
    i3 = np.sqrt(2) * I0 * np.array([1.0, -0.5, -0.5])
    Us = M.tooth_mmf(i3)
    # slotted
    ss = tile_surface(M.Ns, M.Rs, M.bs0, nT, nO, 0.0)
    sr = tile_surface(M.Nr, M.Rr, M.br0, nT, nO, phi)
    A, proj = assemble(M, ss, sr, Nh, basis)
    cond = condense(A, ss, sr, opening, slot_adm)
    U, Phi = solve_linear(M, cond, Us)
    u = fine_potentials(cond, U)
    b_slot = Bg1(M, proj, u, cond)
    # smooth (a)
    ss0 = smooth_surface(M.Ns, 0.0); sr0 = smooth_surface(M.Nr, phi)
    A0, proj0 = assemble(M, ss0, sr0, Nh, basis)
    cond0 = condense(A0, ss0, sr0)
    U0, Phi0 = solve_linear(M, cond0, Us)
    u0 = fine_potentials(cond0, U0)
    b_smooth = Bg1(M, proj0, u0, cond0)
    out = dict(kC=b_smooth / b_slot, Bg1_slot=b_slot, Bg1_smooth=b_smooth,
               Xm_slot=Xm_from_Bg1(M, b_slot, I0), Xm_smooth=Xm_from_Bg1(M, b_smooth, I0),
               ncol=A.shape[0])
    if return_all:
        out.update(A=A, proj=proj, cond=cond, U=U, Phi=Phi, u=u, ss=ss, sr=sr)
    return out


def carter(M):
    """Carter's factor, approximate (8) and exact conformal-mapping forms."""
    def surf(b0, g, tau):
        x = b0 / g
        gam = x ** 2 / (5 + x)
        kap = tau / (tau - gam * g)
        xx = b0 / (2 * g)
        gam_ex = (4 / np.pi) * (xx * np.arctan(xx) - np.log(np.sqrt(1 + xx ** 2)))
        kex = tau / (tau - gam_ex * g)
        return kap, kex
    tau_s = 2 * np.pi * M.Rs / M.Ns; tau_r = 2 * np.pi * M.Rr / M.Nr
    ks, kse = surf(M.bs0, M.g, tau_s)
    kr, kre = surf(M.br0, M.g, tau_r)
    return dict(ks=ks, kr=kr, k=ks * kr, ks_exact=kse, kr_exact=kre, k_exact=kse * kre,
                tau_s=tau_s, tau_r=tau_r)


if __name__ == '__main__':
    M = Machine()
    print(f"Rs={M.Rs*1e3:.5f} mm Rr={M.Rr*1e3:.5f} mm g={M.g*1e3:.5f} mm X={M.X:.6e} L={M.L*1e3:.4f} mm")
    f1, theo = M.fundamental_check()
    print(f"MMF fundamental: staircase {f1:.4f} A  vs theory {theo:.4f} A  (ratio {f1/theo:.5f})")
    print("Carter:", carter(M))


# ------------------------------------------------ cavity-augmented condensation
def condense_cavity(A, surf_s, surf_r, Qs, Qr, nO):
    """Condensation onto the tooth nodes with the slot-cavity admittances added
    on the opening columns (see cavity.py).  Qs: (nO+2)^2 stator cavity matrix,
    dofs = [mouth columns, U_left tooth, U_right tooth]; Qr: (nO+2)^2 rotor
    cavity matrix with the wall dof at index nO (index nO+1 unused/zero)."""
    Ns, Nr = surf_s['N'], surf_r['N']
    ms = len(surf_s['th']); nt = Ns + Nr
    isf = np.concatenate([surf_s['isface'], surf_r['isface']])
    tooth = np.concatenate([surf_s['tooth'], np.where(surf_r['tooth'] >= 0, surf_r['tooth'] + Ns, -1)])
    T = np.where(isf)[0]; O = np.where(~isf)[0]
    P = np.zeros((len(T), nt)); P[np.arange(len(T)), tooth[T]] = 1.0
    ATT = A[np.ix_(T, T)]; ATO = A[np.ix_(T, O)]; AOO = A[np.ix_(O, O)].copy()
    nO_tot = len(O)
    # map: global opening column index -> position in O
    pos = -np.ones(A.shape[0], int); pos[O] = np.arange(nO_tot)
    Qtt = np.zeros((nt, nt)); QtO = np.zeros((nt, nO_tot)); QOO = np.zeros((nO_tot, nO_tot))
    # stator slots: opening after tooth k lies between tooth k and k+1; its columns follow the face columns of tooth k
    def add_slot(cols, left, right, Q):
        idx_O = pos[cols]
        d = list(idx_O)
        QOO[np.ix_(idx_O, idx_O)] += Q[:nO, :nO]
        for a, t in enumerate([left, right]):
            if t is None:
                continue
            QtO[t, idx_O] += Q[nO + a, :nO]
            for b, t2 in enumerate([left, right]):
                if t2 is None:
                    continue
                Qtt[t, t2] += Q[nO + a, nO + b]
    # stator
    th_s = surf_s['th']
    for k in range(Ns):
        cols = np.where((~surf_s['isface']) & (np.arange(ms) // (len(surf_s['th']) // Ns) == k))[0]
        add_slot(cols, k, (k + 1) % Ns, Qs)
    # rotor
    mr = len(surf_r['th'])
    for k in range(Nr):
        cols_local = np.where((~surf_r['isface']) & (np.arange(mr) // (mr // Nr) == k))[0]
        cols = cols_local + ms
        add_slot(cols, Ns + k, None, Qr)
    AOOq = AOO + QOO
    Bmat = ATO.T @ P + QtO.T            # (nO_tot x nt)
    X = np.linalg.solve(AOOq, Bmat)
    Yc = P.T @ ATT @ P + Qtt - Bmat.T @ X
    Yc = 0.5 * (Yc + Yc.T)
    return dict(Yc=Yc, T=T, O=O, P=P, Xo=X, tooth=tooth, ms=ms, mr=mr, cavity=True)


def fine_potentials_cavity(cond, U):
    T, O, P, X = cond['T'], cond['O'], cond['P'], cond['Xo']
    u = np.zeros(cond['ms'] + cond['mr'])
    u[T] = P @ U
    u[O] = -X @ U
    return u


def slotting_ratio_cavity(M, nT, nO, Nh, Qs, Qr, basis='p1', phi=0.0, I0=1.0, return_all=False):
    i3 = np.sqrt(2) * I0 * np.array([1.0, -0.5, -0.5]); Us = M.tooth_mmf(i3)
    ss = tile_surface(M.Ns, M.Rs, M.bs0, nT, nO, 0.0); sr = tile_surface(M.Nr, M.Rr, M.br0, nT, nO, phi)
    A, proj = assemble(M, ss, sr, Nh, basis)
    cond = condense_cavity(A, ss, sr, Qs, Qr, nO)
    U, Phi = solve_linear(M, cond, Us)
    u = fine_potentials_cavity(cond, U)
    b_slot = Bg1(M, proj, u, cond)
    ss0 = smooth_surface(M.Ns, 0.0); sr0 = smooth_surface(M.Nr, phi)
    A0, proj0 = assemble(M, ss0, sr0, Nh, basis); cond0 = condense(A0, ss0, sr0)
    U0, _ = solve_linear(M, cond0, Us); u0 = fine_potentials(cond0, U0)
    b_smooth = Bg1(M, proj0, u0, cond0)
    out = dict(kC=b_smooth / b_slot, Bg1_slot=b_slot, Bg1_smooth=b_smooth, ncol=A.shape[0])
    if return_all:
        out.update(A=A, proj=proj, cond=cond, U=U, Phi=Phi, u=u)
    return out
