# -*- coding: utf-8 -*-
"""
Post-processing of the numerical no-load test (T2).

For every exported run (NL_V*_NL_wave.tab: Time, Current, InputVoltage,
InducedVoltage, FluxLinkage of the three phases) and for the reference no-load
run (transitoire/a vide, Winding Plot 1..4), the fundamental phasors of v(t),
i(t), e(t) are formed over the window t in [1.0, 2.0) s (1000 samples of 1 ms,
exactly 50 periods), and the SAME construction as on the network side (RUN_Z7)
is applied:

    E_2D   = V - R_s I - j w L_ext I        (R_s = $Rs = 0.44574 ohm, L_ext = $Ls = 6.251120231328 mH)
    X_m,2D = |E_2D| / |I|

E_2D is the electromotive force of the two-dimensional flux linkage (it contains
the slot and tooth-tip leakage drop and the space-harmonic field), and the
exported InducedVoltage is its independent check.  Output: noload_results.json,
noload_results.txt, and the comparison with the network closures of
Z7_noload_net.mat at the same line voltages.
"""
import numpy as np, re, os, sys, json, glob
import scipy.io as sio
HERE = r"C:\Users\hp\Desktop\claude\T2_ansys"
EXP = os.path.join(HERE, "exports")
REF = r"C:\Users\hp\Desktop\ANSYS résultat 18.5KW\transitoire\a vide"
Z7 = r"C:\Users\hp\AppData\Local\Temp\claude\C--Users-hp-Desktop-claude\fb8f7acf-d703-40c5-b60c-b7b9dd7398fa\scratchpad\repo\doubly-slotted-annulus-IM\outputs\MEC_IM\Z7_noload_net.mat"
if not os.path.exists(Z7):
    raise SystemExit("Z7_noload_net.mat not found: " + Z7)
RS, LS, F = 0.44574, 6.251120231328e-3, 50.0
W = 2 * np.pi * F; XEXT = W * LS
T0, T1 = 1.0, 2.0


def read_tab(f):
    with open(f, encoding='utf-8', errors='replace') as fh:
        names = re.findall(r'"([^"]+)"', fh.readline())
    data = np.loadtxt(f, skiprows=1)
    return names, data


UNIT = {'A': 1.0, 'mA': 1e-3, 'kA': 1e3, 'uA': 1e-6, 'V': 1.0, 'mV': 1e-3, 'kV': 1e3, 'Wb': 1.0, 'mWb': 1e-3, 'uWb': 1e-6,
        'W': 1.0, 'mW': 1e-3, 'kW': 1e3, 'NewtonMeter': 1.0, 'mNewtonMeter': 1e-3, 'rpm': 1.0, 's': 1.0, 'ms': 1e-3}


def col(names, data, key):
    """Column by quantity name, converted to SI from the unit given in the header as 'name [unit]'."""
    for j, n in enumerate(names):
        if n.startswith(key):
            m = re.search(r'\[([^\]]+)\]', n)
            u = m.group(1).strip() if m else ''
            if u and u not in UNIT:
                raise ValueError("unknown unit %r in column %r" % (u, n))
            return data[:, j] * (UNIT[u] if u else 1.0)
    raise KeyError(key)


def phasor(t, y):
    sel = (t >= T0 - 1e-9) & (t < T1 - 1e-9)
    n = int(sel.sum())
    c = (2.0 / n) * np.sum(y[sel] * np.exp(-1j * W * t[sel]))       # peak complex amplitude, y = Re(c e^{jwt})
    return c / np.sqrt(2), n, np.sqrt(np.mean(y[sel] ** 2))           # rms phasor, sample count, true rms


def process(t, I3, V3, E3, P3=None, label=""):
    out = {'label': label, 'phases': {}}
    for ph in ['A', 'B', 'C']:
        I, n, Irms = phasor(t, I3[ph]); V, _, Vrms = phasor(t, V3[ph]); E, _, Erms = phasor(t, E3[ph])
        E2D = V - RS * I - 1j * XEXT * I
        Z = E2D / I
        d = {'n': n, 'I_rms_fund': abs(I), 'I_rms_true': Irms, 'V_rms_fund': abs(V), 'V_rms_true': Vrms,
             'E_exp_rms_fund': abs(E), 'E_exp_rms_true': Erms, 'E2D_rms': abs(E2D),
             'check_E2D_vs_export_pct': 100 * abs(E2D - E) / abs(E),
             'Xm2D': abs(E2D) / abs(I), 'Xm2D_reactive': Z.imag, 'R2D_inphase': Z.real,
             'angle_I_vs_V_deg': np.degrees(np.angle(I / V)), 'thd_I_pct': 100 * np.sqrt(max(Irms ** 2 - abs(I) ** 2, 0)) / abs(I)}
        if P3 is not None:
            P, _, _ = phasor(t, P3[ph])
            d['check_E_vs_jwPsi_pct'] = 100 * abs(E - 1j * W * P) / abs(E)
            d['check_E2D_vs_jwPsi_pct'] = 100 * abs(E2D - 1j * W * P) / abs(1j * W * P)
            d['jwPsi_rms'] = abs(W * P)
            d['lag_export_vs_jwPsi_deg'] = np.degrees(np.angle(E / (1j * W * P)))
        out['phases'][ph] = d
    keys = ['I_rms_fund', 'I_rms_true', 'V_rms_fund', 'E_exp_rms_fund', 'E2D_rms', 'check_E2D_vs_export_pct', 'Xm2D', 'Xm2D_reactive', 'R2D_inphase', 'angle_I_vs_V_deg', 'thd_I_pct']
    out['mean'] = {k: float(np.mean([out['phases'][p][k] for p in 'ABC'])) for k in keys}
    out['spread_I_pct'] = float(100 * (max(out['phases'][p]['I_rms_fund'] for p in 'ABC') - min(out['phases'][p]['I_rms_fund'] for p in 'ABC')) / out['mean']['I_rms_fund'])
    return out


results = []
# --- reference run (free-running rotor, 690 V)
names, data = read_tab(os.path.join(REF, "Winding Plot 4.tab")); t = data[:, 0]
I3 = {p: col(names, data, "Current(Phase_%s)" % p) for p in 'ABC'}
n3, d3 = read_tab(os.path.join(REF, "Winding Plot 3.tab")); V3 = {p: col(n3, d3, "InputVoltage(Phase_%s)" % p) for p in 'ABC'}
n2, d2 = read_tab(os.path.join(REF, "Winding Plot 2.tab")); E3 = {p: col(n2, d2, "InducedVoltage(Phase_%s)" % p) for p in 'ABC'}
n1, d1 = read_tab(os.path.join(REF, "Winding Plot 1.tab")); P3 = {p: col(n1, d1, "FluxLinkage(Phase_%s)" % p) for p in 'ABC'}
ns, ds = read_tab(os.path.join(REF, "la vitesse en fonction du temps.tab")); sp = ds[:, 1]; spd = float(np.mean(sp[(ds[:, 0] >= T0) & (ds[:, 0] < T1)]))
r = process(t, I3, V3, E3, P3, label="reference run 'a vide', 690 V, free rotor"); r['V_line'] = 690.0; r['speed_rpm'] = spd; r['kind'] = 'reference'
results.append(r)
# --- sweep runs
for f in sorted(glob.glob(os.path.join(EXP, "NL_V*_NL_wave.tab")), key=lambda s: -int(re.search(r'NL_V(\d+)', s).group(1))):
    V = int(re.search(r'NL_V(\d+)', f).group(1))
    names, data = read_tab(f); t = data[:, 0]
    I3 = {p: col(names, data, "Current(Phase_%s)" % p) for p in 'ABC'}
    V3 = {p: col(names, data, "InputVoltage(Phase_%s)" % p) for p in 'ABC'}
    E3 = {p: col(names, data, "InducedVoltage(Phase_%s)" % p) for p in 'ABC'}
    try:
        P3 = {p: col(names, data, "FluxLinkage(Phase_%s)" % p) for p in 'ABC'}
    except KeyError:
        P3 = None
    r = process(t, I3, V3, E3, P3, label="NL_V%d, 1500 rpm imposed" % V); r['V_line'] = float(V); r['kind'] = 'sweep'
    fm = f.replace("_NL_wave.tab", "_NL_misc.tab")
    if os.path.exists(fm):
        nm, dm = read_tab(fm); tm = dm[:, 0]; selm = (tm >= T0) & (tm < T1)
        for key in ['CoreLoss', 'Moving1.Torque', 'Moving1.Speed', 'SolidLoss']:
            try:
                r[key + '_mean'] = float(np.mean(col(nm, dm, key)[selm]))
            except KeyError:
                pass
    results.append(r)

lines = []
def P(s=""):
    lines.append(s); print(s)
P("=== T2 : numerical no-load test, FE side (window [%.1f, %.1f) s, R_s = %.5f ohm, X_ext = w L_ext = %.4f ohm) ===" % (T0, T1, RS, XEXT))
P("%-42s %7s %8s %9s %9s %9s %8s %8s %8s %8s %7s" % ('run', 'V_line', 'U_ph', 'I0 (A)', 'E_exp (V)', '|E_2D| (V)', 'chk %', 'X_m,2D', 'Im(Z)', 'Re(Z)', 'phi_I'))
for r in results:
    m = r['mean']
    P("%-42s %7.0f %8.2f %9.4f %9.2f %9.2f %8.3f %8.3f %8.3f %8.3f %7.2f" % (r['label'], r['V_line'], m['V_rms_fund'], m['I_rms_fund'], m['E_exp_rms_fund'], m['E2D_rms'], m['check_E2D_vs_export_pct'], m['Xm2D'], m['Xm2D_reactive'], m['R2D_inphase'], m['angle_I_vs_V_deg']))
    extra = " | phase spread of I0 %.2f %% | THD_I %.1f %%" % (r['spread_I_pct'], m['thd_I_pct'])
    if 'Moving1.Speed_mean' in r: extra += " | speed %.1f rpm | core loss %.1f W | torque %.3f N m" % (r['Moving1.Speed_mean'], r.get('CoreLoss_mean', np.nan), r.get('Moving1.Torque_mean', np.nan))
    if 'speed_rpm' in r: extra += " | speed %.1f rpm" % r['speed_rpm']
    if 'check_E_vs_jwPsi_pct' in r['phases']['A']: extra += " | export E vs jwPsi %.2f %% (lag %.2f deg) | E_2D vs jwPsi %.3f %% | |jwPsi| %.2f V" % (np.mean([r['phases'][p]['check_E_vs_jwPsi_pct'] for p in 'ABC']), np.mean([r['phases'][p]['lag_export_vs_jwPsi_deg'] for p in 'ABC']), np.mean([r['phases'][p]['check_E2D_vs_jwPsi_pct'] for p in 'ABC']), np.mean([r['phases'][p]['jwPsi_rms'] for p in 'ABC']))
    P(" " * 8 + extra)
P("  chk = |E_2D - InducedVoltage| / |InducedVoltage| on the fundamental phasors (the construction check)")

# --- network side (Z7)
if os.path.exists(Z7):
    z = sio.loadmat(Z7, squeeze_me=True)
    VL = np.atleast_1d(z['VL']).astype(float); OUT = z['OUT']; cases = [str(c) for c in np.atleast_1d(z['cases'][:, 0])] if z['cases'].ndim == 2 else [str(c) for c in z['cases']]
    XM0 = np.atleast_1d(z['XM0']); XEXTn = np.atleast_1d(z['XEXT'])
    P("\n=== comparison with the network closures (RUN_Z7 : E_2D = U_ph - R_s I - j X_ext I, X_ext(network) = %.4f ohm, R_s(network) = 0.4302 ohm) ===" % XEXTn[0])
    fe = {r['V_line']: r for r in results if r['kind'] == 'sweep'}
    ref690 = [r for r in results if r['kind'] == 'reference'][0]
    hdr = "%7s | %9s %9s | " % ('V_line', 'FE I0', 'FE X_m,2D') + " | ".join("%-28s" % c for c in cases)
    P(hdr)
    comp = {}
    for v in range(len(VL)):
        vl = VL[v]
        if vl in fe:
            m = fe[vl]['mean']; s = "%7.0f | %9.4f %9.3f | " % (vl, m['I_rms_fund'], m['Xm2D'])
            comp[str(int(vl))] = {'FE_I0': m['I_rms_fund'], 'FE_Xm2D': m['Xm2D']}
        elif vl == 690.0:
            m = ref690['mean']; s = "%7.0f*| %9.4f %9.3f | " % (vl, m['I_rms_fund'], m['Xm2D'])
            comp['690_ref'] = {'FE_I0': m['I_rms_fund'], 'FE_Xm2D': m['Xm2D']}
        else:
            s = "%7.0f | %9s %9s | " % (vl, '-', '-')
        cells = []
        for k in range(len(cases)):
            I0 = OUT[k, v, 1]; Xm2D = OUT[k, v, 7]
            if vl in fe:
                cells.append("I0 %7.4f (%+5.1f %%) X %7.3f (%+5.1f %%)" % (I0, 100 * (I0 - fe[vl]['mean']['I_rms_fund']) / fe[vl]['mean']['I_rms_fund'], Xm2D, 100 * (Xm2D - fe[vl]['mean']['Xm2D']) / fe[vl]['mean']['Xm2D']))
                comp[str(int(vl))][cases[k]] = {'I0': float(I0), 'Xm2D': float(Xm2D)}
            else:
                cells.append("I0 %7.4f            X %7.3f          " % (I0, Xm2D))
        P(s + " | ".join(cells))
    P("  (* reference run with the free rotor at 1499.9 rpm; the sweep run at 690 V has the rotor driven at 1500 rpm)")
    # --- low-flux extrapolation, same method both sides: linear in E_2D^2 through the two lowest voltages of the sweep
    P("\n=== low-flux extrapolation (linear in |E_2D|^2 through the two lowest sweep voltages, same method both sides) ===")
    low = sorted([vl for vl in fe.keys()])[:2]
    if len(low) == 2:
        def extrap(E2, X2):
            a = (X2[1] - X2[0]) / (E2[1] - E2[0]); return X2[0] - a * E2[0]
        Efe = [fe[v]['mean']['E2D_rms'] ** 2 for v in low]; Xfe = [fe[v]['mean']['Xm2D'] for v in low]
        xfe0 = extrap(Efe, Xfe); P("  FE : X_m,2D(%d V) = %.3f, X_m,2D(%d V) = %.3f -> X_m,2D(E -> 0) = %.3f ohm" % (low[0], Xfe[0], low[1], Xfe[1], xfe0))
        comp['extrap'] = {'FE': xfe0}
        for k in range(len(cases)):
            iv = [int(np.where(VL == v)[0][0]) for v in low]
            En = [OUT[k, i, 6] ** 2 for i in iv]; Xn = [OUT[k, i, 7] for i in iv]
            xn0 = extrap(En, Xn)
            P("  %-28s : X_m,2D(%d V) = %.3f, X_m,2D(%d V) = %.3f -> X_m,2D(E -> 0) = %.3f ohm (%+.1f %% vs FE) ; X_m0(0.2 A) of Table 8 = %.3f" % (cases[k], low[0], Xn[0], low[1], Xn[1], xn0, 100 * (xn0 - xfe0) / xfe0, XM0[k]))
            comp['extrap'][cases[k]] = float(xn0)
    P("  the network's own X_m(I_m) rises from 0.2 A to ~2 A (initial permeability of M800-50A) before falling: X_m0 of Table 8 is a low-current definition, not the unsaturated plateau")
    json.dump({'results': results, 'comparison': comp}, open(os.path.join(HERE, 'noload_results.json'), 'w'), indent=1, default=float)
open(os.path.join(HERE, 'noload_results.txt'), 'w', encoding='utf-8').write("\n".join(lines) + "\n")
