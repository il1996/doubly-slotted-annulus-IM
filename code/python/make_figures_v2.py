# -*- coding: utf-8 -*-
"""Figures of the manuscript, all from computed data.  Version 2 (16 Sept. 2026):
paths adapted to the completed archive; Fig. 3(b) adds the graded tiling of
T4; Fig. 5 and Fig. 6 add the cavity-coupled closure (Z5_sweep_cav.mat,
Z4_fields_cav_*.txt) next to the Phi_O = 0 curves."""
import json, numpy as np, os, sys, scipy.io as sio
sys.dont_write_bytecode = True
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Wedge
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.6,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False,
                     'axes.titlesize': 9, 'axes.labelsize': 8.5, 'legend.fontsize': 7.5})
ARCH = r'<home>\AppData\Local\Temp\claude\C--Users-hp-Desktop-claude\fb8f7acf-d703-40c5-b60c-b7b9dd7398fa\scratchpad\repo\doubly-slotted-annulus-IM'
OUT = os.path.join(ARCH, 'code', 'article', 'figures')
W = os.path.join(ARCH, 'outputs', 'python')
MECO = os.path.join(ARCH, 'outputs', 'MEC_IM')        # MATLAB transcripts and .mat files of the archive
REF = os.path.join(ARCH, 'reference', 'ANSYS_18_5kW')
os.makedirs(OUT, exist_ok=True)
C = {'fe': '#1b4f72', 'neu': '#c0392b', 'op': '#e67e22', 'op2': '#f0b27a', 'cav': '#1e8449', 'carter': '#5d6d7e', 'grey': '#95a5a6', 'grad': '#7d3c98'}
KC_CARTER, KC_CARTER_EX = 1.266466, 1.267932

fem = json.load(open(f'{W}/prod_fem_results.json'))
op = json.load(open(f'{W}/prod_op_results.json'))
til = json.load(open(f'{W}/tiling_sweep_inf_iron.json'))
t4 = json.load(open(f'{W}/t4_graded_results.json'))


def save(fig, name):
    fig.savefig(f'{OUT}/{name}.png', dpi=400, bbox_inches='tight')
    fig.savefig(f'{OUT}/{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def fig1():
    """One slot pitch, with the true slot profiles of Table 1.

    Stator slot : opening 2.000 x 0.500 mm, wedge to 5.236 mm over 2.500 mm,
                  trapezoidal body to 8.472 mm over 24.724 mm.
    Rotor slot  : isthmus 2.000 x 1.000 mm, then a 5.804 mm circle centred
                  3.902 mm below the bore, joined by tangents of half-angle
                  6.44 deg to a 2.038 mm circle centred 16.790 mm deeper;
                  total depth 21.711 mm.
    The slot cross-sections are to scale in both directions and are cut by
    the edge of the window; the air gap is drawn wider than scale.
    """
    TAU = 10.719
    BS0, HS0, BS1, HS1, BS2, HS2 = 2.0, 0.5, 5.236, 2.5, 8.472, 24.724
    BR0, HR0 = 2.0, 1.0
    DR1, DR2, HR2 = 5.80425613768, 2.038213646116, 16.78981901905
    GAP, HS_WIN, HR_WIN, XOFF = 1.0, 4.6, 7.4, 2.2
    IRON, TILE = '#d5d8dc', '#aeb6bf'

    def stator_half(hmax):
        w = BS1 / 2 + (BS2 - BS1) / 2 * (hmax - HS0 - HS1) / HS2
        return [(BS0 / 2, 0.0), (BS0 / 2, HS0), (BS1 / 2, HS0 + HS1), (w, hmax)]

    def rotor_half(hmax, narc=70):
        r1, r2 = DR1 / 2, DR2 / 2
        c1 = HR0 + r1
        alpha = np.arcsin((r1 - r2) / HR2)
        t0, t1 = np.arcsin((BR0 / 2) / r1), np.pi / 2 + alpha
        pts = [(BR0 / 2, 0.0), (BR0 / 2, HR0)]
        for t in np.linspace(t0, t1, narc):
            y = c1 - r1 * np.cos(t)
            if y > hmax:
                break
            pts.append((r1 * np.sin(t), y))
        ytan, xtan = c1 + r1 * np.sin(alpha), r1 * np.cos(alpha)
        if hmax > ytan:
            pts.append((xtan - (hmax - ytan) * np.tan(alpha), hmax))
        return pts

    def slot(ax, half, xc, sign, y_shift):
        right = [(xc + x, sign * y + y_shift) for x, y in half]
        left = [(xc - x, sign * y + y_shift) for x, y in reversed(half)]
        ax.add_patch(Polygon(right + left, closed=True, fc='white', ec='none', zorder=2))
        ax.plot([p[0] for p in right], [p[1] for p in right], color='k', lw=0.7, zorder=3)
        ax.plot([p[0] for p in left], [p[1] for p in left], color='k', lw=0.7, zorder=3)

    def panel(ax, mode):
        col = C['neu'] if mode == 'insulating' else C['cav']
        x0, x1 = -0.7, TAU + 0.7
        ax.set_aspect('equal'); ax.axis('off')
        ax.add_patch(Rectangle((x0, GAP), x1 - x0, HS_WIN, fc=IRON, ec='k', lw=0.6))
        ax.add_patch(Rectangle((x0, -HR_WIN), x1 - x0, HR_WIN, fc=IRON, ec='k', lw=0.6))
        xs = TAU / 2
        slot(ax, stator_half(HS_WIN), xs, +1, GAP)
        ax.text(xs, GAP + HS_WIN - 0.9, 'stator slot', ha='center', va='center', fontsize=7, zorder=4)
        xr = xs + XOFF
        slot(ax, rotor_half(HR_WIN), xr, -1, 0.0)
        ax.text(xr, -HR_WIN + 1.35, 'rotor bar', ha='center', va='center', fontsize=7, zorder=4)

        h = 0.22

        def tile(xa, xb, y, n, c):
            for i in range(n):
                xx = xa + i * (xb - xa) / n
                ax.add_patch(Rectangle((xx, y), (xb - xa) / n, h, fc=c, ec='k', lw=0.25, zorder=5))
        nT, nO = 4, 3
        tile(x0, xs - BS0 / 2, GAP - h, nT, TILE)
        tile(xs + BS0 / 2, x1, GAP - h, nT, TILE)
        tile(xs - BS0 / 2, xs + BS0 / 2, GAP - h, nO, col)
        tile(x0, xr - BR0 / 2, 0.0, nT, TILE)
        tile(xr + BR0 / 2, x1, 0.0, nT, TILE)
        tile(xr - BR0 / 2, xr + BR0 / 2, 0.0, nO, col)

        for xx in np.linspace(xs - 3.4, xs + 3.4, 9):
            if xs - BS0 / 2 + 0.05 < xx < xs + BS0 / 2 - 0.05:
                if mode == 'insulating':
                    tgt = (xs - BS0 / 2 - 0.10, GAP - h) if xx < xs else (xs + BS0 / 2 + 0.10, GAP - h)
                    rad = 0.5 if xx < xs else -0.5
                else:
                    tgt = (xs - BS0 / 2 - 0.03, GAP + HS0 + 0.5) if xx < xs else (xs + BS0 / 2 + 0.03, GAP + HS0 + 0.5)
                    rad = -0.35 if xx < xs else 0.35
                ax.annotate('', xy=tgt, xytext=(xx, h + 0.04), zorder=6,
                            arrowprops=dict(arrowstyle='->', lw=0.8, color=col,
                                            connectionstyle='arc3,rad=%.2f' % rad))
            else:
                ax.annotate('', xy=(xx, GAP - h), xytext=(xx, h + 0.04), zorder=6,
                            arrowprops=dict(arrowstyle='->', lw=0.8, color=C['fe']))

        ax.plot([x0 + 0.10, x0 + 0.46], [h, h], color='k', lw=0.6, zorder=6)
        ax.plot([x0 + 0.10, x0 + 0.46], [GAP - h, GAP - h], color='k', lw=0.6, zorder=6)
        ax.plot([x0 + 0.28, x0 + 0.28], [h, GAP - h], color='k', lw=0.6, zorder=6)
        ax.text(x0 + 0.56, GAP / 2, r'$g$', ha='left', va='center', fontsize=8, zorder=6)
        ax.text(xs, GAP + 0.06, r'$b_0$', ha='center', va='bottom', fontsize=7.5, zorder=4)
        ax.text(x0, GAP + HS_WIN + 0.18, r'$n_T$ columns per tooth face', ha='left', va='bottom', fontsize=6.6)
        ax.text(x1, GAP + HS_WIN + 0.18, r'$n_O$ columns per opening', ha='right', va='bottom',
                fontsize=6.6, color=col)
        ax.text(x1 - 0.15, GAP + HS_WIN - 0.22, 'stator', ha='right', va='top', fontsize=7.5)
        ax.text(x0 + 0.15, -HR_WIN + 0.22, 'rotor', ha='left', va='bottom', fontsize=7.5)
        t = (r'(a) condition $\Phi_O=0$: the opening is a flux barrier' if mode == 'insulating'
             else r'(b) cavity admittance $\mathbf{Q}$: flux enters the slot walls')
        ax.set_title(t, fontsize=7.6, pad=14)
        ax.set_xlim(x0 - 0.05, x1 + 0.05)
        ax.set_ylim(-HR_WIN - 0.15, GAP + HS_WIN + 0.9)

    fig, axs = plt.subplots(1, 2, figsize=(6.8, 3.75))
    panel(axs[0], 'insulating')
    panel(axs[1], 'cavity')
    fig.subplots_adjust(wspace=0.10)
    save(fig, 'fig1_opening_conditions')


def fig2():
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    def curve(d, **kw):
        x = np.array([float(k) for k in d.keys()]); y = np.array(list(d.values()))
        o = np.argsort(x); x = x[o]; y = y[o]
        x = np.append(x, 1.0); y = np.append(y, y[0])
        ax.plot(x, y, **kw)
    curve(fem['pos_real'], color=C['fe'], marker='o', ms=3.5, lw=1.2, label='FE, real slots (reference)')
    curve(op['pos_cav'], color=C['cav'], marker='s', ms=3, lw=1.0, label='operator + slot cavities (33, 16)')
    curve(fem['pos_neu'], color=C['neu'], marker='^', ms=3.5, lw=1.2, ls='--', label='FE, openings as flux barriers')
    curve(op['pos_op_33_16'], color=C['op'], marker='v', ms=3, lw=1.0, ls='--', label=r'operator, $\Phi_O=0$ (33, 16)')
    curve(op['pos_op_17_4'], color=C['op2'], marker='D', ms=2.5, lw=0.9, ls=':', label=r'operator, $\Phi_O=0$ (17, 4)')
    ax.axhline(KC_CARTER_EX, color=C['carter'], lw=1.0, ls='-.')
    ax.text(0.02, KC_CARTER_EX + 0.006, "Carter (exact conformal form) 1.2679", fontsize=7, color=C['carter'])
    ax.set_xlabel(r'rotor position $\varphi / \tau_r$ (one rotor slot pitch)')
    ax.set_ylabel(r'slotting ratio $k_C = B_{g1}^{\rm smooth}/B_{g1}^{\rm slotted}$')
    ax.set_xlim(0, 1); ax.set_ylim(1.2, 1.65)
    ax.legend(loc='upper center', ncol=1, fontsize=6.8, handlelength=2.2)
    save(fig, 'fig2_kc_vs_position')


def fig3():
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.8), sharey=False)
    ax = axs[0]
    nOs = [2, 4, 8, 16]
    for nT, mk, col in zip([9, 17, 33, 65], ['o', 's', '^', 'D'], ['#f5b7b1', '#ec7063', '#cb4335', '#78281f']):
        y = [til[f'p1_{nT}_{n}']['kC'] for n in nOs]
        ax.plot(nOs, y, marker=mk, ms=3.5, lw=1, color=col, label=rf'$n_T={nT}$, symmetric hat')
    for nT, mk in zip([33, 65], ['^', 'D']):
        y = [til[f'p1a_{nT}_{n}']['kC'] for n in nOs]
        ax.plot(nOs, y, marker=mk, ms=3.5, lw=1, ls='--', mfc='none', color='#7b241c', label=rf'$n_T={nT}$, asymmetric hat')
    for nT, mk in zip([33, 65], ['^', 'D']):
        y = [t4['neumann_graded'][f'{nT}_{n}']['kC'] for n in nOs]
        ax.plot(nOs, y, marker=mk, ms=3.5, lw=1, ls=':', mfc='none', color=C['grad'], label=rf'$n_T={nT}$, graded tiling, asym. hat')
    ax.axhline(1.460, color=C['neu'], lw=1, ls='-.'); ax.text(2.1, 1.468, 'FE limit, flux-barrier openings (1.46)', fontsize=7, color=C['neu'])
    ax.axhline(1.2934, color=C['fe'], lw=1, ls='-.'); ax.text(2.1, 1.300, 'FE, real slots (1.293)', fontsize=7, color=C['fe'])
    ax.set_xscale('log', base=2); ax.set_xticks(nOs); ax.set_xticklabels(nOs)
    ax.set_xlabel(r'opening columns $n_O$'); ax.set_ylabel(r'$k_C$ at $\varphi = 0$, infinite iron')
    ax.set_title(r'(a) operator with $\Phi_O = 0$', fontsize=8.5); ax.set_ylim(1.25, 2.05)
    ax.legend(fontsize=5.6, loc='upper left', ncol=2, columnspacing=0.8, handlelength=1.8)
    ax = axs[1]
    nOs2 = [4, 8, 16, 32]
    for nT, mk, col in zip([17, 33, 65], ['s', '^', 'D'], ['#7dcea0', '#229954', '#145a32']):
        y = [op['cav_sweep'].get(f'{nT}_{n}') for n in nOs2]
        xx = [n for n, v in zip(nOs2, y) if v is not None]; yy = [v for v in y if v is not None]
        ax.plot(xx, yy, marker=mk, ms=3.5, lw=1, color=col, label=rf'$n_T={nT}$, symmetric hat')
    for nT, mk in zip([33, 65], ['^', 'D']):
        y = [t4['cav_graded'].get(f'{nT}_{n}') for n in nOs2]
        xx = [n for n, v in zip(nOs2, y) if v is not None]; yy = [v['kC'] for v in y if v is not None]
        ax.plot(xx, yy, marker=mk, ms=3.5, lw=1, ls=':', mfc='none', color=C['grad'], label=rf'$n_T={nT}$, graded tiling, asym. hat')
    ax.axhline(1.2934, color=C['fe'], lw=1, ls='-.'); ax.text(4.1, 1.2965, 'FE, real slots (1.293)', fontsize=7, color=C['fe'])
    ax.axhline(KC_CARTER_EX, color=C['carter'], lw=1, ls=':'); ax.text(4.1, 1.258, 'Carter 1.268', fontsize=7, color=C['carter'])
    ax.set_xscale('log', base=2); ax.set_xticks(nOs2); ax.set_xticklabels(nOs2)
    ax.set_xlabel(r'opening columns $n_O$'); ax.set_ylabel(r'$k_C$ at $\varphi = 0$, infinite iron')
    ax.set_title('(b) operator with slot cavities', fontsize=8.5); ax.set_ylim(1.25, 1.40)
    ax.legend(fontsize=5.8, loc='upper right')
    save(fig, 'fig3_tiling_convergence')


def fig4():
    d = np.load(f'{W}/field_waveforms.npz')
    th, fe, o, cv = d['th'], d['Br_fem'], d['Br_op'], d['Br_cav']
    deg = np.degrees(th)
    ker = np.ones(200) / 200
    sm = np.convolve(np.abs(fe), ker, mode='same')
    c = deg[np.argmax(sm)]
    c = 7.5 * np.round(c / 7.5)
    sgn = np.sign(fe[np.argmin(np.abs(deg - c))])
    x = deg - c
    x = np.where(x > 180, x - 360, x); x = np.where(x < -180, x + 360, x)
    o_ = np.argsort(x); sel = (x[o_] > -1.55 * 7.5) & (x[o_] < 1.55 * 7.5)
    fig, ax = plt.subplots(figsize=(6.4, 2.8))
    ax.plot(x[o_][sel], sgn * fe[o_][sel], color=C['fe'], lw=1.4, label='FE, real slot geometry')
    ax.plot(x[o_][sel], sgn * cv[o_][sel], color=C['cav'], lw=1.0, ls='--', label='operator + slot cavities (33, 16)')
    ax.plot(x[o_][sel], sgn * o[o_][sel], color=C['neu'], lw=1.0, ls=':', label=r'operator, $\Phi_O = 0$ (33, 16)')
    w = np.degrees(2.0e-3 / 81.8887)
    for k in [-1, 0, 1]:
        ax.axvspan(3.75 + 7.5 * k - w / 2, 3.75 + 7.5 * k + w / 2, color='#eaeded', zorder=0)
    ymax = 1.15 * np.max(sgn * fe[o_][sel])
    ax.text(3.75, 0.96 * ymax, 'stator opening', ha='center', va='top', fontsize=6.5, color='#7f8c8d')
    ax.set_ylim(min(0, 1.1 * np.min(sgn * fe[o_][sel])), ymax)
    ax.set_xlabel(r'mechanical angle from the pole centre (degrees), rotor at $\varphi = 0$')
    ax.set_ylabel(r'$B_r$ at mid-gap (T per A of $I_m$)')
    ax.set_xlim(-1.55 * 7.5, 1.55 * 7.5); ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), fontsize=7, ncol=3)
    save(fig, 'fig4_midgap_waveform_linear')


def fig5():
    T = np.loadtxt(f'{REF}/carat#U00e9ristique en fonction glissement/Torque Plot 2.tab', skiprows=1)
    I = np.loadtxt(f'{REF}/carat#U00e9ristique en fonction glissement/Winding Plot 1.tab', skiprows=1)
    b10 = sio.loadmat(f'{MECO}/B10_b1_skewoff.mat', squeeze_me=True)
    z5 = sio.loadmat(f'{MECO}/Z5_sweep_cav.mat', squeeze_me=True)
    sl, Tm, Im = b10['sl'], b10['T'], b10['I']
    sl2, Tc, Ic = z5['sl'], z5['T'], z5['I']
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.7))
    ax = axs[0]
    clean = T[:, 0] <= 0.125
    ax.plot(T[~clean, 0], T[~clean, 1], '.', color=C['grey'], ms=3, label='FE parametric sweep, noisy bands')
    ax.plot(T[clean, 0], T[clean, 1], '.', color=C['fe'], ms=4, label='FE parametric sweep, s ≤ 0.125')
    ax.plot([0.01883, 1.0], [121.63, 104.31], 'o', mfc='none', mec=C['fe'], ms=6, mew=1.2, label='FE dedicated 2 s transients')
    ax.plot(sl, Tm, '-', color=C['neu'], lw=1.2, label=r'network, operator $\Phi_O=0$ (17, 4)')
    ax.plot(sl2, Tc, '--', color=C['cav'], lw=1.2, label='network, operator + cavities (33, 16)')
    ax.set_xlabel('slip $s$'); ax.set_ylabel('electromagnetic torque (N m)'); ax.set_xlim(0, 1); ax.set_ylim(0, 400)
    ax.legend(fontsize=6.2, loc='upper right'); ax.set_title('(a)', fontsize=8.5)
    ax = axs[1]
    ax.plot(I[:, 0], I[:, 1], '.', color=C['fe'], ms=3, label='FE parametric sweep')
    ax.plot([0.0, 0.01883, 1.0], [8.499, 19.72, 104.4], 'o', mfc='none', mec=C['fe'], ms=6, mew=1.2, label='FE dedicated transients')
    ax.plot(sl, Im, '-', color=C['neu'], lw=1.2, label=r'network, $\Phi_O=0$ (17, 4)')
    ax.plot(sl2, Ic, '--', color=C['cav'], lw=1.2, label='network, cavities (33, 16)')
    ax.set_xlabel('slip $s$'); ax.set_ylabel('stator current (A rms)'); ax.set_xlim(0, 1); ax.set_ylim(0, 120)
    ax.legend(fontsize=6.2, loc='lower right'); ax.set_title('(b)', fontsize=8.5)
    save(fig, 'fig5_characteristics')


def fig6():
    fig, axs = plt.subplots(2, 2, figsize=(6.8, 4.4))
    for i, (cas, lab) in enumerate([('avide', 'no load'), ('charge', 'rated load')]):
        d = np.loadtxt(f'{MECO}/Z2_fields_{cas}.txt') if os.path.exists(f'{MECO}/Z2_fields_{cas}.txt') else np.loadtxt(f'{ARCH}/outputs/MEC_IM/Z2_fields_{cas}.txt')
        dc = np.loadtxt(f'{MECO}/Z4_fields_cav_{cas}.txt')
        th, Bre, Bte, Brm, Btm = d.T
        thc, Brec, Btec, Brc, Btc = dc.T
        deg = np.degrees(th)
        for j, (ye, ym, yc, comp) in enumerate([(Bre, Brm, Brc, r'$B_r$'), (Bte, Btm, Btc, r'$B_t$')]):
            ax = axs[j, i]
            sel = deg <= 90
            ax.plot(deg[sel], ye[sel], color=C['fe'], lw=0.9, label='FE (transient, single slice)')
            ax.plot(deg[sel], ym[sel], color=C['neu'], lw=0.7, ls=':', label=r'network, operator $\Phi_O=0$ (17, 4)')
            ax.plot(deg[sel], yc[sel], color=C['cav'], lw=0.8, ls='--', label='network, operator + cavities (33, 16)')
            ax.set_xlim(0, 90)
            if j == 1: ax.set_xlabel('mechanical angle (degrees)')
            ax.set_ylabel(comp + ' at mid-gap (T)')
            if j == 0: ax.set_title(f'({"ab"[i]}) {lab}', fontsize=8.5)
            if i == 0 and j == 0: ax.legend(fontsize=6.2, loc='lower left')
    fig.tight_layout()
    save(fig, 'fig6_machine_fields')


if __name__ == '__main__':
    which = sys.argv[1:] or ['1', '2', '3', '4', '5', '6']
    for k in which:
        globals()['fig' + k]()
    print('figures written to', OUT)
