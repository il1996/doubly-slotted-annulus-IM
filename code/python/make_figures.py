# -*- coding: utf-8 -*-
"""Figures of the revised manuscript, all from computed data."""
import json, numpy as np, os, scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Wedge, FancyArrowPatch
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.6,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False,
                     'axes.titlesize': 9, 'axes.labelsize': 8.5, 'legend.fontsize': 7.5})
OUT = '/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/figs'
os.makedirs(OUT, exist_ok=True)
W = '/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/work'
MECO = '/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/mec/MEC_IM/octave_out'
REF = '/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/repo/reference/ANSYS_18_5kW'
C = {'fe': '#1b4f72', 'neu': '#c0392b', 'op': '#e67e22', 'op2': '#f0b27a', 'cav': '#1e8449', 'carter': '#5d6d7e', 'grey': '#95a5a6'}
KC_CARTER, KC_CARTER_EX = 1.266466, 1.267932

fem = json.load(open(f'{W}/prod_fem_results.json'))
op = json.load(open(f'{W}/prod_op_results.json'))
til = json.load(open(f'{W}/tiling_sweep_inf_iron.json'))


def save(fig, name):
    fig.savefig(f'{OUT}/{name}.png', dpi=400, bbox_inches='tight')
    fig.savefig(f'{OUT}/{name}.pdf', bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- Figure 1: schematic
def fig1():
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 3.0))
    for ax, mode in zip(axs, ['insulating', 'cavity']):
        ax.set_aspect('equal'); ax.axis('off')
        g = 1.0; tau = 5.0; b0 = 1.2; face = tau - b0
        col = C['neu'] if mode == 'insulating' else C['cav']
        # iron
        ax.add_patch(Rectangle((-0.5, -2.6), tau + 1.0, 2.6, fc='#d5d8dc', ec='k', lw=0.6))
        ax.add_patch(Rectangle((-0.5, g), tau + 1.0, 2.8, fc='#d5d8dc', ec='k', lw=0.6))
        # stator slot centred at tau/2 : isthmus + wedge + body
        x0 = face / 2
        ax.add_patch(Polygon([(x0, g), (x0 + b0, g), (x0 + b0, g + 0.45), (x0 + b0 + 1.0, g + 1.4),
                              (x0 + b0 + 1.0, g + 2.8), (x0 - 1.0, g + 2.8), (x0 - 1.0, g + 1.4), (x0, g + 0.45)],
                             closed=True, fc='white', ec='k', lw=0.6))
        ax.text(tau / 2, g + 2.1, 'stator slot', ha='center', va='center', fontsize=7)
        # rotor slot, shifted
        xr = x0 + 0.9
        ax.add_patch(Rectangle((xr, -0.7), b0, 0.7, fc='white', ec='k', lw=0.6))
        ax.add_patch(Wedge((xr + b0 / 2, -1.55), 0.85, 0, 360, fc='white', ec='k', lw=0.6))
        ax.text(tau + 0.35, -2.35, 'rotor', ha='right', fontsize=7.5)
        ax.text(tau + 0.35, g + 2.6, 'stator', ha='right', va='top', fontsize=7.5)
        # tiling boxes just below the stator bore (in the gap), drawn thin
        nT, nO = 5, 3
        for i in range(nT):
            xx = -face / 2 + i * face / nT
            if xx + face / nT > -0.5:
                ax.add_patch(Rectangle((max(xx, -0.5), g - 0.22), face / nT - max(0, -0.5 - xx), 0.22, fc='#aeb6bf', ec='k', lw=0.3))
        for i in range(nT):
            xx = face / 2 + b0 + i * face / nT
            if xx < tau + 0.5:
                ax.add_patch(Rectangle((xx, g - 0.22), min(face / nT, tau + 0.5 - xx), 0.22, fc='#aeb6bf', ec='k', lw=0.3))
        for i in range(nO):
            ax.add_patch(Rectangle((x0 + i * b0 / nO, g - 0.22), b0 / nO, 0.22, fc=col, ec='k', lw=0.3))
        # field lines from the rotor surface
        xs = np.linspace(x0 - 0.9, x0 + b0 + 0.9, 8)
        for xx in xs:
            if x0 + 0.05 < xx < x0 + b0 - 0.05:
                if mode == 'insulating':
                    tgt = (x0 - 0.12, g - 0.22) if xx < x0 + b0 / 2 else (x0 + b0 + 0.12, g - 0.22)
                    rad = 0.4 if xx < x0 + b0 / 2 else -0.4
                else:
                    tgt = (x0 - 0.02, g + 0.75) if xx < x0 + b0 / 2 else (x0 + b0 + 0.02, g + 0.75)
                    rad = -0.3 if xx < x0 + b0 / 2 else 0.3
                ax.annotate('', xy=tgt, xytext=(xx, 0.0), arrowprops=dict(arrowstyle='->', lw=0.8, color=col, connectionstyle='arc3,rad=%.2f' % rad))
            else:
                ax.annotate('', xy=(xx, g - 0.22), xytext=(xx, 0.0), arrowprops=dict(arrowstyle='->', lw=0.8, color=C['fe']))
        # labels
        ax.text(-0.45, g / 2 - 0.12, r'$g$', ha='right', va='center', fontsize=8)
        ax.text(x0 + b0 / 2, g + 0.12, r'$b_0$', ha='center', va='bottom', fontsize=7.5)
        ax.text(-0.1, g + 0.35, r'$n_T$ columns per face', ha='left', va='bottom', fontsize=7)
        ax.text(x0 + b0 / 2, -0.95, r'$n_O$ columns per opening', ha='center', va='top', fontsize=7, color=col)
        if mode == 'insulating':
            ax.set_title(r'(a) condition $\Phi_O = 0$ : the opening is a flux barrier', fontsize=8)
        else:
            ax.set_title(r'(b) cavity admittance $\mathbf{Q}$ : flux enters the slot walls', fontsize=8)
        ax.set_xlim(-0.6, tau + 0.6); ax.set_ylim(-2.7, g + 2.9)
    save(fig, 'fig1_opening_conditions')


# ---------------------------------------------------------------- Figure 2: position dependence
def fig2():
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    def curve(d, **kw):
        x = np.array([float(k) for k in d.keys()]); y = np.array(list(d.values()))
        o = np.argsort(x); x = x[o]; y = y[o]
        x = np.append(x, 1.0); y = np.append(y, y[0])            # periodic
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


# ---------------------------------------------------------------- Figure 3: tiling convergence
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
    ax.axhline(1.460, color=C['neu'], lw=1, ls='-.'); ax.text(2.1, 1.468, 'FE limit, flux-barrier openings (1.46)', fontsize=7, color=C['neu'])
    ax.axhline(1.2934, color=C['fe'], lw=1, ls='-.'); ax.text(2.1, 1.300, 'FE, real slots (1.293)', fontsize=7, color=C['fe'])
    ax.set_xscale('log', base=2); ax.set_xticks(nOs); ax.set_xticklabels(nOs)
    ax.set_xlabel(r'opening columns $n_O$'); ax.set_ylabel(r'$k_C$ at $\varphi = 0$, infinite iron')
    ax.set_title(r'(a) operator with $\Phi_O = 0$', fontsize=8.5); ax.set_ylim(1.25, 1.85)
    ax.legend(fontsize=6.5, loc='upper right')
    ax = axs[1]
    nOs2 = [4, 8, 16, 32]
    for nT, mk, col in zip([17, 33, 65], ['s', '^', 'D'], ['#7dcea0', '#229954', '#145a32']):
        y = [op['cav_sweep'].get(f'{nT}_{n}') for n in nOs2]
        xx = [n for n, v in zip(nOs2, y) if v is not None]; yy = [v for v in y if v is not None]
        ax.plot(xx, yy, marker=mk, ms=3.5, lw=1, color=col, label=rf'$n_T={nT}$')
    ax.axhline(1.2934, color=C['fe'], lw=1, ls='-.'); ax.text(4.1, 1.2965, 'FE, real slots (1.293)', fontsize=7, color=C['fe'])
    ax.axhline(KC_CARTER_EX, color=C['carter'], lw=1, ls=':'); ax.text(4.1, 1.258, 'Carter 1.268', fontsize=7, color=C['carter'])
    ax.set_xscale('log', base=2); ax.set_xticks(nOs2); ax.set_xticklabels(nOs2)
    ax.set_xlabel(r'opening columns $n_O$'); ax.set_ylabel(r'$k_C$ at $\varphi = 0$, infinite iron')
    ax.set_title('(b) operator with slot cavities', fontsize=8.5); ax.set_ylim(1.25, 1.40)
    ax.legend(fontsize=6.5, loc='upper right')
    save(fig, 'fig3_tiling_convergence')


# ---------------------------------------------------------------- Figure 4: field waveforms (infinite iron)
def fig4():
    d = np.load(f'{W}/field_waveforms.npz')
    th, fe, o, cv = d['th'], d['Br_fem'], d['Br_op'], d['Br_cav']
    deg = np.degrees(th)
    # window of two stator slot pitches centred on the tooth carrying the largest flux (pole centre)
    ker = np.ones(200) / 200
    sm = np.convolve(np.abs(fe), ker, mode='same')
    c = deg[np.argmax(sm)]
    c = 7.5 * np.round(c / 7.5)                       # snap to a tooth centre
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


# ---------------------------------------------------------------- Figure 5: machine characteristics
def fig5():
    T = np.loadtxt(f'{REF}/carat#U00e9ristique en fonction glissement/Torque Plot 2.tab', skiprows=1)
    I = np.loadtxt(f'{REF}/carat#U00e9ristique en fonction glissement/Winding Plot 1.tab', skiprows=1)
    b10 = sio.loadmat('/mnt/user-data/uploads/MEC/MEC_IM/B10_b1_skewoff.mat', squeeze_me=True)
    sl, Tm, Im = b10['sl'], b10['T'], b10['I']
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.7))
    ax = axs[0]
    clean = T[:, 0] <= 0.125
    ax.plot(T[~clean, 0], T[~clean, 1], '.', color=C['grey'], ms=3, label='FE parametric sweep, noisy bands')
    ax.plot(T[clean, 0], T[clean, 1], '.', color=C['fe'], ms=4, label='FE parametric sweep, s ≤ 0.125')
    ax.plot([0.01883, 1.0], [121.63, 104.31], 'o', mfc='none', mec=C['fe'], ms=6, mew=1.2, label='FE dedicated 2 s transients')
    ax.plot(sl, Tm, '-', color=C['cav'], lw=1.2, label=r'network, operator $\Phi_O=0$ (17, 4)')
    ax.set_xlabel('slip $s$'); ax.set_ylabel('electromagnetic torque (N m)'); ax.set_xlim(0, 1); ax.set_ylim(0, 400)
    ax.legend(fontsize=6.5, loc='upper right'); ax.set_title('(a)', fontsize=8.5)
    ax = axs[1]
    ax.plot(I[:, 0], I[:, 1], '.', color=C['fe'], ms=3, label='FE parametric sweep')
    ax.plot([0.0, 0.01883, 1.0], [8.499, 19.72, 104.4], 'o', mfc='none', mec=C['fe'], ms=6, mew=1.2, label='FE dedicated transients')
    ax.plot(sl, Im, '-', color=C['cav'], lw=1.2, label='network')
    ax.set_xlabel('slip $s$'); ax.set_ylabel('stator current (A rms)'); ax.set_xlim(0, 1); ax.set_ylim(0, 120)
    ax.legend(fontsize=6.5, loc='lower right'); ax.set_title('(b)', fontsize=8.5)
    save(fig, 'fig5_characteristics')


# ---------------------------------------------------------------- Figure 6: machine mid-gap field
def fig6():
    fig, axs = plt.subplots(2, 2, figsize=(6.8, 4.4))
    for i, (cas, lab) in enumerate([('avide', 'no load'), ('charge', 'rated load')]):
        d = np.loadtxt(f'{MECO}/Z2_fields_{cas}.txt')
        th, Bre, Bte, Brm, Btm = d.T
        deg = np.degrees(th)
        for j, (ye, ym, comp) in enumerate([(Bre, Brm, r'$B_r$'), (Bte, Btm, r'$B_t$')]):
            ax = axs[j, i]
            sel = deg <= 90
            ax.plot(deg[sel], ye[sel], color=C['fe'], lw=0.9, label='FE (transient, single slice)')
            ax.plot(deg[sel], ym[sel], color=C['cav'], lw=0.8, ls='--', label=r'network, operator $\Phi_O=0$')
            ax.set_xlim(0, 90)
            if j == 1: ax.set_xlabel('mechanical angle (degrees)')
            ax.set_ylabel(comp + ' at mid-gap (T)')
            if j == 0: ax.set_title(f'({"ab"[i]}) {lab}', fontsize=8.5)
            if i == 0 and j == 0: ax.legend(fontsize=6.5, loc='lower left')
    fig.tight_layout()
    save(fig, 'fig6_machine_fields')


if __name__ == '__main__':
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
    print('figures written to', OUT)
