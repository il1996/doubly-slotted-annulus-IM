# -*- coding: utf-8 -*-
"""Fig. 3 of the body (file fig3_b0g_sweep ; the tiling-convergence figure, file fig3_tiling_convergence, is Fig. S1 of the supplement since 20 Sept. 2026) -- sweep of the opening-to-gap ratio b_0/g on the linear chain
(prod_sweep_b0g.py -> outputs/python/prod_sweep_b0g_results.json).

One panel, 3.375 in wide, <= 2.0 in high.  Abscissa b_0/g ; ordinate the
deviation from the finite-element reference (formulation A, real slots,
twelve-position mean) of Carter's product (exact conformal form), of the
condensed operator with Phi_O = 0 and of the operator coupled to the slot
cavities, tiling (33, 16), symmetric hat, N_h = 8192.  The bars at zero are the
min-max modulation of the reference over the twelve rotor positions, relative
to its mean.  The FE rows of b0/g = 12 and 16 are at the halved mesh lc 0.0225 (G5 ;\nprod_sweep_b0g_meshrow.py).  Colours and fonts are those of make_figures_v2.py ; identity is
also carried by marker shape and direct labels.
Run from code/python/ ; writes code/article/figures/fig3_b0g_sweep.{pdf,png}."""
import os, sys, json
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.6,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False,
                     'axes.labelsize': 8, 'legend.fontsize': 6.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
HERE = os.path.dirname(os.path.abspath(__file__))
ARCH = os.path.normpath(os.path.join(HERE, '..', '..'))
W = os.path.join(ARCH, 'outputs', 'python'); OUT = os.path.join(ARCH, 'code', 'article', 'figures')
C = {'fe': '#1b4f72', 'op': '#e67e22', 'cav': '#1e8449', 'carter': '#5d6d7e'}

RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(W, 'prod_sweep_b0g_results.json')     # optional: another results file
if len(sys.argv) > 2: OUT = sys.argv[2]                                                         # optional: another output folder
R = json.load(open(RES))
ALL = dict(R['sweep']); ALL['nominal'] = R['nominal']
pts = sorted(ALL.values(), key=lambda b: b['b0_over_g'])
x = np.array([b['b0_over_g'] for b in pts]); m = np.array([b['fe']['mean'] for b in pts])
dev = lambda k: 100 * (np.array(k) / m - 1)
y_carter = dev([b['carter']['k_exact'] for b in pts]); y_neu = dev([b['op_neu']['mean'] for b in pts]); y_cav = dev([b['op_cav']['mean'] for b in pts])
lo = dev([b['fe']['min'] for b in pts]); hi = dev([b['fe']['max'] for b in pts])
inom = int(np.argmin(np.abs(x - R['nominal']['b0_over_g'])))

fig, ax = plt.subplots(figsize=(3.375, 2.0))
ax.axhline(0, color=C['fe'], lw=0.8, zorder=1)
ax.errorbar(x, np.zeros_like(x), yerr=[-lo, hi], fmt='none', ecolor=C['fe'], elinewidth=0.8, capsize=2.5, capthick=0.8, zorder=2,
            label='FE reference: mean = 0, bars = min–max over position')
ax.plot(x, y_neu, color=C['op'], marker='v', ms=3.5, lw=1.0, ls='--', zorder=3, label=r'operator, $\Phi_O = 0$, (33, 16)')
ax.plot(x, y_cav, color=C['cav'], marker='s', ms=3.2, lw=1.0, zorder=4, label='operator + slot cavities, (33, 16)')
ax.plot(x, y_carter, color=C['carter'], marker='o', ms=3.2, lw=1.0, ls='-.', zorder=3, label="Carter's product, exact conformal form")
ymin = min(lo.min(), y_carter.min(), y_cav.min()) - 1.5; ymax = max(y_neu.max(), hi.max()) + 1.5
ax.set_ylim(ymin, ymax)
ax.axvline(x[inom], color='#333333', lw=0.6, ls=(0, (1.5, 2.5)), zorder=1)                      # the machine of the paper, b0/g = 8.23
ax.text(x[inom] + 0.25, ymin + 0.35, 'this machine', fontsize=6.2, color='#333333', ha='left', va='bottom')
ax.set_xlabel(r'slot opening to gap ratio $b_0/g$ (linear scale)')
ax.set_ylabel('deviation from FE mean (%)')
ax.set_xlim(0, max(x) + 1.5)
ax.set_xticks([2, 4, 8.2, 12, 16] if min(x) < 3 else [4, 8.2, 12, 16]); ax.set_xticklabels(['2', '4', '8.2', '12', '16'] if min(x) < 3 else ['4', '8.2', '12', '16'])
ax.set_xticks(np.arange(1, int(max(x)) + 2), minor=True); ax.tick_params(axis='x', which='minor', length=2, width=0.5)
ax.grid(axis='y', color='#dddddd', lw=0.4, zorder=0)
ax.legend(loc='center left', bbox_to_anchor=(0.0, 0.58), handlelength=2.0, borderaxespad=0.3, labelspacing=0.25)   # the band between the near-zero curves and the Phi_O = 0 curve
fig.tight_layout(pad=0.2)
os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, 'fig3_b0g_sweep.png'), dpi=400)
fig.savefig(os.path.join(OUT, 'fig3_b0g_sweep.pdf'))
print('fig3 written: %.3f x %.3f in ; points b0/g = %s' % (fig.get_figwidth(), fig.get_figheight(), np.round(x, 3)))
for b, yc, yn, yk in zip(pts, y_cav, y_neu, y_carter):
    print('  b0/g %7.4f : Carter ex %+7.3f %%  Phi_O=0 %+7.3f %%  cavities %+7.3f %%  FE band %+6.3f..%+6.3f %%' % (b['b0_over_g'], yk, yn, yc, 100 * (b['fe']['min'] / b['fe']['mean'] - 1), 100 * (b['fe']['max'] / b['fe']['mean'] - 1)))
