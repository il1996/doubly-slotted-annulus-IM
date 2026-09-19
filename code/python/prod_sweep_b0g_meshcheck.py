# -*- coding: utf-8 -*-
"""Production: mesh check of the b0/g sweep at the gaps that guard G5 of
prod_sweep_b0g.py does not cover.  G5 validates the finite-element mean at the
finest gap only (b0/g = 16, where the production mesh lc_gap = 0.045 mm spans
2.8 elements across the gap and had to be halved).  This script halves lc once
at the other computed gaps and reports the shift of the ratio (5):
  b0/g = 12 and 4 : the twelve rotor positions, mean shift and per-position shift;
  b0/g = 2        : phi = 0 only (325 000 nodes at lc 0.045 ; the halved mesh
                    exceeds a million nodes and is run once), declared as such.
Same construction as prod_sweep_b0g.py (fem_slots.run_case, closed-form
comparator recomputed at each g).  Reads outputs/python/prod_sweep_b0g_results.json
for the production-mesh values.  Writes prod_sweep_b0g_meshcheck_results.json and
prod_sweep_b0g_meshcheck_out.txt.  Run from code/python/."""
import os, sys, time, json, datetime
sys.dont_write_bytecode = True
import numpy as np
from dsop import Machine
import fem_slots as F
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.normpath(os.path.join(HERE, '..', '..', 'outputs', 'python'))
B0 = 2.0e-3; LC = 0.045
def now(): return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
def say(s): print(s, flush=True)
def machine_at_gap(g):
    M = Machine(); M.g = g; M.Rr = M.Rs - g; M.Dr = 2 * M.Rr; M.X = np.log(M.Rs / M.Rr); M.Rm = 0.5 * (M.Rs + M.Rr); return M
R = json.load(open(os.path.join(W, 'prod_sweep_b0g_results.json'))); S = R['sweep']
M0 = Machine(); tau_r = 2 * np.pi / M0.Nr; T0 = time.time()
PLAN = [(12.0, 'all'), (4.0, 'all'), (2.0, 'phi0')]
out = dict(date=now(), lc_production=LC, lc_halved=LC / 2, checks={}); L = ['prod_sweep_b0g_meshcheck.py -- run of %s ; lc_gap %.3f -> %.4f mm at the gaps not covered by G5 of prod_sweep_b0g.py' % (now(), LC, LC / 2)]
for r, mode in PLAN:
    key = '%g' % r
    if key not in S: say('  b0/g = %g absent from the sweep results -- skipped' % r); continue
    g = B0 / r; M = machine_at_gap(g); prod = S[key]['fe']; lc_prod = prod['lc']
    if abs(lc_prod - LC) > 1e-12: say('  b0/g = %g : production row already at lc %.5f (G5) -- skipped' % (r, lc_prod)); continue
    fr = sorted(float(k) for k in prod['per_position']); fr = fr if mode == 'all' else [0.0]
    t = time.time(); new = {}; nodes = []
    for f in fr:
        rc = F.run_case(M, f * tau_r, lc_gap=LC / 2); new['%.6f' % f] = float(rc['kC']); nodes.append(int(rc['nnodes']))
        say('      FE  b0/g=%g lc=%.5f phi/tau_r=%.4f : kC %.6f (lc %.3f: %.6f, shift %+.4f %%)  nodes %d  (%.0f s)' % (r, LC / 2, f, rc['kC'], LC, prod['per_position']['%.6f' % f], 100 * (rc['kC'] / prod['per_position']['%.6f' % f] - 1), rc['nnodes'], time.time() - t))
    shifts = {k: 100 * (new[k] / prod['per_position'][k] - 1) for k in new}
    rec = dict(mode=mode, positions=fr, per_position_halved=new, per_position_production={k: prod['per_position'][k] for k in new}, shift_pct=shifts, nodes_min=min(nodes), nodes_max=max(nodes), duration_s=time.time() - t)
    if mode == 'all':
        m_new = float(np.mean(list(new.values()))); rec.update(mean_halved=m_new, mean_production=prod['mean'], mean_shift_pct=100 * (m_new / prod['mean'] - 1))
        L.append('b0/g = %6.3f (g = %.6f mm) : twelve positions ; FE mean lc %.3f = %.6f, lc %.4f = %.6f, shift %+.4f %% ; per-position shifts %+.4f .. %+.4f %% ; nodes %d-%d ; %.0f s' % (
            r, g * 1e3, LC, prod['mean'], LC / 2, m_new, rec['mean_shift_pct'], min(shifts.values()), max(shifts.values()), rec['nodes_min'], rec['nodes_max'], rec['duration_s']))
    else:
        L.append('b0/g = %6.3f (g = %.6f mm) : phi = 0 only ; FE lc %.3f = %.6f, lc %.4f = %.6f, shift %+.4f %% ; nodes %d ; %.0f s (the twelve-position check at this gap would cost twelve such solves and was not run)' % (
            r, g * 1e3, LC, prod['per_position']['0.000000'], LC / 2, new['0.000000'], shifts['0.000000'], rec['nodes_min'], rec['duration_s']))
    out['checks'][key] = rec; json.dump(out, open(os.path.join(W, 'prod_sweep_b0g_meshcheck_results.json'), 'w'), indent=1)
g16 = R['guards']['G5']
L.append('for the record, G5 of prod_sweep_b0g.py at b0/g = %g : lc %s -> means %s, shifts %s %%, validated lc %s' % (g16['b0_over_g'], g16['lcs'], ['%.6f' % m for m in g16['means']], ['%+.4f' % s for s in g16['shifts_pct']], g16['validated_lc']))
L.append('nominal (Table 2 of the paper): lc 0.045 -> 0.028 moves k_C at phi = 0 by %+.4f %% (conv_real)' % (100 * (R['guards']['G1']['phi0_finest_mesh_lc0028'] / R['guards']['G1']['phi0_archived'] - 1)))
L.append('end %s ; total %.0f s' % (now(), time.time() - T0))
open(os.path.join(W, 'prod_sweep_b0g_meshcheck_out.txt'), 'w', encoding='utf-8').write('\n'.join(L) + '\n'); say('\n'.join(L))
