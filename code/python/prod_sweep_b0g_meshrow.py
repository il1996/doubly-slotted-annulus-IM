# -*- coding: utf-8 -*-
"""Production: promote the FE row of b0/g = 12 of the sweep to the halved mesh.

prod_sweep_b0g_meshcheck.py showed that halving lc_gap at b0/g = 12 moves the
twelve-position FE mean by -0.123 %, beyond the 0.1 % criterion that guard G5
applies at the finest gap.  This script applies to b0/g = 12 what G5 applied
to b0/g = 16: the twelve positions at lc 0.0225 (checks['12'] of
prod_sweep_b0g_meshcheck_results.json, same construction, same comparator)
become the FE row of that gap, and the production-mesh row is kept under
'fe_production_lc'.  It rewrites prod_sweep_b0g_results.json in place and
Tables 1 and 2 of prod_sweep_b0g_out.txt with an added line for b0/g = 12 at
lc 0.0225, the former line kept and marked 'production mesh'; a note is added
to the transcript header.  It then re-reads the bounds quoted in the text of
Section IV-C against the promoted table, digit by digit.
Run from code/python/ after prod_sweep_b0g.py and prod_sweep_b0g_meshcheck.py."""
import os, sys, json, datetime
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.normpath(os.path.join(HERE, '..', '..', 'outputs', 'python'))
RJ = os.path.join(W, 'prod_sweep_b0g_results.json'); RT = os.path.join(W, 'prod_sweep_b0g_out.txt'); MJ = os.path.join(W, 'prod_sweep_b0g_meshcheck_results.json')
now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
R = json.load(open(RJ)); Mc = json.load(open(MJ)); c12 = Mc['checks']['12']
assert c12['mode'] == 'all' and len(c12['per_position_halved']) == 12
S12 = R['sweep']['12']
if 'fe_production_lc' in S12:
    print('b0/g = 12 already promoted (fe at lc %.5f) -- nothing to do' % S12['fe']['lc']); sys.exit(0)
per = c12['per_position_halved']; v = np.array([per[k] for k in sorted(per, key=float)])
S12['fe_production_lc'] = S12['fe']
S12['fe'] = dict(lc=Mc['lc_halved'], per_position=per, Bg1_smooth_closed=S12['fe_production_lc']['Bg1_smooth_closed'], nodes_min=c12['nodes_min'], nodes_max=c12['nodes_max'],
                 mean=float(v.mean()), min=float(v.min()), max=float(v.max()), promoted_by='prod_sweep_b0g_meshrow.py ' + now, source='prod_sweep_b0g_meshcheck_results.json checks[12]')
R.setdefault('notes', []).append('%s : FE row of b0/g = 12 promoted to lc %.4f (prod_sweep_b0g_meshrow.py, from prod_sweep_b0g_meshcheck_results.json) ; production-mesh row kept under fe_production_lc' % (now, Mc['lc_halved']))
json.dump(R, open(RJ, 'w'), indent=1)
# ---- transcript: header note, Table 1 and Table 2 rows -------------------------------
L = open(RT, encoding='utf-8').read().split('\n')
dev = lambda a, b: 100 * (a / b - 1)
b = S12; p = S12['fe_production_lc']
row1 = '%8.4f %10.6f %8.5f | %10.6f %10.6f %10.6f | %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %12.6f %12.6f' % (
    b['b0_over_g'], b['g_mm'], b['fe']['lc'], b['fe']['mean'], b['fe']['min'], b['fe']['max'], b['carter']['k'], b['carter']['k_exact'], b['op_neu']['mean'], b['op_neu']['min'], b['op_neu']['max'],
    b['op_cav']['mean'], b['op_cav']['min'], b['op_cav']['max'], b['fe']['Bg1_smooth_closed'] * 1e3, b['op_neu']['Bg1_smooth_assembled'] * 1e3)
m = b['fe']['mean']
row2 = '%8.4f %10.6f | %+12.3f %+12.3f %+12.3f %+12.3f | %+12.3f %+12.3f | %+12.3f %+12.3f' % (
    b['b0_over_g'], b['g_mm'], dev(b['carter']['k'], m), dev(b['carter']['k_exact'], m), dev(b['op_neu']['mean'], m), dev(b['op_cav']['mean'], m),
    dev(b['fe']['min'], m), dev(b['fe']['max'], m), dev(b['op_cav']['min'], b['op_cav']['mean']), dev(b['op_cav']['max'], b['op_cav']['mean']))
out = []; t = None
for ln in L:
    if ln.startswith('TABLE 1'): t = 1
    elif ln.startswith('TABLE 2'): t = 2
    elif ln.startswith('TABLE 3'): t = 3
    if ln.startswith(' 12.0000') and t in (1, 2):
        out.append((row1 if t == 1 else row2) + '   <- lc 0.0225 (promoted, prod_sweep_b0g_meshrow.py)'); out.append(ln + '   <- production mesh lc 0.045, superseded'); continue
    out.append(ln)
    if ln.startswith('durations (s):'):
        out.append('note %s : the FE row of b0/g = 12 is promoted to lc 0.0225 (twelve positions of prod_sweep_b0g_meshcheck_results.json, shift of the mean -0.1227 %%), as G5 did for b0/g = 16 ; the production-mesh row is kept below it, marked' % now)
open(RT, 'w', encoding='utf-8').write('\n'.join(out))
print('promoted: b0/g = 12 FE mean %.6f (lc 0.045) -> %.6f (lc %.4f), min %.6f max %.6f' % (p['mean'], b['fe']['mean'], b['fe']['lc'], b['fe']['min'], b['fe']['max']))
# ---- the bounds quoted in Section IV-C, digit by digit, before and after -------------
ALL = dict(R['sweep']); ALL['nominal'] = R['nominal']; pts = sorted(ALL.values(), key=lambda x: x['b0_over_g'])
def bounds(use_prod12):
    rows = []
    for x in pts:
        fe = x['fe_production_lc'] if (use_prod12 and x['b0_over_g'] == 12.0) else x['fe']; mm = fe['mean']
        rows.append(dict(r=x['b0_over_g'], cex=dev(x['carter']['k_exact'], mm), capx=dev(x['carter']['k'], mm), neu=dev(x['op_neu']['mean'], mm), cav=dev(x['op_cav']['mean'], mm), lo=dev(fe['min'], mm), hi=dev(fe['max'], mm)))
    return rows
for lab, rows in [('before (12 at lc 0.045)', bounds(True)), ('after  (12 at lc 0.0225)', bounds(False))]:
    cav = [r['cav'] for r in rows]; neu = [r['neu'] for r in rows]
    print('%s : max|Carter ex| %.3f %% ; max|Carter apx| %.3f %% ; cavities %+.3f %% (b0/g=2) -> %+.3f %% (16), monotone %s, inside FE band at every g %s ; FE band +-%.3f %% (2) -> +-%.3f %% (16) ; Phi_O=0 min %+.3f %% (b0/g=%g) max %+.3f %% (b0/g=%g)' % (
        lab, max(abs(r['cex']) for r in rows), max(abs(r['capx']) for r in rows), cav[0], cav[-1], all(cav[i] < cav[i + 1] for i in range(len(cav) - 1)), all(r['cav'] < r['hi'] for r in rows),
        max(abs(rows[0]['lo']), rows[0]['hi']), max(abs(rows[-1]['lo']), rows[-1]['hi']), min(neu), rows[int(np.argmin(neu))]['r'], max(neu), rows[int(np.argmax(neu))]['r']))
    print('   row 12 : Carter ex %+.3f  Carter apx %+.3f  Phi_O=0 %+.3f  cavities %+.3f  FE band %+.3f..%+.3f' % tuple(next(r for r in rows if r['r'] == 12.0)[k] for k in ['cex', 'capx', 'neu', 'cav', 'lo', 'hi']))
print('text bounds quoted: 0.11 % (Carter exact), 0.33 % (approximate), +0.27 % -> +0.87 % (cavities), +-0.6 % -> +-3.2 % (FE band), 10.8 to 16.8 % (Phi_O = 0, largest at b0/g = 4)')
