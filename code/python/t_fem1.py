import numpy as np, time, sys
import fem_slots as F
from dsop import Machine, MU0
M=Machine()
# 1) smooth check: thin current sheet, no rotor slots
F.HS0_save=(F.BS0,F.HS0,F.BS1,F.HS1,F.BS2,F.HS2)
def set_thin():
    F.BS0,F.HS0,F.BS1,F.HS1,F.BS2,F.HS2 = 2.0,0.02,2.0,0.02,2.0,0.06
def restore():
    F.BS0,F.HS0,F.BS1,F.HS1,F.BS2,F.HS2 = F.HS0_save
set_thin()
for lc in [0.12,0.06]:
    r=F.run_case(M,0.0,lc_gap=lc,rotor_slots=False)
    print(f"SMOOTH CHECK lc={lc}: nodes {r['nnodes']} Bg1 FEM {r['Bg1']:.6f} vs analytic {r['Bg1_smooth']:.6f} -> ratio {r['Bg1']/r['Bg1_smooth']:.5f}  ({r['t']:.0f}s)", flush=True)
restore()
# 2) mesh convergence, slotted, phi=0
for lc in [0.12,0.08,0.06,0.045]:
    r=F.run_case(M,0.0,lc_gap=lc)
    print(f"SLOTTED lc={lc}: nodes {r['nnodes']} tri {r['ntri']}  Bg1 {r['Bg1']:.6f}  kC {r['kC']:.5f}  ({r['t']:.0f}s)", flush=True)
