import numpy as np
from dsop import *
M=Machine(); tau_r=2*np.pi/M.Nr
print("operator (Phi_O=0), (17,4), Nh=8192, p1, infinite iron, vs rotor position:")
vals=[]
for frac in np.arange(0,1,0.125):
    r=slotting_ratio(M,17,4,8192,'p1',phi=frac*tau_r); vals.append(r['kC'])
    print(f"  phi/tau_r={frac:.3f}: kC={r['kC']:.5f}  Bg1_slot={r['Bg1_slot']:.5f}")
print("  mean over positions:", np.mean(vals), " min/max", min(vals), max(vals))
print("operator (33,8) at 3 positions:")
for frac in [0,0.25,0.5]:
    r=slotting_ratio(M,33,8,8192,'p1',phi=frac*tau_r); print(f"  phi/tau_r={frac:.3f}: kC={r['kC']:.5f}")
