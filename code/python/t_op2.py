import numpy as np, time
from dsop import *
M=Machine()
t=time.time()
r=slotting_ratio(M,17,4,8192,'p1')
print(f"(17,4) Nh=8192 p1: kC={r['kC']:.6f}  Xm_smooth={r['Xm_smooth']:.4f}  Xm_slot={r['Xm_slot']:.4f}  Bg1 smooth {r['Bg1_smooth']:.5f} slot {r['Bg1_slot']:.5f}  cols {r['ncol']}  t={time.time()-t:.1f}s")
for Nh in [512,1024,2048,4096]:
    r=slotting_ratio(M,17,4,Nh,'p1'); print(f"  Nh={Nh}: kC={r['kC']:.6f} Xm_slot={r['Xm_slot']:.4f} Xm_smooth={r['Xm_smooth']:.4f}")
print("p0 basis, truncation sweep (17,4):")
for Nh in [512,1024,2048,4096,8192]:
    r=slotting_ratio(M,17,4,Nh,'p0'); print(f"  Nh={Nh}: kC={r['kC']:.6f} Xm_slot={r['Xm_slot']:.4f} Xm_smooth={r['Xm_smooth']:.4f}")
