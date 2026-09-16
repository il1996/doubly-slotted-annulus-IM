import numpy as np, time, json
from dsop import *
M=Machine()
res={}
for basis in ['p1','p1a']:
    print("=== basis",basis)
    for nT in [9,17,33,65]:
        for nO in [2,4,8,16]:
            t=time.time()
            r=slotting_ratio(M,nT,nO,8192,basis)
            res[f"{basis}_{nT}_{nO}"]=dict(kC=r['kC'],Xm_slot=r['Xm_slot'],Xm_smooth=r['Xm_smooth'],ncol=r['ncol'])
            print(f"  nT={nT:3d} nO={nO:2d} cols={r['ncol']:5d}  kC={r['kC']:.6f}  Xm_slot={r['Xm_slot']:.4f}  Xm_smooth={r['Xm_smooth']:.4f}  ({time.time()-t:.1f}s)", flush=True)
json.dump(res,open('tiling_sweep_inf_iron.json','w'),indent=1)
