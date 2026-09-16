import numpy as np, time
from dsop import *
M=Machine()
# quick checks at small tiling: invariants, limits, k_C
for basis in ['p0','p1','p1a']:
    t=time.time()
    ss=tile_surface(M.Ns,M.Rs,M.bs0,6,2); sr=tile_surface(M.Nr,M.Rr,M.br0,6,2)
    A,proj=assemble(M,ss,sr,2048,basis)
    inv=invariants(A,proj,M)
    ev=np.linalg.eigvalsh(A)
    print(f"{basis}: cols {A.shape[0]}  I1 {inv['I1']:.2e}  I4 {inv['I4']:.2e} (Lam {inv['Lam_meas']:.6e} vs {inv['Lam0']:.6e})  sym {inv['sym']:.1e}  eig min {ev[0]:.2e} max {ev[-1]:.3e}  t={time.time()-t:.1f}s")
# uniform grid test for p1
for basis in ['p0','p1','p1a']:
    ss=tile_surface(M.Ns,M.Rs,0.0,8,0); sr=tile_surface(M.Nr,M.Rr,0.0,8,0)   # uniform
    A,proj=assemble(M,ss,sr,2048,basis); inv=invariants(A,proj,M)
    print(f"{basis} uniform: I1 {inv['I1']:.2e} I4 {inv['I4']:.2e}")
