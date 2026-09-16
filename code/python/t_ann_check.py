import numpy as np, time
import fem_annulus as FA
from dsop import Machine, MU0
M=Machine()
i3=np.sqrt(2)*np.array([1.0,-0.5,-0.5]); Us=M.tooth_mmf(i3)
for lc in [0.08,0.05]:
    t=time.time()
    xyz,tri,arcs,fs,fr=FA.build_annulus_mesh(M,0.0,lc)
    # Dirichlet everywhere: faces at U_j, stator openings at the mean of the two neighbours, rotor all 0
    Uop=0.5*(Us+np.roll(Us,-1))
    phi,R=FA.solve_scalar(M,xyz,tri,arcs,fs,fr,Us,0.0,opening='dirichlet',Uopen_s=Uop,Uopen_r=np.zeros(M.Nr))
    Bn=FA.midgap_harmonics_scalar(M,xyz,tri,phi)
    # analytic: piecewise constant on the 96 stator arcs
    tau=2*np.pi/M.Ns; a_o=M.bs0/M.Rs; a_f=tau-a_o
    edges=[]; vals=[]
    for k in range(M.Ns):
        c=tau*k; edges+= [c-a_f/2, c+a_f/2]; vals+=[Us[k], Uop[k]]
    nn=np.arange(1,200)
    cs=FA.piecewise_constant_fourier(edges,vals,nn); cr=np.zeros_like(cs)
    Ban=FA.annulus_Bn_from_surface_coeffs(M,cs,cr,nn.astype(float))
    print(f"lc={lc} nodes={len(xyz)} ({time.time()-t:.0f}s)   n:  FEM  analytic  ratio")
    for n in [2,10,14,22,26,46,50]:
        print(f"   {n:3d}  {Bn[n]:.6f}  {Ban[n-1]:.6f}  {Bn[n]/Ban[n-1]:.5f}")
    # Neumann openings, same mesh
    phi,R=FA.solve_scalar(M,xyz,tri,arcs,fs,fr,Us,0.0,opening='neumann')
    BnN=FA.midgap_harmonics_scalar(M,xyz,tri,phi)
    # smooth reference: pure staircase (jumps at slot centres) analytic
    th_slots=2*np.pi*(np.arange(M.Ns)+0.5)/M.Ns
    edges2=list(th_slots-tau); vals2=list(Us)   # pitch j from slot j-1 centre to slot j centre
    cs2=FA.piecewise_constant_fourier(edges2,vals2,nn); Ban2=FA.annulus_Bn_from_surface_coeffs(M,cs2,cr,nn.astype(float))
    print(f"   NEUMANN openings: Bg1 = {BnN[2]:.6f}; smooth staircase analytic {Ban2[1]:.6f}; kC(Neumann) = {Ban2[1]/BnN[2]:.5f}")
