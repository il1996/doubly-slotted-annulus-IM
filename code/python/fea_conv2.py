import numpy as np, os
for sub in ["transitoire/en charge","transitoire/a vide"]:
    REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'reference', 'ANSYS_18_5kW', *sub.split('/'))
    load=lambda p: np.loadtxt(p, skiprows=1)
    V=load(os.path.join(REF,"Winding Plot 3.tab")); I=load(os.path.join(REF,"Winding Plot 4.tab")); E=load(os.path.join(REF,"Winding Plot 2.tab"))
    t=I[:,0]; dt=t[1]-t[0]; m=(t>=1.0)&(t<=2.0)
    print("====",sub)
    print(" sum of currents rms:", np.sqrt(np.mean(np.sum(I[m,1:],axis=1)**2)))
    vn=np.mean(V[:,1:]-E[:,1:],axis=1)   # candidate neutral shift
    print(" neutral shift rms:", np.sqrt(np.mean(vn[m]**2)))
    R=0.4457
    for k in range(3):
        dI=np.gradient(I[:,1+k],dt)
        y=V[m,1+k]-vn[m]-E[m,1+k]
        A=np.vstack([I[m,1+k],dI[m]]).T
        coef,_,_,_=np.linalg.lstsq(A,y,rcond=None); fit=A@coef
        r2=1-np.sum((y-fit)**2)/np.sum((y-y.mean())**2)
        print(f" phase col{k}: R={coef[0]:.4f}  L={coef[1]*1e3:.4f} mH  R2={r2:.6f}  resid rms {np.sqrt(np.mean((y-fit)**2)):.3f} V")
    # spectrum of vn
    w=2*np.pi*50; tt=t[m]
    for h in [1,3,5,7,9]:
        c=2*np.mean(vn[m]*np.exp(-1j*h*w*tt)); print(f"   vn harmonic {h}: {abs(c)/np.sqrt(2):.2f} V rms")
