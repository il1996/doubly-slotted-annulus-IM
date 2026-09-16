import numpy as np, os
REF = "/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/repo/reference/ANSYS_18_5kW/transitoire/en charge"
def load(p): return np.loadtxt(p, skiprows=1)
V=load(os.path.join(REF,"Winding Plot 3.tab")); I=load(os.path.join(REF,"Winding Plot 4.tab")); E=load(os.path.join(REF,"Winding Plot 2.tab")); PSI=load(os.path.join(REF,"Winding Plot 1.tab"))
t=I[:,0]; dt=t[1]-t[0]
dpsi=np.gradient(PSI[:,1],dt)
m=(t>=1.0)&(t<=2.0)
print("corr(E_A, dpsi_A/dt) =", np.corrcoef(E[m,1],dpsi[m])[0,1], " ratio rms", np.sqrt(np.mean(E[m,1]**2))/np.sqrt(np.mean(dpsi[m]**2)))
# check which relation holds:  V - R I - E = L dI/dt ?
R=0.4457
dI=np.gradient(I[:,1],dt)
for sgn in [+1,-1]:
    resid=V[m,1]-R*I[m,1]-sgn*E[m,1]
    # least squares L
    L=np.sum(resid*dI[m])/np.sum(dI[m]**2); r2=1-np.sum((resid-L*dI[m])**2)/np.sum(resid**2)
    print(f"sign {sgn:+d}: fitted L_ext = {L*1e3:.4f} mH, R2 = {r2:.4f}, resid rms {np.sqrt(np.mean(resid**2)):.2f} V")
# Try V - E = R I + L dI/dt  with fitted R and L
A=np.vstack([I[m,1],dI[m]]).T
for sgn in [+1,-1]:
    y=V[m,1]-sgn*E[m,1]; coef,res,_,_=np.linalg.lstsq(A,y,rcond=None)
    fit=A@coef; r2=1-np.sum((y-fit)**2)/np.sum((y-y.mean())**2)
    print(f"sign {sgn:+d}: V - ({sgn:+d})E = R I + L dI/dt -> R={coef[0]:.4f} ohm, L={coef[1]*1e3:.4f} mH, R2={r2:.5f}")
