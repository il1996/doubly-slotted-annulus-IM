import numpy as np, os
REF = "/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/repo/reference/ANSYS_18_5kW"
def load(p): return np.loadtxt(p, skiprows=1)
for label, sub, cur in [("a vide","transitoire/a vide","Winding Plot 4.tab"),("en charge","transitoire/en charge","Winding Plot 4.tab"),("rotor bloque","transitoire/rotor bloqu#U00e9","Winding Plot 3.tab")]:
    d=os.path.join(REF,sub)
    V=load(os.path.join(d,"Winding Plot 3.tab")) if label!="rotor bloque" else None
    I=load(os.path.join(d,cur)); E=load(os.path.join(d,"Winding Plot 2.tab")); PSI=load(os.path.join(d,"Winding Plot 1.tab"))
    t=I[:,0]; m=(t>=1.0)&(t<=2.0)
    print("="*80); print(label)
    Irms=np.sqrt(np.mean(I[m,1:]**2,axis=0)); print(" I rms", Irms, "mean", Irms.mean())
    Erms=np.sqrt(np.mean(E[m,1:]**2,axis=0)); print(" E(2D induced) rms", Erms, "mean", Erms.mean())
    if V is not None:
        Vrms=np.sqrt(np.mean(V[m,1:]**2,axis=0)); print(" V rms", Vrms, "mean", Vrms.mean())
        P=np.mean(np.sum(V[m,1:]*I[m,1:],axis=1)); S=3*Vrms.mean()*Irms.mean()
        print(f" P_elec (mean sum v.i) = {P:.2f} W ; S = {S:.1f} VA ; PF = {P/S:.4f}")
        # power via induced voltage (air-gap+2D leakage EMF times current) = electromagnetic + core? 
        Pe=np.mean(np.sum(E[m,1:]*I[m,1:],axis=1)); print(f" P through 2D induced voltage (sum e.i) = {Pe:.2f} W")
        Pcu=3*0.4457*Irms.mean()**2; print(f" 3 R I^2 with R=0.4457: {Pcu:.1f} W ; P_elec - Pcu = {P-Pcu:.1f}")
    # phasors from last 20 ms (one cycle) via Fourier at 50 Hz over 1.0-2.0 s (50 cycles)
    w=2*np.pi*50; tt=t[m]
    def ph(x): return 2*np.mean(x*np.exp(-1j*w*tt))
    Ia=ph(I[m,1]); Ea=ph(E[m,1]); Pa=ph(PSI[m,1])
    print(f" phasor |Ia|={abs(Ia)/np.sqrt(2):.4f} A rms, |Ea|={abs(Ea)/np.sqrt(2):.3f} V rms, angle(Ea)-angle(Ia) = {np.degrees(np.angle(Ea)-np.angle(Ia)):.2f} deg")
    print(f" psi phasor {abs(Pa):.5f} Wb peak; w*psi = {w*abs(Pa)/np.sqrt(2):.3f} V rms")
    if V is not None:
        Va=ph(V[m,1]); print(f" |Va|={abs(Va)/np.sqrt(2):.3f} V rms; angle(Va)-angle(Ia)={np.degrees(np.angle(Va)-np.angle(Ia)):.2f} deg -> PF_fund={np.cos(np.angle(Va)-np.angle(Ia)):.4f}")
        # Z_total = V/I
        Z=Va/Ia; print(f" Z=V/I = {Z.real:.4f} + j{Z.imag:.4f} ohm ; E/I: {(Ea/Ia).real:.4f} + j{(Ea/Ia).imag:.4f}")
        # check the sign convention: (V - E)/I should be R_s + jX_ext
        Zext=(Va-Ea)/Ia; print(f" (V-E)/I = {Zext.real:.4f} + j{Zext.imag:.4f} ohm  (expected R_s + j X_end)")
        Zext2=(Va+Ea)/Ia; print(f" (V+E)/I = {Zext2.real:.4f} + j{Zext2.imag:.4f} ohm")
