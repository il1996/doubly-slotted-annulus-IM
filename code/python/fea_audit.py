import numpy as np, glob, os, sys
REF = "/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/repo/reference/ANSYS_18_5kW"
def load(path):
    with open(path, encoding='utf-8', errors='replace') as f:
        hdr = f.readline().strip().split('\t')
    A = np.loadtxt(path, skiprows=1)
    return hdr, A
def stats(name, t, y, windows):
    out=[]
    for (t0,t1) in windows:
        m=(t>=t0-1e-9)&(t<=t1+1e-9)
        yy=y[m]
        out.append((t0,t1,yy.mean(),yy.std(),yy.min(),yy.max(), np.sqrt(np.mean(yy**2))))
    return out
runs = {"a vide":"transitoire/a vide", "en charge":"transitoire/en charge", "rotor bloque":"transitoire/rotor bloqu#U00e9"}
for label, sub in runs.items():
    d=os.path.join(REF,sub)
    print("="*100); print("RUN:",label)
    # speed
    for fn in ["la vitesse en fonction du temps.tab","Speed Plot 1.tab"]:
        p=os.path.join(d,fn)
        if os.path.exists(p):
            h,A=load(p); t=A[:,0]; w=A[:,1]
            print(" speed file:",fn," t range",t[0],t[-1],"n",len(t))
            for (t0,t1) in [(1.0,2.0),(1.5,2.0),(1.9,2.0)]:
                m=(t>=t0)&(t<=t1); print(f"   speed mean {t0}-{t1}: {w[m].mean():.4f} rpm, min {w[m].min():.3f} max {w[m].max():.3f}")
    # torque
    for fn in ["Torque Plot 1.tab","Plot 1.tab"]:
        p=os.path.join(d,fn)
        if os.path.exists(p):
            h,A=load(p); t=A[:,0]
            for j in range(1,A.shape[1]):
                y=A[:,j]
                print(" torque col:",h[j])
                for (t0,t1) in [(1.0,2.0),(1.5,2.0),(1.9,2.0),(1.0,1.5)]:
                    m=(t>=t0)&(t<=t1); print(f"   mean {t0}-{t1}: {y[m].mean():.5f}  std {y[m].std():.4f}  pk-pk {y[m].max()-y[m].min():.4f}")
    # currents
    p=os.path.join(d,"Winding Plot 4.tab")
    if not os.path.exists(p): p=os.path.join(d,"Winding Plot 3.tab")
    h,A=load(p); t=A[:,0]
    if 'Current' in h[1]:
        print(" current file:",os.path.basename(p), h[1:])
        for (t0,t1) in [(1.0,2.0),(1.5,2.0),(1.9,2.0)]:
            m=(t>=t0)&(t<=t1)
            rms=[np.sqrt(np.mean(A[m,j]**2)) for j in range(1,4)]
            print(f"   rms {t0}-{t1}: A {rms[0]:.4f}  C {rms[1]:.4f}  B {rms[2]:.4f}  mean3 {np.mean(rms):.4f}")
    # induced voltage (EMF)
    p=os.path.join(d,"Winding Plot 2.tab"); h,A=load(p); t=A[:,0]
    print(" induced voltage:",h[1:])
    for (t0,t1) in [(1.0,2.0),(1.5,2.0),(1.9,2.0)]:
        m=(t>=t0)&(t<=t1)
        rms=[np.sqrt(np.mean(A[m,j]**2)) for j in range(1,4)]
        print(f"   rms {t0}-{t1}: A {rms[0]:.3f}  C {rms[1]:.3f}  B {rms[2]:.3f}  mean3 {np.mean(rms):.3f}")
    # flux linkage
    p=os.path.join(d,"Winding Plot 1.tab"); h,A=load(p); t=A[:,0]
    for (t0,t1) in [(1.0,2.0),(1.9,2.0)]:
        m=(t>=t0)&(t<=t1)
        pk=[A[m,j].max() for j in range(1,4)]
        print(f"   flux linkage peak {t0}-{t1}: {pk}")
    # losses
    for fn in ["Loss Plot 1.tab","Loss Plot 2.tab","End Connection Plot 3.tab","End Connection Plot 1.tab"]:
        p=os.path.join(d,fn)
        if os.path.exists(p):
            h,A=load(p); t=A[:,0]
            for j in range(1,A.shape[1]):
                y=A[:,j]
                print(" ",h[j], end=": ")
                for (t0,t1) in [(1.0,2.0),(1.5,2.0),(1.9,2.0)]:
                    m=(t>=t0)&(t<=t1); print(f"[{t0}-{t1}] mean {y[m].mean():.5f} rms {np.sqrt(np.mean(y[m]**2)):.5f}", end="  ")
                print()
    # powers
    for fn in ["puissance entre et sortie.tab","puissance #U00e9lectromagn#U00e9tique.tab","Output Variables Plot 3.tab","Output Variables Plot 2.tab","Output Variables Plot 1.tab","puissance de perte m#U00e9canique.tab","Output Variables Plot 4.tab"]:
        p=os.path.join(d,fn)
        if os.path.exists(p):
            h,A=load(p); t=A[:,0]
            for j in range(1,A.shape[1]):
                y=A[:,j]
                print(" ",fn,"|",h[j], end=": ")
                for (t0,t1) in [(1.0,2.0),(1.5,2.0),(1.9,2.0)]:
                    m=(t>=t0)&(t<=t1); print(f"[{t0}-{t1}] mean {y[m].mean():.5f}", end="  ")
                print()
