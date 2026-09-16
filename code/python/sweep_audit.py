import numpy as np, os
REF = "/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/repo/reference/ANSYS_18_5kW/carat#U00e9ristique en fonction glissement"
T = np.loadtxt(os.path.join(REF,"Torque Plot 2.tab"), skiprows=1)
I = np.loadtxt(os.path.join(REF,"Winding Plot 1.tab"), skiprows=1)
s=T[:,0]; tq=T[:,1]; cur=I[:,1]
print("n points", len(s), "s range", s.min(), s.max(), "step", np.unique(np.round(np.diff(s),5))[:5])
# local roughness: relative deviation from 3-point moving average / or second difference
def rough(y):
    r=np.zeros_like(y); 
    for i in range(1,len(y)-1):
        loc=0.5*(y[i-1]+y[i+1]); r[i]=abs(y[i]-loc)/abs(loc)
    return r
r=rough(tq)
print("torque: max value %.3f at s=%.3f"%(tq.max(), s[np.argmax(tq)]))
for lo,hi in [(0,0.05),(0.05,0.13),(0.13,0.2),(0.2,0.5),(0.5,0.9),(0.9,1.0)]:
    m=(s>lo)&(s<=hi); 
    print(f" band {lo}-{hi}: n={m.sum()}  roughness mean {100*r[m].mean():.2f}%  max {100*r[m].max():.2f}%  torque {tq[m].min():.1f}..{tq[m].max():.1f}")
# print table around rated and breakdown
for i in range(len(s)):
    if s[i]<=0.2 or s[i]>=0.85:
        print(f"  s={s[i]:.3f}  T={tq[i]:8.3f}  I={cur[i]:8.3f}  rough={100*r[i]:.2f}%")
# value at s=0.0188 by interpolation
print("interp T at 0.0188:", np.interp(0.0188,s,tq), " I:", np.interp(0.0188,s,cur))
print("T at s=0.02:", tq[s==0.02], "I:", cur[s==0.02])
# std about local mean in s>=0.9
m=s>=0.9
p=np.polyfit(s[m],tq[m],1); res=tq[m]-np.polyval(p,s[m]); print("s>=0.9: std of residual about linear fit / mean = %.1f%%"%(100*res.std()/tq[m].mean()))
m=(s>0.13)&(s<=0.5)
p=np.polyfit(s[m],tq[m],2); res=tq[m]-np.polyval(p,s[m]); print("0.13<s<=0.5: std residual about quadratic / mean = %.1f%%"%(100*res.std()/tq[m].mean()))
