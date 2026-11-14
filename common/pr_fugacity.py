import numpy as np
R=8.314462618
# CH4, CO2
Tc=np.array([190.564,304.1282]); Pc=np.array([4.5992e6,7.3773e6]); w=np.array([0.01142,0.22394])
kij=np.array([[0,0.095],[0.095,0]])
def params(T):
    k=0.37464+1.54226*w-0.26992*w**2
    a=0.45724*R**2*Tc**2/Pc*(1+k*(1-np.sqrt(T/Tc)))**2
    b=0.07780*R*Tc/Pc
    return a,b
def lnphi(x,T,P,root):
    a,b=params(T)
    aij=np.sqrt(np.outer(a,a))*(1-kij)
    am=x@aij@x; bm=x@b
    A=am*P/(R*T)**2; B=bm*P/(R*T)
    Z=np.roots([1,-(1-B),A-3*B**2-2*B,-(A*B-B**2-B**3)])
    Z=np.real(Z[np.abs(np.imag(Z))<1e-9]); Z=Z[Z>B]
    Z=Z.max() if root=='v' else Z.min()
    s=2*(aij@x)/am
    lp=b/bm*(Z-1)-np.log(Z-B)-A/(2*np.sqrt(2)*B)*(s-b/bm)*np.log((Z+(1+np.sqrt(2))*B)/(Z+(1-np.sqrt(2))*B))
    return lp,Z
def best(x,T,P):
    out=[]
    for r in 'vl':
        lp,Z=lnphi(x,T,P,r); g=x@(np.log(x)+lp); out.append((g,lp,Z))
    return min(out,key=lambda t:t[0])
def stable(z,T,P):
    _,lpz,_=best(z,T,P); d=np.log(z)+lpz
    for K in [np.array([10,0.1]),np.array([0.1,10])]:
        Y=z*K
        for _ in range(500):
            y=Y/Y.sum(); _,lpy,_=best(y,T,P)
            Yn=np.exp(d-lpy)
            if np.max(abs(Yn-Y))<1e-12: break
            Y=Yn
        if Y.sum()>1+1e-6 and np.max(abs(Y/Y.sum()-z))>1e-4: return False
    return True
z=np.array([0.5,0.5])
print("T(K) P(MPa)  stable  Z    rho(1/nm3)  phi_CH4 phi_CO2  f_CH4(MPa) f_CO2(MPa)")
for T,P in [(275,1),(275,3),(275,5),(275,7),(275,10),(275,20),(275,50),(295,10),(315,10),(295,5),(315,5),(265,10),(295,1),(315,1)]:
    g,lp,Z=best(z,T,P*1e6); rho=P*1e6/(Z*R*T)*6.02214e23/1e27
    ph=np.exp(lp)
    print(f"{T:4d} {P:6.1f}  {str(stable(z,T,P*1e6)):6s} {Z:5.3f}  {rho:6.2f}   {ph[0]:.3f}  {ph[1]:.3f}   {ph[0]*0.5*P:.3f}  {ph[1]*0.5*P:.3f}")
