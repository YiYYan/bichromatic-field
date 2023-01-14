function pav=myfunpavchrw(w0,A)
global w1 w2 phi1 phi2 r sz sx Delta Dmx Cmn Cmp Imx alpha0 sn sp kxv
delta1t=w0*besselj(0,A/w1*xifun1(w0,A,r))*besselj(0,r*A/w2*xifun2(w0,A,r))-w1;
A0t=w0*besselj(1,A/w1*xifun1(w0,A,r))*besselj(1,r*A/w2*xifun2(w0,A,r));
A1t=2*A*(1-xifun1(w0,A,r));
A2t=2*r*A*(1-xifun2(w0,A,r));
dphi21=phi2-phi1;
Hchrwt=kron(Imx,0.5*(delta1t*sz+0.5*A1t*sx))+kron(Dmx,Delta*eye(2))+kron(Cmp,0.5*(0.5*A2t*sn-A0t*sz)*exp(1i*dphi21))+kron(Cmn,0.5*(0.5*A2t*sp-A0t*sz)*exp(-1i*dphi21));
[u,v]=eig(Hchrwt);
[xppz,xppp]=xppfun(u(1:2:end-1,alpha0),u(2:2:end,alpha0));
pav=0.5*(1-sum(real(besselj(kxv,A/w1*xifun1(w0,A,r)).*besselj(-kxv,r*A/w2*xifun2(w0,A,r)).*xppz.*exp(-1i*kxv*dphi21))+2*real(exp(-1i*kxv*dphi21).*besselj(1-kxv,A/w1*xifun1(w0,A,r)).*besselj(kxv,r*A/w2*xifun2(w0,A,r)).*xppp))^2);
end