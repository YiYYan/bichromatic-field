% \bar{P} function for rwa
function pav=myfunpavrwa(w0,A)
global w1 phi1 phi2 r sz sx Delta Dmx Cmx Imx alpha0
    A1=A;A2=r*A1;
    h1=0.25*A2*exp(1i*(phi1-phi2))*[0 1;0 0];
    Delta1=w0-w1;
    Hdiag=kron(Dmx*Delta,eye(2))+kron(Imx,0.5*Delta1*sz+0.25*A1*sx);
    Hoffdiag=kron(Cmx,h1);
    Hf=Hdiag+Hoffdiag+Hoffdiag';
    [u,v]=eig(Hf);
    pav=0.5*(1-sum(abs(u(1:2:end-1,alpha0)).^2-abs(u(2:2:end,alpha0)).^2)^2);%sum(abs(u(1:2:end-1,alpha0)).^2-abs(u(2:2:end,alpha0)).^2);%
end