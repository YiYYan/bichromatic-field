function pav=myfunpavgft(w0,A)
global r phi1 phi2 sx sz w1 w2 Imx Imy Dmx Dmy Cmx Cmy alpha0
A1=A;A2=r*A;
h10=0.25*A1*exp(1i*phi1)*sx;
h01=0.25*A2*exp(1i*phi2)*sx;
Hdiag=kron(Dmx*w1,kron(Imy,eye(2)))+kron(Imx,kron(Dmy*w2,eye(2)))+kron(Imx,kron(Imy,0.5*w0*sz));
Hoffdiag=kron(Cmx',kron(Imy,h10))+kron(Imx,kron(Cmy',h01));
Hf=Hdiag+Hoffdiag+Hoffdiag';
[u,v]=eig(Hf);
pav=sum(abs(u(2:2:end,:)).^2*abs((u(alpha0,:)')).^2);
end