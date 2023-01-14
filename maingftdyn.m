%script for the calculation of the GFT dynamics
clc;
clear all;
tic;
kx=45;% for the first mode
ky=45;% for the second mode, if consider single mode drive, set ky=0; 
kxv=-kx:kx;kyv=-ky:ky;
Nkx=length(kxv);Nky=length(kyv);
Nk=Nkx*Nky;
Nh=2*Nkx*Nky;
[Imx,Dmx,Cmx]=genelemat(kx);
[Imy,Dmy,Cmy]=genelemat(ky);
%%=================parameters========================
w0=1.0;w1=1.00;w2=1.005;
wa=0.5*(w1+w2);Delta=w2-w1;r=1.0;
A1=0.2;A2=r*A1;phi1=0;phi2=0;
%=========================================
sz=[1 0;0 -1];
sx=[0 1;1 0];
h10=0.25*A1*exp(1i*phi1)*sx;
h01=0.25*A2*exp(1i*phi2)*sx;
alpha0=find(kron(kxv==0,kron(kyv==0,[1,0]))==1);
Hdiag=kron(Dmx*w1,kron(Imy,eye(2)))+kron(Imx,kron(Dmy*w2,eye(2)))+kron(Imx,kron(Imy,0.5*w0*sz));
Hoffdiag=kron(Cmx',kron(Imy,h10))+kron(Imx,kron(Cmy',h01));
Hf=Hdiag+Hoffdiag+Hoffdiag';
[u,v]=eig(Hf);
vd=diag(v).';
kxy=(kron(kxv*w1,ones(1,length(kyv)))+kron(ones(1,length(kxv)),w2*kyv)).';
expt=@(t) exp(1i*kxy*t);
tl=linspace(0,2000,4000);
nx=length(tl);
pt=zeros(1,nx);
for ii=1:nx
    t=tl(ii);
    pt(ii)=sum(conj(u(alpha0+1,:)).*exp(-1i*vd*t).*sum(repmat(expt(t),1,Nh).*u(1:2:end-1,:)));
end
plot(tl,abs(pt).^2),hold on
toc;
% dlmwrite('pav_gft_A1_r1_D0p2.dat',[wxl;pav].','precision',10)