% This script is for the calculation of time-averaged transition probability versus \omega_0
clc;
clear all;
global w1 w2 kxv Nkx kx
kx=15;% the largest Fourier mode to be retained
kxv=-kx:kx;
Nkx=length(kxv);
Nk=Nkx;
Nh=2*Nkx;
[Imx,Dmx,Cmx]=genelemat(kx);
Cmn=Cmx;Cmp=Cmx.';
%%=================parameters========================
w1=1.0;w2=1.2;
Delta=w2-w1;r=1.0/2.0;
A=0.5;% A strands for A_1
phi1=0;phi2=0;
%%=========================================
sz=[1 0;0 -1];
sx=[0 1;1 0];
sp=[0 1;0 0];
sn=[0 0;1 0];
alpha0=find(kron(kxv==0,[1,0])==1);
% wa=0.5*(w1+w2);
% wxl=[linspace(0.95,0.9544,50) linspace(0.9547,0.9566,300) linspace(0.9569,0.9681,100) linspace(0.9683,0.9735,300) linspace(0.9739,1.0245,200) linspace(1.025,1.0315,200) linspace(1.032,1.0405,100) linspace(1.041,1.0445,300) linspace(1.045,1.05,50)];
% wxl=linspace(0.95,1.05,3001);
% wxl=[linspace(0.25,0.4370,100) linspace(0.4371,0.4378,300) linspace(0.4380,0.6686,100) linspace(0.6688,0.6901,300) linspace(0.691,1.496,500) linspace(1.497,1.522,300) linspace(1.523, 1.732,100) linspace(1.733,1.734,300) linspace(1.735,1.85,100)];
% wxl=[linspace(0.25,0.4521,100) linspace(0.4525,0.4557,300) linspace(0.456,1.673,1000) linspace(1.676,1.692,300) linspace(1.695,1.85,100)];
wxl=[linspace(0.25,0.2780,50) linspace(0.2781,0.2784,500) linspace(0.2790,1.760,1000) linspace(1.767,1.784,300) linspace(1.786,1.85,100)];
% [X,Y]=meshgrid(wxl,Al);
nx=length(wxl);
pav=zeros(1,nx);
%% ==============================================
for ii=1:nx
    w0=wxl(ii);
    delta1t=w0*besselj(0,A/w1*xifun1(w0,A,r))*besselj(0,r*A/w2*xifun2(w0,A,r))-w1;
    A0t=w0*besselj(1,A/w1*xifun1(w0,A,r))*besselj(1,r*A/w2*xifun2(w0,A,r));
    A1t=2*A*(1-xifun1(w0,A,r));
    A2t=2*r*A*(1-xifun2(w0,A,r));
    dphi21=phi2-phi1;
    Hchrwt=kron(Imx,0.5*(delta1t*sz+0.5*A1t*sx))+kron(Dmx,Delta*eye(2))+kron(Cmp,0.5*(0.5*A2t*sn-A0t*sz)*exp(1i*dphi21))+kron(Cmn,0.5*(0.5*A2t*sp-A0t*sz)*exp(-1i*dphi21));
    [u,v]=eig(Hchrwt);
    [xppz,xppp,xppn]=xppfun(u(1:2:end-1,alpha0),u(2:2:end,alpha0));
    pav(ii)=0.5*(1-sum(besselj(kxv,A/w1*xifun1(w0,A,r)).*besselj(-kxv,r*A/w2*xifun2(w0,A,r)).*xppz.*exp(-1i*kxv*dphi21)+exp(-1i*kxv*dphi21).*besselj(1-kxv,A/w1*xifun1(w0,A,r)).*besselj(kxv,r*A/w2*xifun2(w0,A,r)).*xppp+exp(-1i*kxv*dphi21).*besselj(1+kxv,A/w1*xifun1(w0,A,r)).*besselj(-kxv,r*A/w2*xifun2(w0,A,r)).*xppn)^2);
end
figure(1);plot(wxl,pav),hold on