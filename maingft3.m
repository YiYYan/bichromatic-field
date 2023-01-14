% script for the GFT calculation of the time-averaged transition
% probability versus A_1/\omega_1 and \omega_0/\omega_1
clc;
clear all;
%%=================parameters========================
% kx=8;% for the first mode
% ky=8;% for the second mode
% kxv=-kx:kx;kyv=-ky:ky;
% Nkx=length(kxv);Nky=length(kyv);
% Nk=Nkx*Nky;
% Nh=2*Nkx*Nky;
% alpha0=find(kron(kyv==0,kron(kxv==0,[1,0]))==1);
% [Imx,Dmx,Cmx]=genelemat(kx);
% [Imy,Dmy,Cmy]=genelemat(ky);
w1=1.0;w2=1.2;
Delta=0.2;r=1.0;
A1=0.4;A2=0.02;phi1=0;phi2=0;wa=0.5*(w1+w2);
%=========================================
sz=[1 0;0 -1];
sx=[0 1;1 0];
Al=linspace(0.001,1.0,80);
wxl=linspace(0.38,1.82,600);
[X,Y]=meshgrid(wxl,Al);
nx=numel(X);
pav=zeros(1,nx);
for ii=1:nx
%     w1=X(ii)-0.5*Delta;
%     w2=w1+Delta;
    A1=Y(ii);
    A2=r*A1;
    w0=X(ii);
    if A1<=0.3
        kk=8;
    elseif A1<=0.4&&A1>0.3
        kk=10;
    elseif A1>0.4&&A1<=0.45
        kk=10;
    else
        kk=15;
    end
    kx=kk;% for the first mode
    ky=kk;% for the second mode
    kxv=-kx:kx;kyv=-ky:ky;
    Nkx=length(kxv);Nky=length(kyv);
    Nk=Nkx*Nky;
    Nh=2*Nkx*Nky;
    alpha0=find(kron(kyv==0,kron(kxv==0,[1,0]))==1);
    [Imx,Dmx,Cmx]=genelemat(kx);
    [Imy,Dmy,Cmy]=genelemat(ky);
    h10=0.25*A1*exp(1i*phi1)*sx;
    h01=0.25*A2*exp(1i*phi2)*sx;
    Hdiag=kron(Dmx*w1,kron(Imy,eye(2)))+kron(Imx,kron(Dmy*w2,eye(2)))+kron(Imx,kron(Imy,0.5*w0*sz));
    Hoffdiag=kron(Cmx',kron(Imy,h10))+kron(Imx,kron(Cmy',h01));
    Hf=Hdiag+Hoffdiag+Hoffdiag';
    [u,v]=eig(Hf);
    pav(ii)=sum(abs(u(2:2:end,:)).^2*abs((u(alpha0,:)')).^2);
end
Z=reshape(pav,size(X));
contourf(X,Y,Z,50),hold on
% save('fig1_gft_r0p5new.mat','X','Y','Z')