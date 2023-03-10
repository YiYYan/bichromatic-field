% Script for the RWA calculation of the resonance positions when the
% resonance peak is broad: searching minimum -0.5 of -\bar{P} at the resonance
% band
clc;
clear all;
global w1 phi1 phi2 r sz sx Delta Dmx Cmx Imx alpha0
kx=15;% for the first mode
ky=0;% for the second mode, if consider single mode drive, set ky=0; 
kxv=-kx:kx;
Nkx=length(kxv);
Nk=Nkx;
Nh=2*Nkx;
[Imx,Dmx,Cmx]=genelemat(kx);
%%=================在此处修改参数========================
w0=1.0;w1=1.0;w2=1.2;
Delta=w2-w1;r=1.0;
phi1=0;phi2=0;
%=========================================
sz=[1 0;0 -1];
sx=[0 1;1 0];
alpha0=find(kron(kxv==0,[1,0])==1);
wxl=linspace(1.0,1.2,101);
nx=length(wxl);
bslist=zeros(2,nx);
for ii=1:nx
    w0=wxl(ii);
    [xs,fval]=fminsearch(@(x) -myfunpavrwa(w0,x),0.2);
    bslist(:,ii)=[xs;fval];
    disp(fval);
end
figure(1);plot(wxl,bslist(1,:)),hold on