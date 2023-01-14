% Script for the CHRW calculation of the resonance positions when the
% resonance peak is narrow: solving d=0
clc;
clear all;
global w1 phi1 phi2 r sz sx Delta Dmx Imx alpha0 Cmn Cmp w2 sn sp kxv Nkx kx
kx=15;% for the first mode
kxv=-kx:kx;
Nkx=length(kxv);
Nk=Nkx;
Nh=2*Nkx;
[Imx,Dmx,Cmx]=genelemat(kx);
Cmn=Cmx;Cmp=Cmx.';
%%=================parameters========================
w0=1.0;w1=1.0;w2=1.2;
Delta=w2-w1;r=0.5;
phi1=0;phi2=0;
%=========================================
sz=[1 0;0 -1];
sx=[0 1;1 0];
sp=[0 1;0 0];
sn=[0 0;1 0];
alpha0=find(kron(kxv==0,[1,0])==1);
rp_guess=dlmread('band4chrwr0p5L.dat');% 
Axl=rp_guess(:,2).';
wres=rp_guess(:,1);
nx=length(Axl);
bslist=zeros(2,nx);
for ii=1:nx
    A=Axl(ii);
    [xs,fval]=fzero(@(x) myfundpavchrw(x,A),wres(ii));
    bslist(:,ii)=[xs;fval];
    disp(fval);
end
figure(1);plot(bslist(1,:),Axl),hold on