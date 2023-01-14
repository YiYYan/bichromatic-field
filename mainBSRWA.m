% Script for the RWA calculation of the resonance positions when the
% resonance peak is narrow: solving d=0
clc;
clear all;
global w1 phi1 phi2 r sz sx Delta Dmx Cmx Imx alpha0
kx=15;% for the first mode
kxv=-kx:kx;
Nkx=length(kxv);
Nk=Nkx;
Nh=2*Nkx;
[Imx,Dmx,Cmx]=genelemat(kx);
%%=================parameters========================
w0=1.0;w1=1.0;w2=1.2;
Delta=w2-w1;r=1;
phi1=0;phi2=0;
%=========================================
sz=[1 0;0 -1];
sx=[0 1;1 0];
alpha0=find(kron(kxv==0,[1,0])==1);
rp_guess=dlmread('band1chrwr1.dat');% import the initial guess of the resonance positions
wxl=rp_guess(:,2).';
wres=rp_guess(:,1).';
nx=length(wxl);
bslist=zeros(2,nx);
% options=optimset('TolFun',1e-6,'TolX',1e-6);
for ii=1:nx
    A1=wxl(ii);
%    [xs,fval]=fminbnd(@(x) -myfunpavrwa(w0,x),1e-8,0.3,options);
%     [xs,fval]=fminsearch(@(x) -myfunpavrwa(w0,x),0.5);
    [xs,fval]=fzero(@(x) myfundpavrwa(x,A1),wres(ii));
    disp(fval);
    bslist(:,ii)=[xs;fval];
end
figure(1);plot(bslist(1,:),wxl),hold on
