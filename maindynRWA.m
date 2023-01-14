% script is for the RWA calculation of the transition probability versus time
% clc;
clear all;
tic;
global w1 w2 kxv Nkx kx Delta
kx=35;% for the first mode
kxv=-kx:kx;
Nkx=length(kxv);
Nk=Nkx;
Nh=2*Nkx;
[Imx,Dmx,Cmx]=genelemat(kx);
Cmn=Cmx;Cmp=Cmx.';
%%=================parameters========================
w0=1.0;w1=1.00;w2=1.005;
Delta=w2-w1;r=1.0;
A1=0.2;A2=r*A1;phi1=0;phi2=0;
%%===============================================
sz=[1 0;0 -1];
sx=[0 1;1 0];
sp=[0 1;0 0];
sn=[0 0;1 0];
alpha0=find(kron(kxv==0,[1,0])==1);
tt=linspace(0.0,2000,5001);
nx=length(tt);
pt=zeros(1,nx);
%% ==============================================
h1=0.25*A2*exp(1i*(phi1-phi2))*sp;
Delta1=w0-w1;
Hdiag=kron(Dmx*Delta,eye(2))+kron(Imx,0.5*Delta1*sz+0.25*A1*sx);
Hoffdiag=kron(Cmx,h1);
Hf=Hdiag+Hoffdiag+Hoffdiag';
[u,v]=eig(Hf);
ut10=utfun(u(:,alpha0).',0);
ut20=utfun(u(:,alpha0+1).',0);
sz11=ut10'*sz*ut10;
sz21=ut20'*sz*ut10;
sz12=ut10'*sz*ut20;
er12=v(alpha0,alpha0)-v(alpha0+1,alpha0+1);
for ii=1:nx
    pt(ii)=0.5-0.5*(utfun(u(:,alpha0).',tt(ii))'*sz*utfun(u(:,alpha0).',tt(ii))*sz11+real(exp(1i*er12*tt(ii))*utfun(u(:,alpha0).',tt(ii))'*sz*utfun(u(:,alpha0+1).',tt(ii))*sz21));
end
toc;
figure(1);plot(tt,pt),hold on