% script is for the CHRW calculation of the transition probability versus time
% clc;
clear all;
tic;
global w1 w2 kxv Nkx kx Delta
kx=15;% for the first mode
kxv=-kx:kx;
Nkx=length(kxv);
Nk=Nkx;
Nh=2*Nkx;
[Imx,Dmx,Cmx]=genelemat(kx);
Cmn=Cmx;Cmp=Cmx.';
%%=================parameters========================
w0=1.0;w1=1.00;w2=1.2;
Delta=w2-w1;r=1.0;
A=0.5;phi1=0;phi2=0;
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
delta1t=w0*besselj(0,A/w1*xifun1(w0,A,r))*besselj(0,r*A/w2*xifun2(w0,A,r))-w1;
A0t=w0*besselj(1,A/w1*xifun1(w0,A,r))*besselj(1,r*A/w2*xifun2(w0,A,r));
z1=A/w1*xifun1(w0,A,r);
z2=r*A/w2*xifun2(w0,A,r);
A1t=2*A*(1-xifun1(w0,A,r));
A2t=2*r*A*(1-xifun2(w0,A,r));
dphi21=phi2-phi1;
Hchrwt=kron(Imx,0.5*(delta1t*sz+0.5*A1t*sx))+kron(Dmx,Delta*eye(2))+kron(Cmp,0.5*(0.5*A2t*sn-A0t*sz)*exp(1i*dphi21))+kron(Cmn,0.5*(0.5*A2t*sp-A0t*sz)*exp(-1i*dphi21));
[u,v]=eig(Hchrwt);
fz=@(t) cos(z1*sin(w1*t+phi1)+z2*sin(w2*t+phi2));
fp=@(t) -1i*sin(z1*sin(w1*t+phi1)+z2*sin(w2*t+phi2))*exp(1i*(w1*t+phi1));
ut10=utfun(u(:,alpha0).',0);
ut20=utfun(u(:,alpha0+1).',0);
sz11=ut10'*sz*ut10;
sz21=ut20'*sz*ut10;
sz12=ut10'*sz*ut20;
er12=v(alpha0,alpha0)-v(alpha0+1,alpha0+1);
for ii=1:nx
    pt(ii)=0.5-0.5*(fz(tt(ii))*(utfun(u(:,alpha0).',tt(ii))'*sz*utfun(u(:,alpha0).',tt(ii))*sz11+real(exp(1i*er12*tt(ii))*utfun(u(:,alpha0).',tt(ii))'*sz*utfun(u(:,alpha0+1).',tt(ii))*sz21))...
        +real(fp(tt(ii))*(2*utfun(u(:,alpha0).',tt(ii))'*sp*utfun(u(:,alpha0).',tt(ii))*sz11+exp(1i*er12*tt(ii))*utfun(u(:,alpha0).',tt(ii))'*sp*utfun(u(:,alpha0+1).',tt(ii))*sz21+exp(-1i*er12*tt(ii))*utfun(u(:,alpha0+1).',tt(ii))'*sp*utfun(u(:,alpha0).',tt(ii))*sz12)));
end
toc;
figure(1);plot(tt,pt),hold on