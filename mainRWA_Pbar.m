% script for the calculation of the time-averaged transition probability
% versus \omega_0
clc;
clear all;
kx=20;% for the first mode
kxv=-kx:kx;
Nkx=length(kxv);
Nk=Nkx;
Nh=2*Nkx;
[Imx,Dmx,Cmx]=genelemat(kx);
%%=================parameters========================
w0=1.0;w1=1.0;w2=1.2;
Delta=w2-w1;r=1.0;
A1=0.5;A2=r*A1;phi1=0;phi2=0;wa=0.5*(w1+w2);
%=========================================
sz=[1 0;0 -1];
sx=[0 1;1 0];
h1=0.25*A2*exp(1i*(phi1-phi2))*[0 1;0 0];
alpha0=find(kron(kxv==0,[1,0])==1);
% beta0=find(kron(kyv==0,kron(kxv==0,[0,1]))==1);
wxl=linspace(0.38,1.82,3001);
% wxl=[linspace(0.95,0.958,100) linspace(0.9585,0.9604,300) linspace(0.9606,0.9712,100) linspace(0.9716,0.9766,300) linspace(0.9769,1.027,200) linspace(1.0271,1.033,200) linspace(1.034,1.0435,100) linspace(1.044,1.046,300) linspace(1.0462,1.05,100)];
% wxl=[linspace(0.25,0.4635,100) linspace(0.4636,0.4652,300) linspace(0.4655,0.688,200) linspace(0.6884,0.7191,300) linspace(0.720,1.504,500) linspace(1.505,1.536,300) linspace(1.537,1.745,100) linspace(1.746,1.748,400) linspace(1.749,1.88,200)]; 
% wxl=[linspace(0.25,0.4051,100) linspace(0.4052,0.4132,300) linspace(0.4135,1.540,1000) linspace(1.546,1.584,200) linspace(1.588,1.82,100) linspace(1.824,1.828,300) linspace(1.830,1.85,100)];
% wxl=[linspace(-0.675,-0.61,200) linspace(-0.60992,-0.60984,300) linspace(-0.6094,-0.4423,100) linspace(-0.4421,-0.4313,300) linspace(-0.431,0.2517,500) linspace(0.2518,0.293,300) linspace(0.294,0.4495,100) linspace(0.4496,0.4563,300) linspace(0.4566,0.618,100) linspace(0.61888,0.61896,300) linspace(0.6190,0.675,100)];
% wxl=[linspace(-0.675,-0.4948,200) linspace(-0.4947,-0.4944,300) linspace(-0.4943,-0.3356,100) linspace(-0.3355,-0.3193,300) linspace(-0.319,0.3334,500) linspace(0.3335,0.3456,300) linspace(0.3458,0.4995,100) linspace(0.4995,0.49975,300) linspace(0.5,0.675,200)];
% wxl=[linspace(-0.675,-0.3703,200) linspace(-0.37018,-0.3700,300) linspace(-0.3699,-0.2276,100) linspace(-0.2275,-0.2039,300) linspace(-0.2038,0.2112,500) linspace(0.2113,0.2241,300) linspace(0.2243,0.3709,100) linspace(0.37115,0.3712,300) linspace(0.3714,0.675,200)];
% wxl=[linspace(-0.675,-0.5793,200) linspace(-0.5792,-0.5742,300) linspace(-0.5741,-.4105,100) linspace(-0.4104,-0.3749,300) linspace(-0.375,0.375,500) linspace(0.3749,0.4104,300) linspace(0.4105,0.5741,100) linspace(0.5742,0.5792,300) linspace(0.5793,0.675,200)];
% wxl=[linspace(-0.675,-0.3682,200) linspace(-0.3681,-0.3679,300) linspace(-0.3680,-0.2247,100) linspace(-0.2248,-0.2018,300) linspace(-0.2017,0.2011,500) linspace(0.2018,0.2248,300) linspace(0.225,0.3677,100) linspace(0.3679,0.3681,300) linspace(0.3683,0.675,200)];
% wxl=[linspace(0.5,0.6287,100) linspace(0.62875,0.62885,300) linspace(0.6289, 0.7753,100) linspace(0.7754,0.7894,200) linspace(0.79,1.208,300) linspace(1.209,1.221,300) linspace(1.222,1.37,100) linspace(1.3701,1.3702,400) linspace(1.3703,1.5,100)];  
% wxl=[linspace(0.5,0.50009,10) linspace(0.5001,0.50041,300) linspace(0.5006,0.656,100) linspace(0.6561,0.6635,300) linspace(0.6637,1.322,500) linspace(1.3222,1.333,300) linspace(1.3333,1.4941,100) linspace(1.4942,1.495,300) linspace(1.4955,1.5,100)];
% wxl=[linspace(0.5,0.547,100) linspace(0.5471,0.551,300) linspace(0.552,1.432,1000) linspace(1.4325,1.44,300) linspace(1.4402,1.5,100)];
% wxl=[linspace(0.5,0.5194,20) linspace(0.5195,0.5208,300) linspace(0.521,0.6833,50) linspace(0.6834,0.7117,200) linspace(0.7118,1.291,400) linspace(1.292,1.329,200) linspace(1.330,1.478,50) linspace(1.4781,1.4805,300) linspace(1.481,1.5,50)];
% wxl=[linspace(0.5,0.631,20) linspace(0.6315,0.6323,400) linspace(0.6326,0.778,30) linspace(0.779,0.8,200) linspace(0.802,1.2,200) linspace(1.202,1.22,200) linspace(1.222,1.367,50) linspace(1.3678,1.3682,400) linspace(1.370,1.5,20);];
nx=length(wxl);
pav=zeros(1,nx);
wal=zeros(1,nx);
for ii=1:nx
    w0=wxl(ii);
    Delta1=w0-w1;
    Hdiag=kron(Dmx*Delta,eye(2))+kron(Imx,0.5*Delta1*sz+0.25*A1*sx);
    Hoffdiag=kron(Cmx,h1);
    Hf=Hdiag+Hoffdiag+Hoffdiag';
    [u,v]=eig(Hf);
    pav(ii)=sum(abs(u(2:2:end,:)).^2*abs((u(alpha0,:)')).^2);
end
figure(1);plot(wxl,pav),hold on
% dlmwrite('pav_rwa_A1_r0p5_D0p2new.dat',[wxl;pav].','precision',10)