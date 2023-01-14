function dy=fundydt(t,y)
global A1 A2 w1 w2 w0 kappa
ft=0.5*(A1*cos(w1*t)+A2*cos(w2*t));
Mat=[1i*w0-0.5*kappa 0 -1i*ft;
    0 -1i*w0-0.5*kappa 1i*ft;
    -2*1i*ft 2*1i*ft -kappa];
Bv=[0;0;-kappa];
dy=Mat*y+Bv;
end