function xi=xifun2(w0,A,r)
global w1 w2
xi=w2/(w0+w2)*(1+A^2*w0/(8*(w0+w2)^3)*(r^2+2*(w0+w2)^2/(w0+w1)^2));
end