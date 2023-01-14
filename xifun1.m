function xi=xifun1(w0,A,r)
global w1 w2
xi=w1/(w0+w1)*(1+A^2*w0/(8*(w0+w1)^3)*(1+2*r^2*(w0+w1)^2/(w0+w2)^2));
end