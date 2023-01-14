function [xppz,xppp,xppn]=xppfun(uu,ud)
global kxv Nkx kx
xppz=zeros(1,Nkx);
xppp=zeros(1,Nkx);
for nn=kxv
    xppz(nn+kx+1)=sum(conj(uu(1-nn*(nn<0):end-nn*(nn>0))).*uu(1+nn*(nn>0):end+nn*(nn<0))-conj(ud(1-nn*(nn<0):end-nn*(nn>0))).*ud(1+nn*(nn>0):end+nn*(nn<0)));
    xppp(nn+kx+1)=sum(conj(uu(1-nn*(nn<0):end-nn*(nn>0))).*ud(1+nn*(nn>0):end+nn*(nn<0)));
end
xppn=conj(fliplr(xppp));
end 