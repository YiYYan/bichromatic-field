function ut=utfun(ba,t)
global kxv Delta
ut=sum(repmat(exp(1i*kxv*Delta*t),2,1).*[ba(1:2:end-1);ba(2:2:end)],2);
end