function [myrank,mycond,rel_err,rel_resid] = cond_fun(A,b,bnew,x,delx)
per_noise=0.02;
myrank=rank(A);
mycond=cond(A);
if nargin < 4 
    x=A\b;
    xnew=A\(bnew);
end
delb=bnew-b;
delx=xnew-x
rel_err=norm(delx)/norm(x);
rel_resid=norm(delb)/norm(b);
end