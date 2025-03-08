function [E] = updateE(X,alpha,Z,Y1,mu,type)
[m,n]=size(X);
if type
   E = 1./(1+alpha)*(X-X*Z-(Y1/mu));
else 
   E = max(abs(X-X*Z-(Y1/mu))-((alpha/mu)*eye(m,n)),0);
end