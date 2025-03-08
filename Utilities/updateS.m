
function [A,S] = updateS(X, distX,idx,k,gamma,beta,S,C)

num = size(X,2);
A = zeros(num);
for i=1:num
     di = distX(i,2:k+2)-(4/beta)*S(i,idx(2:k+2));
    sumk=0;
   for j=1:k
    sumk=sumk+di(j);
   end
    for j=1:k
        o=isnan((di(k+1)-di(j))/(k*di(k+1)-sumk));
        if o==0
 	 A(i,idx(i,j+1))=(di(k+1)-di(j))/(k*di(k+1)-sumk);
        else
      A(i,idx(i,j+1))=C(i,idx(i,j+1));
        end
    end
end


