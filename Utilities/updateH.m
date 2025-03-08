function [Z] = updateT(S_tensor,sX,mu,viewnum)
    s = S_tensor(:);
    [z, objV] = wshrinkObj(s,1/mu,sX,0,3);
    Z_tensor = reshape(z, sX);
    for i=1:viewnum
    Z{i}=Z_tensor(:,:,i);
    end