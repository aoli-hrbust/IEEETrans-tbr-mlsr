function XX = DataNormal(X)
XX = X./(repmat(sqrt(sum(X.^2,1)),size(X,1),1)+10e-10);
