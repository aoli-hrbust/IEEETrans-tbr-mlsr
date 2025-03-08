function Z = updateZ(X,beta,mu,Y2,Y1,L,E,J)
    n=size(X,2);

     
    lyB =beta*L;

    lyA = (X'*X+10e-10*eye(n))\((eye(n)+(X)'*(X)));

     lyC = -(X'*X+10e-10*eye(n))\((J+(Y2/mu)+(X'*X)-X'*E-X'*(Y1/mu)));
    Z = lyap(lyA,lyB,lyC);
  
