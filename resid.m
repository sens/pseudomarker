function r = resid( y, X )
% RESID Calculate residuals from regression
%
% R=RESID(Y,X)
% R is the residual from regressing Y on the columns of the matrix X.

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/07/26 14:44:18 $	

  [n,p] = size(X);
  df = n - p;			
  [Q, R]=qr(X,0);
  
  % which ones have big contribution to rank
  idx = ( abs( diag(R) ) > 1e-10 );
  r = sum( idx ); 
  % if not full rank compute a generalized inverse
  if(r<p)
    X = X(:,find(idx>0));
    [Q, R]=qr(X,0);
  end
  
  b = R\(Q'*y);
  yhat = X*b;		
  r = (y-yhat);

  
