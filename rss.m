function [ss,df] = rss( y, X )
% RSS Residual sum of squares or determinant of RSS for univariate or
%     multivariate normal distributions.
%
% RSS(Y,X) calculates the residual sum of squares from the model
% specified by the model matrix X applied to the data vector Y.
%
% The matrix X should be of full rank.  If not, the program will delete
% columns from the model matrix, so that it is of full rank.  The non-full
% rank case of this program has not been properly tested.

%  [ss,df]=rss4(y,X);
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
  ss = (y-yhat)'*(y-yhat);
  ss = det(ss);

  
