function [ss, df] = rss4( y, X )
% RSS Residual sum of squares or determinant of RSS for univariate or
%     multivariate normal distributions.
%
% RSS(Y,X) calculates the residual sum of squares from the model
% specified by the model matrix X applied to the data vector Y.
%
% The matrix X should be of full rank.  If not, the program will delete
% columns from the model matrix, so that it is of full rank.  The non-full
% rank case of this program has not been properly tested.
  
  [n,p] = size(X);
  [n,q] = size(y);
  df = n - p;			
  %[Q, R]=qr(X,0);
  P = [ X y ];
  Z = P'*P;
  [Q,R] = qr(Z,0);

  % which ones have big contribution to rank
  idx = ( abs( diag(R) ) > 1e-10 );
  r = sum( idx ); 
  % if not full rank compute a generalized inverse
  if(r<p+q)
    Z = Z(find(idx>0),find(idx>0));
    [Q,R]=qr(Z);
  end

  t = R\Q';
  % idx = r+1:r+q;
  ss = 1/det(t(r-q+1:r,r-1+1:r));  
  % ss = 1/det(t(p+1:end,p+1:end));
  %ss = t*t;
  %ss = det(ss);

  
