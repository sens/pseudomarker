function [b,rss,v] = bls( y, X, b0, v0 )
% BLS Fit a Bayesian least squares regression model
%
  [n,p] = size(X);
  df = n - p;
  
  [a,b] = size(v0);
  w = sqrt(v0);
  
  XX = zeros(n+p,p);
  XX(1:n,:) = X;
  XX((n+1):(n+p),:) = diag(w);
  
  yy = zeros(n+p,1);
  yy(1:n,:) = y;
  a = w.*b0;
  yy((n+1):(n+p),1) = a';
  
  [Q, R]=qr(XX,0);
  r = sum( abs( diag(R) ) > 1e-10 ); % rank
  if(r<p)
    warning( 'Model matrix not of full rank.' );
    % R=R(1:r,:);
    % Q=Q(:,1:r);
  end
  
  b = R\(Q'*yy);
  
  if( nargout>1 )
    yhat = XX*b;
    rss = (yy-yhat)'*(yy-yhat);
  end

  if( nargout>2 )
    v = R'*R;
  end
  