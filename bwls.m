function [b,rss,v] = bwls( y, X, A, b0, v0 )
% BWLS Fit a Bayesian weighted least squares regression model
%
  [n,p] = size(X);
  df = n - p;
  
  [a,b] = size(v0);
  w1 = sqrt(A);
  w2 = sqrt(v0);
  
  XX = zeros(n+p,p);
  XX(1:n,:) = repmat(w1,1,p).*X;
  XX((n+1):(n+p),:) = diag(w2);
  
  yy = zeros(n+p,1);
  yy(1:n,:) = w1.*y;
  tmp = w2(:).*b0(:);
  yy((n+1):(n+p),1) = tmp;
  
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
  