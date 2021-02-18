function [b,rss,v] = ols( y, X )
% OLS Fit a least squares regression model
%  
% B=OLS(Y,X)
% [B,RSS]=OLS(Y,X)
% [B,RSS,V]=OLS(Y,X)  
%
% Y = respose vector
% X = design matrix
%
% B = parameter estimates
% RSS = residual sum of squares
% V = unscaled precision matrix
  
  [n,p] = size(X);
  df = n - p;			
  [Q, R]=qr(X,0);
  r = sum( abs( diag(R) ) > 1e-10 ); % rank
  if(r<p)
    warning( 'Model matrix not of full rank.' );
    % R=R(1:r,:);
    % Q=Q(:,1:r);
  end
  
  b = R\(Q'*y);
  
  if( nargout>1 )
    yhat = X*b;
    rss = (y-yhat)'*(y-yhat);
  end

  if( nargout>2 )
    v = R'*R;
  end
  