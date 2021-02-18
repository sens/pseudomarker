function [m,v,rss] = parameters( y, X )
% PARAMETERS Parameter estimates and their variances.
%

  %lastwarn('');
  [n,t] = size(y);  
  [n,p] = size(X);
  df = n - p;			
  [Q, R]=qr(X,0);
  nrow=size(R);
  r = sum( abs( diag(R) ) > 1e-10 );
  if( r<p )
    %b=pinv(R)*(Q'*y)
    %vvv = pinv( R'*R );
    %v = diag(vvv)
    yhat = mean(y);
    rss = (y-yhat)'*(y-yhat);
    m = 0;
    v = 0;
    %lastwarn('');
  else
    b = R\(Q'*y);
    vvv = inv( R'*R );
    v = diag(vvv);

    yhat = X*b;		
    rss = (y-yhat)'*(y-yhat);
    m = b;

    if( t==1 )
      sigma2hat = rss/n;
    else
      sigma2hat = 1 / diag( inv(rss/n) );
    end
    
    %sigma2hat = rss/chi2rnd(n);    
    v = repmat(sigma2hat,1,p).*repmat(v,1,t)';
  end

