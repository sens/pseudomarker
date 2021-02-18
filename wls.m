function [b,rss,v] = wls( y, X, A )
%
  [n,p] = size(X);
  w=sqrt(A);
  XX = repmat(w,1,p).*X;
  yy = w.*y;
  
  % for shrinkage option
  % if(p>1)
  %  Xnew = Xnew'*Xnew + ones(p,p)
  % end
  
  [Q, R]=qr(XX,0);

  r = sum( abs( diag(R) ) > 1e-10 ); % rank
  if(r<p)
    warning( 'Model matrix not of full rank.' );
    % R=R(1:r,:);
    % Q=Q(:,1:r);
  end
  
   b = R\(Q'*yy);
  
  if( nargout>1)
    yhat = XX*b;		
    rss = sum((yy-yhat).*(yy-yhat));
  end

  if( nargout>2 )
    v = R'*R;
  end
   
