function [b,dev,v,fitted] = blogitfit( y, X, m, b0, v0 )
% BLOGITFIT Fit logistic regression model the Bayesian way
% 
% [B,DEV,V,FITTED]=BLOGITFIT(Y,X,M,B0,V0)
%     Y = 0-1 response added
%     X = covariate
%     M = index of binomial distribution
%     B0 = prior mean
%     V0 = prior unscaled variance
%  
  [n,p] = size(X);
  df = n - p;			
  
  if( nargin==2 )
    m = ones( length(y), 1 );
  end
  
%  I just changed this  
%  y = repmat( 1/(2*n), n, 1 ) + y;
%  m = repmat( 1/n, n, 1 ) + m;
  
  %mmm = sum(y.*m)/sum(m);
  %bprior = [ invlogit(mmm) zeros(1,p-1)]';
  bold = b0;
  
  change=1;
  i=0;

  while( change > 1e-6 )

    %i=i+1;
    if( change==1)
      eta = (y+0.5)./(m+1);
    else
      eta = X*bold;
    end
    
    prob = invlogit(eta);
    z = eta + (y-m.*prob)./(m.*(prob.*(1-prob)));
    w = m.*prob.*(1-prob);
    [bnew,rss,v]=bwls(z,X,w,b0,v0);
    change=sum(abs(bold-bnew))/p;
    bold = bnew;  
  end

  b = bnew;
  tmp = 1./sqrt(v0);
  tmp = log( normpdf( b, b0, tmp' ) );
  % FOLLOWING LINE CHANGED FROM ORIGINAL
  %  dev = sum( y.*log(prob) + (1-y).*log(1-prob) ) + sum(tmp);
  dev = sum( y.*log(prob) + (1-y).*log(1-prob) );

  
  varargout = cell(1,nargout-2);
  
  if(nargout>3)
    fitted = {m.*prob};
  end
  

%%%
function p = invlogit(x)
  p = exp(x)./(1+exp(x));