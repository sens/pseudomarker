function [dev,b,varargout] = logitfit( y, X, m )
% LOGITFIT Fit logistic regression model
%  
  [n,p] = size(X);
  df = n - p;			
  
  if( nargin==2 )
    m = ones( length(y), 1 );
  end
  
  y = repmat( 1/(2*n), n, 1 ) + y;
  m = repmat( 1/n, n, 1 ) + m;
  
  mmm = sum(y.*m)/sum(m);
  bprior = [ invlogit(mmm) zeros(1,p-1)]';
  
  bold = bprior;
  change=1;
  i=0;

  while( change > 1e-10 )

    %i=i+1;
    if( change==1)
      eta = (y+0.5)./(m+1);
    else
      eta = X*bold;
    end
    
    prob = invlogit(eta);
    z = eta + (y-m.*prob)./(m.*(prob.*(1-prob)));
    w = m.*prob.*(1-prob);
    bnew=wls(z,X,w);
    change=sum(abs(bold-bnew));
    bold = bnew;
  end
  
  dev = sum( y.*log(prob) + (1-y).*log(1-prob) );
  b = bnew;
  %b'

  varargout = cell(1,nargout-2);
  
  if(nargout==3)
    nargout;
    varargout(1) = {m.*prob};
  end
  

%%%
function p = invlogit(x)
  p = exp(x)./(1+exp(x));

