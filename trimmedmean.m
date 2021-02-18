function [mm,vv] = trimmedmean(x)
% TRIMMEDMEAN Trimmed mean and variance.

%  mm = mean(x); % raw mean
%  vv = var(x); % raw variance
  
  [n,k] = size(x);
  x = sort(x);
  
  t = floor(log(n)/log(2)/2);
  xx = x( (t+1):(n-t), : );
  mm = mean(xx);
  vv = var(xx);
  
%  if(n>=4)
%    xmin = min(x);
%    xmax = max(x);
%    [n,mm,vv]=delobs(xmax,n,mm,vv);
%    [n,mm,vv]=delobs(xmin,n,mm,vv);    
%  end
  
  
%%%%
function [nout,mout,vout]=delobs(x,n,m,v)
  
  if( n==0 )
    error( 'No observation' );
  elseif( n==1 )
    delta = 0;
  else 
    delta = ( n/(n-1) ) * ( m - x ).^2;
  end

  nout = n-1;
  mout = (n*m-x)/nout;
  vout = ( (n-1)*v - delta ) / (nout-1);
