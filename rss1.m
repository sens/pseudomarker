function r = rss1( y, z, x )
% RSS1 Calculate type I sum of squares  
%
% RSS1(Y,Z)  
% RSS1(Y,Z,X)    
%
  
  n = size(y,1);
  nterms = size(z,2);
  r = zeros(nterms,1);
  if( nargin==2 )
    x = ones(n,1);
  end

  w0 = [];
  w1 = [];
  for( i=1:nterms )
    w1(:,i)=z(:,i);
    r(i) =  rss( y, [ w0 x ] ) - rss( y, [ w1 x ] ) ;
    w0=w1;
  end
  