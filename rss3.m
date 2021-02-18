function r = rss3( y, z, x )
% RSS3 Calculate type III sum of squares  
%
% RSS3(Y,Z)  
% RSS3(Y,Z,X)    
%
  
  n = size(y,1);
  nterms = size(z,2);
  r = zeros(nterms,1);
  if( nargin==2 )
    x = ones(n,1);
  end
  
  for( i=1:nterms )
    w = z;
    w(:,i)=[];
    r(i) = rss( y, [ w x ] ) - rss( y, [ z x ] );
  end
  