function f = asfactor2( x )
% ASFACTOR2 Converts a *column* vector to a model matrix.
%
% ASFACTOR2(X) converts the vector X into a model matrix. 
% X = numeric vector
% F = model matrix 
%
% See also ADD_MODEL, PROD_MODEL.
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/08/08 15:18:37 $	

  [levels,xxx,levelid] = unique(x);
  
  if(nargin==1)
    high = max(levelid);
  end
  
  n = length( x );
  f = zeros( n, high );
  
  for( i=1:high )
    f(:,i) = [ levelid==i ];
  end
  