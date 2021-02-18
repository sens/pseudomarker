function f = asfactor( x )
% ASFACTOR Converts a *column* vector to a model matrix.
%
% ASFACTOR(X) converts the vector X into a model matrix. 
% X = numeric vector
% F = model matrix like contr.treatment in S.
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
  f = zeros( n, high-1 );
  
  for( i=1:high-1 )
    f(:,i) = [ levelid==i+1 ];
  end
  