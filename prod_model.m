function m = prod_model( m1, m2 )
% PROD_MODEL Calculates model matrix corresponding to interactions.
%  
% PROD_MODEL(M1,M2) calculates the model matrix corresponding to the
% interaction terms defined by the model matrices M1 and M2.  
% M1 = first model matrix
% M2 = second model matrix  
%
% See also  ADD_MODEL, ASFACTOR.
  
  [ n, p1 ] = size( m1 );
  [ n, p2 ] = size( m2 );

  p = p1*p2;
  m = zeros( n, p );
  
  for( i=1:p1 )
    col1 = (i-1)*p2 + 1;
    col2 = (i-1)*p2 + p2;
    m( :, col1:col2 ) = repmat( m1(:,i), 1, p2 ) .* m2;
  end
  