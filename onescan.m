function [lod,bf] = onescan( y, x, z, igeno, cross )
% ONESCAN Cmpute the approximate posterior distribution of QTL assuming
% there is only one of them. Works for MULTIVARIATE phenotypes.
% 
% y = vector of trait values  
% x = matrix of additive covariates, use [] if none
% z = matrix of interacting covariates, use [] if none  
% igeno = array of imputed genotypes; the number of "pages" in this
%         array determines the accuracy of the "lod" score  
% cross = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2001/08/08 19:29:38 $	
  
[ n, ntr ] = size(y);  
[ n, p ] = size(x);
[ n, m, npages ] = size( igeno );

lod = zeros( npages, m );

%[rawss,df] = rss( y, ones(n,1) );
[rawss,df] = rss( y, [ ones(n,1) x ] );

if( cross == 'bc' )
  modelmat = ones( n, 2 );
elseif ( cross == 'f2' )
  modelmat = ones( n, 3 );
end

covintmodelmat = [];

for( i=1:m )
  for( j=1:npages )

    if( cross=='bc' )
      modelmat( :, 2 ) = bcmodel( igeno(:,i,j) );
    elseif( cross=='f2' )
      modelmat( :, [ 2 3 ] ) = f2model( igeno(:,i,j) );
    end
    
    if( ~isempty(z) )
      covintmodelmat = prod_model( modelmat( :, 2:end ), z );
    end
    [ rss1, df ] = rss( y, [ modelmat x covintmodelmat ] );    
    lod( j, i ) = - (n/2) * log( rss1/rawss );    
  end 
end

if( size(lod,1) > 1 )
  const = max( lod(:) );
  %lod = mean( exp( lod - const ) );
  %lod = lognormalmean( lod - const );
  lod = wtaverage( lod - const );  
  lod = log( lod ) + const;
  lod = squeeze( lod );
end

tmp = lod ( find( lod~=0 ) );
bf = mean( exp( tmp ) );
