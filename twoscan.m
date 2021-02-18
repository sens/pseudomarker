function [lod,bf] = twoscan( y, x, z, igeno, cross, model )
% TWOSCAN Cmpute the approximate posterior distribution of QTL assuming
% there is only two of them. Handles MULTIVARIATE data.
% 
% y = vector of trait values  
% x = matrix of additive covariates, use [] if none
% z = matrix of interacting covariates, use [] if none  
% igeno = array of imputed genotypes; the number of "pages" in this
%         array determines the accuracy of the "lod" score  
% cross = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS
%
% model = model size works for BAKCCROSS ONLY
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/08/08 19:32:35 $	

[ n, ntr ] = size(y);  
[ n, p ] = size(x);  
[ n, m, npages ] = size( igeno );

covintmodelmat = [];

lod = zeros( npages, m, m );
%[rawss,df] = rss( y, ones(n,1) );
[rawss,df] = rss( y, [ ones(n,1) x ] );

if( cross == 'bc' )
  switch model
   case 'a+b'
    modelmat = ones( n, 3 );    
   otherwise
    modelmat = ones( n, 4 );
  end

elseif ( cross == 'f2' )
  switch model
   case 'a+b'
    modelmat = ones( n, 5 );
   otherwise
    modelmat = ones( n, 9 );
  end
end

for( i=1:m )
  for( j=(i+1):m )

    for( k=1:npages )

    if( cross=='bc' )
      switch model
       case 'a+b'       
	modelmat( :, 2 ) = bcmodel( igeno(:,i,k) );
	modelmat( :, 3 ) = bcmodel( igeno(:,j,k) );
       otherwise
	modelmat( :, 2 ) = bcmodel( igeno(:,i,k) );
	modelmat( :, 3 ) = bcmodel( igeno(:,j,k) );
	modelmat( :, 4 ) = modelmat(:,2).*modelmat(:,3);
      end
      
    elseif( cross=='f2' )
      switch model
       case 'a+b'
	modelmat( :, [ 2 3 ] ) = f2model( igeno(:,i,k) );
	modelmat( :, [ 4 5 ] ) = f2model( igeno(:,j,k) );
       otherwise
	modelmat( :, [ 2 3 ] ) = f2model( igeno(:,i,k) );
	modelmat( :, [ 4 5 ] ) = f2model( igeno(:,j,k) );
	modelmat( :, 6 ) = modelmat(:,2).*modelmat(:,4);  
	modelmat( :, 7 ) = modelmat(:,2).*modelmat(:,5);  
	modelmat( :, 8 ) = modelmat(:,3).*modelmat(:,4);  
	modelmat( :, 9 ) = modelmat(:,3).*modelmat(:,5);  	
      end
    end
 
    if( ~isempty(z) )
      covintmodelmat = prod_model( modelmat( :, 2:end ), z );
    end
    
    [ rss1, df ] = rss( y, [ modelmat x covintmodelmat ] );    
    lod( k, i, j ) = - (n/2) * log( rss1/rawss );    

    end
    
  end
end

if( size(lod,1) > 1 )
  const = max( lod(:) );
  %lod = mean( exp( lod - const ) );
  %lod = lognormalmean( lod - const );
  lod = wtaverage( lod - const ); 
  lod = log( lod ) + const;
end

lod = squeeze( lod );
tmp = lod ( find( lod~=0 ) );
bf = mean( exp( tmp ) );


