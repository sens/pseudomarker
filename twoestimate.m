function [mmm,vvv] = twoestimate( y, x, igeno, cross, model )
% TWOESTIMATE Cmpute the approximate posterior distribution of genetic model
% there is only two of them. 
% 
% y = vector of trait values  
% x = matrix of covariates, use [] if none
% igeno = array of imputed genotypes; the number of "pages" in this
%         array determines the accuracy of the "lod" score  
% cross = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS
%
% model = model size works for BAKCCROSS ONLY
  
[ n, ntr ] = size(y);  
[ n, p ] = size(x);  
[ n, m, npages ] = size( igeno );

lod = zeros( npages, m, m );

[rawss,df] = rss( y, ones(n,1) );

if( cross == 'bc' )
  switch model
   case 'a+b'
    modelmat = ones( n, 3 );    
    modelsize = 3;
   otherwise
    modelmat = ones( n, 4 );
    modelsize = 4;
  end

elseif ( cross == 'f2' )
  switch model
   case 'a+b'
    modelsize = 5;
    modelmat = ones( n, 5 );
   otherwise
    modelsize = 9;
    modelmat = ones( n, 9 );
  end
end

b = zeros( npages, m, m, modelsize );
v = zeros( npages, m, m, modelsize );

for( i=1:m )
  for( j=(i+1):m )

    for( k=1:npages )

    if( cross=='bc' )
      switch model
       case 'a+b'       
	%modelmat( :, 2 ) = igeno(:,i,k);	
	%modelmat( :, 3 ) = igeno(:,j,k);
	modelmat( :, 2 ) = 2*igeno(:,i,k)-1;
	modelmat( :, 3 ) = 2*igeno(:,j,k)-1;	
       otherwise
	%modelmat( :, 2 ) = igeno(:,i,k);
	%modelmat( :, 3 ) = igeno(:,j,k);
	modelmat( :, 2 ) = 2*igeno(:,i,k)-1;
	modelmat( :, 3 ) = 2*igeno(:,j,k)-1;
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
    
    [ b(k,i,j,:), v(k,i,j,:), rss ] = parameters( y, [ modelmat x ] );
    %lod( k, i, j ) = - ((df-ntr+1)/2) * log( rss/rawss );
    lod( k, i, j ) = - (n/2) * log( rss/rawss );    

    end
    
  end
end

const = max( lod(:) );
prob = exp( lod - const );
total = sum(prob(:));
prob = prob/total;

nrep = length( prob(:) );
[mmm,vvv] = unconditional( reshape( b, nrep, modelsize), ...
    reshape( v, nrep, modelsize), prob(:) );


