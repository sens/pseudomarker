function [lod,bf] = threescan3( y, x, igenoA, igenoB, igenoC, cross, model )
% THREESCANCAN3 Cmpute the approximate posterior distribution of QTL assuming
% there are only three of them on three different chromosomes.
% 
% y = vector of trait values  
% x = matrix of covariates, use [] if none
% igenoA = array of imputed genotypes for first chromosome; the number of
% "pages" in this array determines the accuracy of the "lod" score
% igenoB = array of imputed genotypes for second chromosome 
% igenoC = array of imputed genotypes for third chromosome 
% cross = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS
% model = model works for BACKCROSS ONLY
  
[ n, ntr ] = size(y);  
[ n, p ] = size(x);  
[ n, mA, npages ] = size( igenoA );
[ n, mB, npages ] = size( igenoB );
[ n, mC, npages ] = size( igenoC );

lod = zeros( npages, mA, mB, mC );
%[rawss,df] = rss( y, ones(n,1) );
[rawss,df] = rss( y, [ ones(n,1) x ] );

if( cross == 'bc' )
  switch model
   case 'a+b+c'
    modelmat = ones( n, 4 );    
   case 'a+b+c+ab+bc+ca'
    modelmat = ones( n, 7 );
    case 'a*b*c'
    modelmat = ones( n, 8 );
  end

elseif ( cross == 'f2' )
  modelmat = ones( n, 27 );
  switch model
   case 'a+b+c'
    modelmat = ones( n, 7 );    
   case 'a+b+c+ab+bc+ca'
    modelmat = ones( n, 19 );
   case 'a*b*c'
    modelmat = ones( n, 27 );
  end
end


for( i=1:mA )
  for( j=1:mB )
    for( k=1:mC )
      for( l=1:npages )

	if( cross=='bc' )
	  
	  switch model
	   case 'a+b+c'
	    modelmat( :, 2 ) = igenoA(:,i,l); % A
	    modelmat( :, 3 ) = igenoB(:,j,l); % B
	    modelmat( :, 4 ) = igenoC(:,k,l); % C
	   case 'a+b+c+ab+bc+ca'
	    modelmat( :, 2 ) = igenoA(:,i,l); % A
	    modelmat( :, 3 ) = igenoB(:,j,l); % B
	    modelmat( :, 4 ) = igenoC(:,k,l); % C
	    modelmat( :, 5 ) = modelmat(:,2).*modelmat(:,3); % AB
	    modelmat( :, 6 ) = modelmat(:,2).*modelmat(:,4); % AC
	    modelmat( :, 7 ) = modelmat(:,3).*modelmat(:,4); % BC
	   case 'a*b*c'
	    modelmat( :, 2 ) = igenoA(:,i,l); % A
	    modelmat( :, 3 ) = igenoB(:,j,l); % B
	    modelmat( :, 4 ) = igenoC(:,k,l); % C
	    modelmat( :, 5 ) = modelmat(:,2).*modelmat(:,3); % AB
	    modelmat( :, 6 ) = modelmat(:,2).*modelmat(:,4); % AC
	    modelmat( :, 7 ) = modelmat(:,3).*modelmat(:,4); % BC
	    modelmat( :, 8 ) = modelmat(:,5).*modelmat(:,4); % ABC     
	  end
	  
	elseif( cross=='f2' )

	  modelmat( :, [ 2 3 ] ) = f2model( igenoA(:,i,l) ); %A
	  modelmat( :, [ 4 5 ] ) = f2model( igenoB(:,j,l) ); %B
	  modelmat( :, [ 6 7 ] ) = f2model( igenoC(:,k,l) ); %C   
	  
	  switch model
	   case 'a+b+c'
	    % do nothing
	   case 'a+b+c+ab+bc+ca'
	    % AB, BC, CA
	    modelmat( :, 8:11 ) = prod_model( modelmat( :, [ 2 3 ] ), ...
					      modelmat( :, [ 4 5 ] ) );
	    modelmat( :, 12:15 ) = prod_model( modelmat( :, [ 4 5 ] ), ...
					       modelmat( :, [ 6 7 ] ) );
	    modelmat( :, 16:19 ) = prod_model( modelmat( :, [ 6 7 ] ), ...
					       modelmat( :, [ 2 3 ] ) );
	   case 'a*b*c'
	    % AB, BC, CA, ABC
	    modelmat( :, 8:11 ) = prod_model( modelmat( :, [ 2 3 ] ), ...
					      modelmat( :, [ 4 5 ] ) );
	    modelmat( :, 12:15 ) = prod_model( modelmat( :, [ 4 5 ] ), ...
					       modelmat( :, [ 6 7 ] ) );
	    modelmat( :, 16:19 ) = prod_model( modelmat( :, [ 6 7 ] ), ...
					       modelmat( :, [ 2 3 ] ) );
	    modelmat( :, 20:27 ) = prod_model( modelmat( :, [ 8:11 ] ), ...
					       modelmat( :, [ 6 7  ] ) );
	  end
	  
	  
	end
	
	[ rss1, df ] = rss( y, [ modelmat x ] );
	%lod( l, i, j, k ) = - ((df-ntr+1)/2) * log( rss/rawss );
	lod( l, i, j, k ) = - (n/2) * log( rss1/rawss );
	
      end
    
    end
  end
end


if( size(lod,1) > 1 )
  const = max( lod(:) );
  %lod = mean( exp( lod - const ) );
  lod = lognormalmean( lod - const );
  lod = log( lod ) + const;
end

lod = squeeze( lod );
tmp = lod ( find( lod~=0 ) );
bf = mean( exp( tmp ) );

