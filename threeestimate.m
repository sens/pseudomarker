function [mmm,vvv] = threeestimate( y, x, igeno, cross, model )
% THREEESTIMATE Cmpute the approximate posterior distribution of QTL assuming
% there are only three of them on the same chromosome.
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
[ n, m, npages ] = size( igeno );

lod = zeros( npages, m, m, m );
[rawss,df] = rss( y, ones(n,1) );

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

modelsize = size(modelmat,2);

b = zeros( npages, m, m, m, modelsize );
v = zeros( npages, m, m, m, modelsize );

for( i=1:m )
  for( j=(i+1):m )
    for( k=(j+1):m )
      for( l=1:npages )

	if( cross=='bc' )
	  
	  switch model
	   case 'a+b+c'
	    modelmat( :, 2 ) = 2*igeno(:,i,l)-1; % A
	    modelmat( :, 3 ) = 2*igeno(:,j,l)-1; % B
	    modelmat( :, 4 ) = 2*igeno(:,k,l)-1; % C
	   case 'a+b+c+ab+bc+ca'
	    modelmat( :, 2 ) = 2*igeno(:,i,l)-1; % A
	    modelmat( :, 3 ) = 2*igeno(:,j,l)-1; % B
	    modelmat( :, 4 ) = 2*igeno(:,k,l)-1; % C
	    modelmat( :, 5 ) = modelmat(:,2).*modelmat(:,3); % AB
	    modelmat( :, 6 ) = modelmat(:,2).*modelmat(:,4); % AC
	    modelmat( :, 7 ) = modelmat(:,3).*modelmat(:,4); % BC
	   case 'a*b*c'
	    modelmat( :, 2 ) = 2*igeno(:,i,l)-1; % A
	    modelmat( :, 3 ) = 2*igeno(:,j,l)-1; % B
	    modelmat( :, 4 ) = 2*igeno(:,k,l)-1; % C
	    modelmat( :, 5 ) = modelmat(:,2).*modelmat(:,3); % AB
	    modelmat( :, 6 ) = modelmat(:,2).*modelmat(:,4); % AC
	    modelmat( :, 7 ) = modelmat(:,3).*modelmat(:,4); % BC
	    modelmat( :, 8 ) = modelmat(:,5).*modelmat(:,4); % ABC     
	  end
	  
	elseif( cross=='f2' )
	  

	  modelmat( :, [ 2 3 ] ) = f2model( igeno(:,i,l) ); %A
	  modelmat( :, [ 4 5 ] ) = f2model( igeno(:,j,l) ); %B
	  modelmat( :, [ 6 7 ] ) = f2model( igeno(:,k,l) ); %C   
	  
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
	
	[ b(l,i,j,k,:) v(l,i,j,k,:), rss ] = parameters( y, [ modelmat x ] );
	%lod( l, i, j, k ) = - ((df-ntr+1)/2) * log( rss/rawss );
	lod( l, i, j, k ) = - (n/2) * log( rss/rawss );
	
      end
    
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




