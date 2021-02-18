function [mmm,vvv] = threeestimate3( y, x, igenoA, igenoB, igenoC, cross,...
				     model )
% THREEESTIMATE3 Cmpute the approximate posterior distribution of QTL assuming
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

b = zeros( npages, mA, mB, mC, modelsize );
v = zeros( npages, mA, mB, mC, modelsize );

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
	   case 'a+b+c+ab+ac+bc'
	    modelmat( :, 2 ) = 2*igenoA(:,i,l)-1; % A
	    modelmat( :, 3 ) = 2*igenoB(:,j,l)-1; % B
	    modelmat( :, 4 ) = 2*igenoC(:,k,l)-1; % C
	    modelmat( :, 5 ) = modelmat(:,2).*modelmat(:,3); % AB
	    modelmat( :, 6 ) = modelmat(:,2).*modelmat(:,4); % AC
	    modelmat( :, 7 ) = modelmat(:,3).*modelmat(:,4); % BC
	   otherwise
	    modelmat( :, 2 ) = 2*igenoA(:,i,l)-1; % A
	    modelmat( :, 3 ) = 2*igenoB(:,j,l)-1; % B
	    modelmat( :, 4 ) = 2*igenoC(:,k,l)-1; % C
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


% ------------------------------------
function u = f2( gt )
  u = [ (gt==1) (gt==2) ]
  
