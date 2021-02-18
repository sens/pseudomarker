function [mmm,vvv,varargout] = oneestimate( y, x, igeno, cross )
% ONEESTIMATE Cmpute the approximate posterior distribution of genetic model
% there is only one of them.
% 
% [M,V]=ONEESTIMATE(Y,X,IGENO,CROSS)
% [M,V,LOD,BF]=ONEESTIMATE(Y,X,IGENO,CROSS)  
% 
% The first form just gives the estimates of the effects, M, and the posterior
% variances, V.  The second form returns also returns the BLOD
% corresponding to the pseudomarker locations in IGENO in the vector LOD
% and the unnormalized Bayes factor in BF.
%  
% Y = vector of trait values  
% X = matrix of covariates, use [] if none
% IGENO = array of imputed genotypes; the number of "pages" in this
%         array determines the accuracy of the "lod" score  
% CROSS = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS
%
  
[ n, ntr ] = size(y);  
[ n, p ] = size(x);
[ n, m, npages ] = size( igeno );

lod = zeros( npages, m );
if( cross=='bc' )
  modelsize=2;
  b = zeros( npages, m, 2 );
  v = zeros( npages, m, 2 );
elseif( cross=='f2' )
  modelsize=3;
  b = zeros( npages, m, 3 );
  v = zeros( npages, m, 3 );
end

rawss = rss( y, ones(n,1) );

if( cross == 'bc' )
  modelmat = ones( n, 2 );
elseif ( cross == 'f2' )
  modelmat = ones( n, 3 );
end

for( i=1:m )
  for( j=1:npages )

    if( cross=='bc' )
      modelmat( :, 2 ) = bcmodel( igeno(:,i,j) );      
    elseif( cross=='f2' )
      % for usual contrasts
      %      modelmat( :, [ 2 3 ] ) = f2model( igeno(:,i,j) );
      % for allelic effects
      modelmat = f2model2( igeno(:,i,j) );
    end
    
    [ b(j,i,:), v(j,i,:), rss ] = parameters( y, [ modelmat x ] );
    % lod( j, i ) = - ((df-ntr+1)/2) * log( rss/rawss );
    lod( j, i ) = - (n/2) * log( rss/rawss );    
  end 
end

const = max( lod(:) );

prob = exp( lod - const );
total = sum(prob(:));
prob = prob/total;

%lod = mean( exp( lod - const ) );
lod = lognormalmean( lod - const );
lod = log( lod ) + const;
lod = squeeze( lod );
tmp = lod ( find( lod~=0 ) );
bf = mean( exp( tmp ) );

mmm = b;
vvv = v;
ppp = prob;
nrep = length( prob(:) );
[mmm,vvv] = unconditional( reshape( b, nrep, modelsize), ...
    reshape( v, nrep, modelsize), prob(:) );

varargout = cell(1,2);

if( nargout==3 )
  varargout(1) = {lod};
elseif( nargout==4 )
    varargout(1) = {lod};
    varargout(2) = {bf};
end

