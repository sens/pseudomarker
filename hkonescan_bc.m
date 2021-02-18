function [lod,bf] = hkonescan_bc( y, x, geno, mpos, qpos )
% HKONESCAN Cmpute the approximate posterior distribution of QTL assuming
% there is only one of them. Uses the Haley-Knott algorithm
% 
% y = vector of trait values  
% x = matrix of covariates, use [] if none
% geno = marker genotypes; can be missing
% mpos = marker positions
% qpos = positions where the LOD score has to be calculated  
% cross = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS

[ n, ntr ] = size(y);  
[ n, p ] = size(x);
[ n, m ] = size( geno );
q = length( qpos );

lod = zeros( 1, q );

%[rawss,df] = rss( y, ones(n,1) );
[mmm,sss] = ols( y, [ ones(n,1) x ] );
sigma = sqrt(sss/n);
a = normpdf( y, mmm, sigma );
baselod = sum( log( a ) )

  npar = 2;
  modelmat = ones( n, npar );
  
% start the EM algorithm
% first we have to have all the probabilities for the QTL genotypes at
% the hypothesized positions; then we calculate the a priori and the a
% posteriori positions

aprior = zeros(n,q);
apost = aprior;

for( i=1:n )
  idx = find( ~isnan( geno( i, : ) ) ); % which indexes are not missing
  mgeno = geno( i, idx );
  haspos = mpos( idx );
  aprior( i, : ) = getprobs( mgeno, haspos, qpos );
end
apost = aprior;
% ----- end get genotype probs



for( i=1:q )

%  i
  beta = [ mmm zeros( npar-1, 1 ) ]';
  absdiff = 1;
  rss = sss;
  
  % Haley and Knott block
  modelmat( :, 2 ) = bcmodel( aprior(:,i) );
  [ beta, rss ]  = ols( y, modelmat );
  sigma = sqrt(rss/n);
  lod( i ) = sum ( log( normpdf( y, modelmat*beta, sigma ) ) ) - baselod;

end

tmp = lod ( find( lod~=0 ) );
bf = mean( exp( tmp ) );


% ------
function probs = getprobs( mgeno, mlocs, qlocs )
% mgeno = marker genotypes observed
% mlods = positions of the observed marker genotypes
% qlocs = locations at which the probability of genotype 1 is wanted  

  nmk = length( mgeno ); % number of markers with observed genotypes
  q = length( qlocs ); % number of postions wanted
  probs = zeros( 1, q ); % initialize output
  
  if ( nmk==0 ) % if no markers typed
    probs = repmat( 0.5, 1, q );
  else
    % go through the wanted positions
    for( i=1:q )
      % number of typed positions to the left
      toleft = sum( mlocs <= qlocs(i) );
      % number of typed positions to the right
      toright = nmk - toleft;

      
      % if some typed marker to the left
      if( toleft>0 )
	lgt = mgeno( toleft );
	ldt = mlocs( toleft );
	theta1 = recomb( qlocs(i) - ldt );  
	p1 = theta1*(1-lgt) + (1-theta1)*lgt;
      else % if no markers to the left then missing
	p1 = 1;
      end
      
      % if some typed marker to right
      if( toright>0 )
	rgt = mgeno( toleft+1 );
	rdt = mlocs( toleft+1 );
	theta2 = recomb( rdt - qlocs(i) );	
	p2 = theta2*(1-rgt) + (1-theta2)*rgt;	
      else % if no markers to the right then missing
	p2 = 1;
      end

       if( ( toleft > 0 ) & ( toright == 0 ) )      
	 probs(i) = p1;
       elseif( ( toleft == 0 ) & ( toright > 0 ) )      	
	 probs(i) = p2;
       else		
	 p0 = p1*p2 + (1-p1)*(1-p2);
	 probs(i) = p1*p2/p0;
       end
    end
  end

% -----------

% -----------
function p = intervalprob( gt, d )

  if( isnan(gt) )
    p = 0.5;
  else
    theta = recomb( d );
    p = theta * (1-gt) + (1-theta) * gt;
  end
  
% ------------
function p = condprob2( g1, g2, d1, d2 )
  
  p1 = intervalprob( g1, d1 );
  p2 = intervalprob( g2, d2 );

  a = p1*p2;
  b = (1-p1)*(1-p2);
  p = a/(a+b);

  
% ----------
function m = matstar( n, p )
  
  fff = sum( p );
  sss = 2*fff-n;
  m = [ n+10 sss; sss n+10 ];

  
  

