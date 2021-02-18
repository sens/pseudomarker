function [lod,bf] = emonescan_f2( y, x, geno, mpos, qpos )
% EMONESCAN Cmpute the approximate posterior distribution of QTL assuming
% there is only one of them. Uses the EM algorithm.  For F2 only.
% 
% y = vector of trait values  
% x = matrix of covariates, use [] if none
% geno = marker genotypes; can be missing
% mpos = marker positions
% qpos = positions where the LOD score has to be calculated  

  
[ n, ntr ] = size(y);  
[ n, p ] = size(x);
[ n, m ] = size( geno );
q = length( qpos );

lod = zeros( 1, q );

%[rawss,df] = rss( y, ones(n,1) );
[mmm,sss] = ols( y, [ ones(n,1) x ] );
sigma = sqrt(sss/n);
a = normpdf( y, mmm, sigma );
baselod = sum( log( a ) );

npar = 3;
modelmat = ones( n, npar );
  
% start the EM algorithm
% first we have to have all the probabilities for the QTL genotypes at
% the hypothesized positions; then we calculate the a priori and the a
% posteriori positions

aprior = zeros(n,q,3); % rows are individuals, cols are qlocs, pages are
                       % gtypes 
apost = aprior;

for( i=1:n )
%  tmp = [ geno( j, : ); mpos ];
  idx = find( ~isnan( geno( i, : ) ) ); % which indexes are not missing
  mgeno = geno( i, idx );
  haspos = mpos( idx );
%  fprintf( 'Obs %d\n', i );
  tmp = getprobs( mgeno, haspos, qpos );
  aprior( i, :, : ) = tmp;
end
apost = aprior;


for( i=1:q )
%  i
  beta = [ mmm mmm mmm ]';
  absdiff = 1;
  rss = sss;
  
  % EM block
  while( absdiff > 1e-4 )
    % EM steps
    % E step
    % change beta's to the mu's
    % this assumes that the coding is of "contr.treatment"
    mu = beta;
    % estimate of error variance
    sigma = sqrt(rss/n);
    % calculation of posterior expectations
    a = aprior( :, i, 1 ) .* normpdf( y, mu(1), sigma );
    b = aprior( :, i, 2 ) .* normpdf( y, mu(2), sigma );
    c = aprior( :, i, 3 ) .* normpdf( y, mu(3), sigma );    
    d = a+b+c;
    
    apost( :, i, 1 ) = a./d;
    apost( :, i, 2 ) = b./d;
    apost( :, i, 3 ) = c./d;    
    
    % M step  
    % make model matrix
    modelmat = reshape(apost(:,i,:),n,3);
    ttt = sum( modelmat );
    W = diag( ttt );
    beta1 = W\(modelmat'*y);
    rss = y'*( y - 2*modelmat*beta1 ) + beta1'*W*beta1;
    % difference between the two estimates
    absdiff = sum( abs( beta - beta1 ) );
    beta = beta1;
  end
  
  % EM block
  mu = beta;
  a = aprior( :, i, 1 ) .* normpdf( y, mu(1), sigma );
  b = aprior( :, i, 2 ) .* normpdf( y, mu(2), sigma );
  c = aprior( :, i, 3 ) .* normpdf( y, mu(3), sigma );    
  d = a+b+c;

  % Haley and Knott block
% $$$ modelmat( :, 2 ) = bcmodel( aprior(:,i) );
% $$$ [ beta, rss ]  = ols( y, modelmat );
% $$$ sigma = sqrt(rss/n);
% $$$ lod( i ) = sum ( log( normpdf( y, modelmat*beta, sigma ) ) ) - baselod;

  % EM block
  lod( i ) = sum(log(d)) - baselod;

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
  probs = zeros( q, 3 ); % initialize output

  if ( nmk==0 ) % if no markers typed
    probs = repmat( [0.25 0.5 0.25], q, 1 );
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
	p1 = f2trans( lgt, [ 0 1 2 ], theta1 );
      else % if no markers to the left then missing
	%p1 = [ 1 1 1 ];
      end
      
      % if some typed marker to right
      if( toright>0 )
	rgt = mgeno( toleft+1 );
	rdt = mlocs( toleft+1 );
	theta2 = recomb( rdt - qlocs(i) );	
	p2 = f2trans( rgt, [ 0 1 2 ], theta2 );
      else % if no markers to the right then missing
	%p2 = [ 1 1 1 ];
      end

       if( ( toleft > 0 ) & ( toright == 0 ) )      
	 probs(i,:) = p1;
       elseif( ( toleft == 0 ) & ( toright > 0 ) )      	
	 probs(i,:) = p2;
       else		
	 p0 = sum(p1.*p2); 
	 probs(i,:) = (p1.*p2/p0);
       end
    end
  end


  
% --------------------------------------------------------
function p = f2trans( x, y, theta )
% transition function for an intercross

  switch x
   
   case 0
    p = (y==0)*(1-theta)^2 + (y==1)*2*theta*(1-theta) + ...
	(y==2)*theta^2;
   case 1
    p = (y==0)*theta*(1-theta) + (y==1)* ( theta^2 + (1-theta)^2  ) + ...
	(y==2)*theta*(1-theta);
    
   case 2
    p = (y==0)*theta^2 + (y==1)*2*theta*(1-theta) + ...
	(y==2)*(1-theta)^2;
    
  end



