function [lod,bf] = emtwoscan( y, x, geno, mpos, qpos, cross, model )
% EMTWOSCAN Cmpute the approximate posterior distribution of QTL assuming
% there is two of them on the same chromosome. Uses the EM algorithm
% 
% y = vector of trait values  
% x = matrix of covariates, use [] if none
% geno = marker genotypes; can be missing
% mpos = marker positions
% qpos = positions where the LOD score has to be calculated  
% cross = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS
% model = model type, 'a+b' or 'a*b'
% THIS FUNCTION NEEDS TO BE DEBUGGED.  ONLY BACKCROSS IS SUPPORTED NOW.
% ADDITIVE MODEL SUPPORT IS TENTATIVE.
  
[ n, ntr ] = size(y);  % number of traits; not needed now
[ n, p ] = size(x); % number of covariates; not needed now
[ n, m ] = size( geno ); % size of genotype matrix
q = length( qpos ); % number of QTL positions to be tested

lod = zeros( q, q ); % initialize lod matrix

% calculate baselod
[mmm,sss] = ols( y, [ ones(n,1) x ] ); % mean and rss of null model
sigma = sqrt(sss/n); % estimate sd
% calculate individual likelihood contributions and then sum
a = normpdf( y, mmm, sigma ); 
baselod = sum( log( a ) )

% find number of parameters in the model 
if( cross == 'bc' )
  if( model == 'a+b' )
    npar = 3;
  elseif( model =='a*b' )
    npar = 4;
  end
elseif ( cross == 'f2' )
  if( model == 'a+b' )
    npar = 5;
  elseif( model == 'a*b' )
    npar = 9;
  end
end

% make model matrix
modelmat = ones( n, npar );
  
% start the EM algorithm
% first we have to have all the probabilities for the QTL genotypes at
% the hypothesized positions; then we calculate the a priori and the a
% posteriori positions

% apriori probs
Ap1 = zeros(n,q,q); % for left position g1
Ap2 = zeros(n,q,q); % for right position g2
Ap00 = zeros(n,q,q); % p(g1=0,g2=0)
Ap01 = zeros(n,q,q); % p(g1=0,g2=1)
Ap10 = zeros(n,q,q); % p(g1=1,g2=0)
Ap11 = zeros(n,q,q); % p(g1=1,g2=1)

% get the apriori probs by going through the individuals one by one
for( i=1:n )
  idx = find( ~isnan( geno( i, : ) ) ); % which indexes are not missing
  mgeno = geno( i, idx ); % genotypes that are not missing
  haspos = mpos( idx ); % the corresponding marker positions
  [ Ap1(i,:,:), Ap2(i,:,:), Ap11(i,:,:) ] = getprobs( mgeno, haspos, qpos );
end

% calculate the joints from the marginals
Ap10 = Ap1-Ap11;
Ap01 = Ap2-Ap11;
Ap00 = 1-Ap10-Ap01-Ap11;

% ----- end get genotype probs

% matrix relating the regression coefficients (beta) to the group means (mu)
if( model=='a*b' )
  D = [ 1 -1 -1  1; ...
	1  1 -1 -1; ...
	1 -1  1 -1; ...
	1  1  1  1];
elseif( model=='a+b' )
  D = [ 1 -1 -1; ...
	1  1 -1; ...
	1 -1  1; ...
	1  1  1];
end

% $$$ D = [ 1 0 0 0; ...
% $$$       1 1 0 0; ...
% $$$       1 0 1 0; ...
% $$$       1 1 1 1];


% now go through each position that is permissible and calculate the LOD
% score there

for( i=1:q )
  for( j=(i+1):q )

   % [ i j]
     beta = [ mmm; zeros( npar-1, 1 ) ]; % initialize the betas
     absdiff = 1; % difference between successive iterations
     rss = sss; % initialize rss
     
     % initialize a posteriori probs
     Bp00 = Ap00(:,i,j); 
     Bp10 = Ap10(:,i,j);
     Bp01 = Ap01(:,i,j);
     Bp11 = Ap11(:,i,j);
     
 
     % EM block
     % while the difference between iterations is big
     while( absdiff > 1e-5 )
       
       % E step
       % change beta's to the mu's
       mu = D * beta;
       % estimate of error variance
       sigma = sqrt(rss/n);
       
       % calculation of posterior expectations of each genotype class 
       Q00 = Ap00(:,i,j) .* normpdf( y, mu(1), sigma );
       Q10 = Ap10(:,i,j) .* normpdf( y, mu(2), sigma );
       Q01 = Ap01(:,i,j) .* normpdf( y, mu(3), sigma );
       Q11 = Ap11(:,i,j) .* normpdf( y, mu(4), sigma );
       
       S = Q00 + Q01 + Q10 + Q11;
       Bp00 = Q00./S;
       Bp01 = Q01./S;
       Bp10 = Q10./S;
       Bp11 = Q11./S;    
       
       Bp1 = Bp10 + Bp11;
       Bp2 = Bp01 + Bp11;
       
       % M step  
       % make model matrix
       if( model=='a*b' )
	 modelmat( :, 2:4 ) = [ bcmodel(Bp1) bcmodel(Bp2) ...
 		    4*Bp11-2*Bp1-2*Bp2+1 ];
	 W = matstar( n, Bp1, Bp2, Bp11 );       
       elseif( model=='a+b' )
	 modelmat( :, 2:3 ) = [ bcmodel(Bp1) bcmodel(Bp2) ];
	 W = matplus( n, Bp1, Bp2, Bp11 );              
       end

% $$$        modelmat( :, 2:4 ) = [ Bp1 Bp2 Bp11 ];
% $$$        W = matstar( n, Bp1, Bp2, Bp11 );
       beta1 = W\(modelmat'*y);
       rss = y'*( y - 2*modelmat*beta1 ) + beta1'*W*beta1;
       % difference between the two estimates
       absdiff = sum( abs( beta - beta1 ) );
       beta = beta1;
     end;
    
% $$$     % Haley and Knott block
% $$$     modelmat( :, 2:4 ) = [ bcmodel(Ap1(:,i,j)) bcmodel(Ap2(:,i,j)) ...
% $$$ 		    4*Ap11(:,i,j)-2*Ap1(:,i,j)-2*Ap2(:,i,j)+1 ];
% $$$     [ beta, rss ]  = ols( y, modelmat );
% $$$     sigma = sqrt(rss/n);
% $$$     lod( i,j ) = sum ( log( normpdf( y, modelmat*beta, sigma ) ) ) ...
% $$$          - baselod;
    
    % EM block
     lod( i, j ) = sum( log( Q00 + Q01 + Q10 + Q11 ) ) - baselod;
     
  end
end
  
  tmp = lod ( find( lod~=0 ) );
  bf = mean( exp( tmp ) );

% ------
function [ p1, p2, p12 ] = getprobs( mgeno, mlocs, qlocs )
% mgeno = marker genotypes observed
% mlods = positions of the observed marker genotypes
% qlocs = locations at which the probability of genotype 1 is wanted  

  nmk = length( mgeno ); % number of markers with observed genotypes
  q = length( qlocs ); % number of postions wanted
  p1 = zeros( q, q ); % initialize output
  p2 = zeros( q, q ); % initialize output
  p12 = zeros( q, q ); % initialize output
  
  % go through the wanted positions
  for( i=1:q )
    % number of typed positions to the left
    toleft1 = sum( mlocs <= qlocs(i) );
    % number of typed positions to the right
    toright1 = nmk - toleft1;

    if( toleft1 == 0 )
      g1 = NaN;
      d1 =  -Inf;
    else 
      g1 = mgeno( toleft1 );
      d1 = qlocs(i) - mlocs( toleft1 );
    end

    if( toright1 == 0 )
      g2 = NaN;
      d2 = Inf;
    else 
      g2 = mgeno( toright1 );
      d2 = mlocs( toright1 ) - qlocs(i) ;
    end
    
    for( j=(i+1):q )

      % number of typed positions to the left
      toleft2 = sum( mlocs <= qlocs(j) );
      % number of typed positions to the right
      toright2 = nmk - toleft2;

      if( toleft2 == 0 )
	g3 = NaN;
	d3 = -Inf;
      else 
	g3 = mgeno( toleft2 );
	d3 = qlocs(j) - mlocs( toleft2 );
      end
      
      if( toright2 == 0 )
	g4 = NaN;
	d4 = Inf;
      else 
	g4 = mgeno( toright2 );
	d4 = mlocs( toright2 ) - qlocs(j);
      end

      if( toleft1 == toleft2 )
%	mlocs
%	[ d1 d2 d3 d4]
%	[ i j ]
	dmid = qlocs(j)-qlocs(i);
	p1(i,j) = condprob2( g1, g2, d1, d2 );
	p2(i,j) = condprob2( g3, g4, d3, d4 );
	p12(i,j) = condprob3( g1, g2, d1, dmid, d4 );
      else
	p1(i,j) = condprob2( g1, g2, d1, d2 );
	p2(i,j) = condprob2( g3, g4, d3, d4 );	
	p12(i,j) = p1(i,j)*p2(i,j);
      end
      
      % end of j loop
    end
    % end of i loop
  end
  


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

% ------------
function p = condprob3( g1, g3, d1, d2, d3 )

  p1 = intervalprob( g1, d1 );
  p2 = intervalprob( 1, d2 );
  p3 = intervalprob( g3, d3 );

  p11 = p1*p2*p3;
  p01 = (1-p1)*(1-p2)*p3;
  p10 = p1*(1-p2)*(1-p3);
  p00 = (1-p1)*p2*(1-p3);  
  
  p = p11 / ( p11 + p01 + p10 + p00 );
  
% ----------
function m = matplus( n, p1, p2, p12 )

  fff = sum( [ p1 p2 p12 ] );
  sss(1:2) = 2*fff(1:2)-n;
  sss(3) = 4*fff(3)-2*fff(1)-2*fff(2)+n;
  m = [ n sss(1) sss(2) ; ...
	sss(1) n sss(3) ; ...
 	sss(2) sss(3) n ];
  
% ----------
function m = matstar( n, p1, p2, p12 )
  
  fff = sum( [ p1 p2 p12 ] );
  sss(1:2) = 2*fff(1:2)-n;
  sss(3) = 4*fff(3)-2*fff(1)-2*fff(2)+n;
  
  m = [ n sss(1) sss(2) sss(3); ...
 	sss(1) n sss(3) sss(2); ...
 	sss(2) sss(3) n sss(1); ...
 	sss(3) sss(2) sss(1) n ];
  
% $$$   m = [ n      fff(1) fff(2) fff(3);
% $$$ 	fff(1) fff(1) fff(3) fff(3);
% $$$ 	fff(2) fff(3) fff(2) fff(3);
% $$$ 	fff(3) fff(3) fff(3) fff(3)];
  