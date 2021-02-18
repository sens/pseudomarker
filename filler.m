function igeno = filler( geno, mlocs, fillocs, npages, cross )
% FILLER Function to impute missing genotypes at specified marker
%        locations.
%
% FILLER(GENO,MLOCS,FILLOCS,NPAGES,CROSS)  
% GENO = matrix of marker genotypes
% MLOCS = locations of the typed markers
% FILLOCS = locations of the markers to be filled in
% NPAGES = number of sets of imputations  
% CROSS = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS
%  
% This program will impute marker genotypes at specified locations based
% on the typed genotypes and their positions.  The imputation is done
% from their MARGINAL distribution, so use with care.
%
  
% number of markers and number of individuals
[n,m] = size( geno ); 

npseudo = length( fillocs ); % number of pseudomarkers

igeno = zeros( n, npseudo, npages ); % create dummy pseudomarker array

for( i=1:n ) % go through each individual

  left = 0; % left marker number
  leftg = NaN; % genotype missing
  leftd = -Inf; % distance -infinity 

  % find the indices of the non-missing data
  typed = find( ~isnan( geno(i,:) ) ); 
  if( length(typed)>=1 ) % if there are some typed markers
    right = typed(1); % index of first typed marker
    rightg = geno(i,right); % genotype of first typed marker
    rightd = mlocs(right); % location of first typed marker
  else
    right = m+1; % untyped marker
    rightg = NaN; % genotype missing
    rightd = Inf; % located at +infinity
  end
  
  % cycle through the pseudomarkers
  for( j=1:npseudo-1 )
    % simulate the marker genotype given the flanking genotypes and locations
    igeno( i, j, : ) = midintrnd( leftg, rightg, leftd, fillocs(j), rightd, ...
                                  npages );

    if( rightd<fillocs(j+1) ) % if next pseudomarker is to the right of the 
                         % right end of the current marker interval
      left = right;  
      leftd = rightd;
      leftg = rightg;
      typed(1) = []; % remove the first element in the typed marker list
      
      if( length(typed)>=1 )
	right = typed(1);
	rightg = geno(i,right);
	rightd = mlocs(right);
      else
	right = m+1;
	rightg = NaN;
	rightd = Inf;
      end
    end
  end
  igeno( i, npseudo, : ) = midintrnd( leftg, rightg, leftd, ...
				      fillocs(npseudo), rightd, npages );
  
end


% ------------------------------------------------------------------
function gt = midintrnd( leftg, rightg, leftd, loc, rightd, npages, cross )
% function to impute given genotypes of flanking markers; marker
% genotypes may be missing
% works only for the BACKCROSS and INTERCROSS right now

theta1 = recomb( loc - leftd );
theta2 = recomb( rightd - loc );
%theta0 = recomb( rightd - leftd );

if( cross = 'bc' )
  
  if( isnan( leftg ) ) % if left gtype is missing
    p1 = 0.5;
  else
    p1 = theta1*(1-leftg) + (1-theta1)*leftg;
  end
  
  if( isnan( rightg ) ) % if rightt gtype is missing
    p2 = 0.5;
  else
    p2 = theta2*(1-rightg) + (1-theta2)*rightg;
  end
  
  %  p0 = theta0*abs( rightg - leftg ) + (1-theta0)*(1-abs(leftg-rightg));
  % prob of the flanking markers is equal to the sum of the probs of the
  % two cases when the middle genotype is 0 or 1
  p0 = p1*p2 + (1-p1)*(1-p2);
  % conditional prob 
  prob = p1*p2/p0;
  
  gt = binornd( 1, prob, 1, npages );

end

if( cross = 'f2' )

  p = zeros( 1, 3 );
  q = zeros( 1, 3 );

  % conditional probs of left interval
  p(1) = f2trans( leftg, 0, theta1 );
  p(2) = f2trans( leftg, 1, theta1 );
  p(3) = f2trans( leftg, 2, theta1 );

  % conditional probs of right interval
  q(1) = f2trans( 0, rightg, theta2 );
  q(2) = f2trans( 1, rightg, theta2 );
  q(3) = f2trans( 2, rightg, theta2 );
  
  w = p.*q; % weights vector for the three cases

  gt = discreternd( w ) - 1; % simulate the genotype based on the weights
                             % subtract 1 because discreternd simulates the
                             % index
end

% --------------------------------------------------------
function p = f2trans( x, y, theta )
% transition function for an intercross

  switch x
   
   case 0
    switch y
     case 0
      p = (1-theta)^2;
     case 1
      p = 2*theta*(1-theta);
     case 2
      p = theta^2;
     case NaN
      p = 1;
    end    
   
   case 1
    switch y
     case 0
      p = theta*(1-theta);
     case 1
      p = theta^2 + (1-theta)^2;
     case 2
      p = theta*(1-theta);
     case NaN
      p = 1;
    end
   
   case 2
    switch y
     case 0
      p = theta^2;
     case 1
      p = 2*theta*(1-theta);
     case 2
      p = (1-theta)^2;
     case NaN
      p = 1;
    end
    
     case NaN
      switch y
       case 0
	p = 1/4;
       case 1
	p = 1/2;
       case 2
	p = 1/4;
       case NaN
	p = 1;
      end
      
  end
  


