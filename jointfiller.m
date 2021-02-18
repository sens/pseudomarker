function igeno = jointfiller( geno, mlocs, fillocs, npages, cross )
% JOINTFILLER Function to impute missing genotypes at specified pairs of
% marker locations
%
% geno = matrix of marker genotypes
% mlocs = locations of the typed markers
% fillocs = locations of the markers to be filled in
% npages = number of sets of imputations  
% cross = cross type, 'bc, for BACKCROSS and 'f2' for INTERCROSS
%
% igenoA = left pseudomarker
% igenoB = right pseudomarker  
% loc = location of the pair
%  
% This program will impute marker genotypes at specified locations based
% on the typed genotypes and their positions.  The imputation is done
% from their JOINT PAIRWISE distribution, so use with care.
% 
  
% number of markers and number of individuals
[n,m] = size( geno ); 

% number of pseudomarker locations
npseudo = length( fillocs );

igeno = zeros( n, npseudo, npages ); % create dummy pseudomarker array

for( i=1:n ) % go through each individual
  
  % find the indices of the non-missing data
  typed = find( ~isnan( geno(i,:) ) );
  typedlocs = mlocs( typed ); % locations of typed loci
  typedgeno = geno( i, typed ); % genotypes of typed loci
  
  % find the flanking markers for the first pseudomarker
  [ lgt, ldt, rgt, rdt ] = whichone( typedgeno, typedlocs, [], [], ...
				     fillocs(1) );
  % impute and assign the genotypes
  first = midintrnd( lgt, rgt, ldt, fillocs(1), rdt, npages, cross );
  igeno( i, 1, : ) = first;
  
  for( j = 2:npseudo )
    % find the flanking markers for the second pseudomarker
    [ lgt, ldt, rgt, rdt ] = whichone( typedgeno, typedlocs, first, ...
				       fillocs(j-1), fillocs(j) );
    % impute and assign the genotypes
    first = midintrnd( lgt, rgt, ldt, fillocs(j), rdt, ...
		       npages, cross );
    igeno( i, j, : ) = first;
  end

end

% ------------------------------------------------------------------
function gt = midintrnd( leftg, rightg, leftd, loc, rightd, npages, cross )
% function to impute given genotypes of flanking markers; marker
% genotypes may be missing
% works only for the BACKCROSS right now

% either flanking marker may be missing in which case it is coded as NaN
% the corresponding locations are -Inf or Inf according as they are the
% left or right flanking marker  
  
theta1 = recomb( loc - leftd ); % length of left subinterval
theta2 = recomb( rightd - loc ); % length of right subinterval
%theta0 = recomb( rightd - leftd );

% let L, X and R be the genotypes of the left, middle and right ends
% we want to find P(X|L,R) = P(L,X,R)/P(L,R)
% p1 = P(X=1|L)
% p2 = P(R|X=1)
% P(L,X=1,R) = 0.5*p1*p2
% P(L,X=0,R) = 0.5*(1-p1)*(1-p2)
% P(L,R) = P(L,X=1,R) + P(L,X=0,R) = 0.5*p0
% hence P(X=1|L,R) = p1*p2 / p0

if( cross == 'bc' )
  % left subinterval coditional prob
  if( isnan( leftg ) ) % if left gtype is missing
    p1 = 0.5;
  else
    p1 = theta1*(1-leftg) + (1-theta1)*leftg;
  end
  
  % right subinterval conditional prob
  if( isnan( rightg ) ) % if right gtype is missing
    p2 = 0.5;
  else
    p2 = theta2*(1-rightg) + (1-theta2)*rightg;
  end
  
  %  p0 = theta0*abs( rightg - leftg ) + (1-theta0)*(1-abs(leftg-rightg));
  % prob of the flanking markers is equal to the sum of the probs of the
  % two cases when the middle genotype is 0 or 1
  p0 = p1.*p2 + (1-p1).*(1-p2);
  % conditional prob 
  prob = p1.*p2./p0;
  
  % there are two cases - when the left marker is a pseudomarker and when
  % it isn't; if the length of the prob vector is more than one, then the
  % left one is a pseudomarker
  
  if( length(prob)>1 ) % left one is pseudomarker
    gt = binornd( 1, prob ); % length of prob vector tells us how many to
			     % simulate 
  else
    gt = binornd( 1, prob, 1, npages );
  end

end

if( cross == 'f2' )

  if( length(leftg)==1 ) % left one is not pseudomarker
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
    
    gt = discreternd( w, npages ) - 1; % simulate the genotype based on the
				       % weights subtract 1 because
				       % discreternd simulates the index
				       
  else
    npages = length( leftg );
    p = zeros( npages, 3 );
    q = zeros( npages, 3 );
    
    % conditional probs of left interval
    p(:,1) = f2trans( 0, leftg, theta1 );
    p(:,2) = f2trans( 1, leftg, theta1 );
    p(:,3) = f2trans( 2, leftg, theta1 );
    
    % conditional probs of right interval
    q(:,1) = f2trans( 0, rightg, theta2 );
    q(:,2) = f2trans( 1, rightg, theta2 );
    q(:,3) = f2trans( 2, rightg, theta2 );
    
    w = p.*q; % weights vector for the three cases
    
    gt = discreternd2( w' ) - 1; % simulate the genotype based on the
				       % weights subtract 1 because
				       % discreternd simulates the index
    
  end  
end



% ------------------------------------------------------------------
function [ lgt, ldt, rgt, rdt ] = whichone( mgeno, mlocs, pgeno, ploc, loc )
% function to decide which genotypes and which locations flank a
% particular location; the markers and pseudomarkers are assumed to be
% ordered; this function works for ALL CROSSES
% the understanding is that all markers and pseudomarkers are ORDERED
%
% mgeno = marker genotypes for an individual
% mlocs = marker locations
% pgeno = pseudomarker genotype
% ploc = pseudomarker location
% loc = location to be imputed
%
% lgt = left marker genotype
% ldt = left marker location
% rgt = right marker genotype
% rdt = right marker location  
  
  
nmk = length( mlocs ); % number of marker loci
toleft = sum( mlocs<=loc ); % number of markers to the left of the
                            % desired location
toright = nmk - toleft; % number of markers to the left of the desired
                        % location 

if( toleft>0 )
  lgt = mgeno( toleft );
  ldt = mlocs( toleft );
else % if no markers to the left then missing
  lgt = NaN;
  ldt = -Inf;
end

if( toright>0 )
  rgt = mgeno( toleft+1 );
  rdt = mlocs( toleft+1 );
else % if no markers to the right then missing
  rgt = NaN;
  rdt = Inf;
end

if( ~isempty(pgeno) ) % if there is a pseudomarker to the left
  if( ploc > ldt ) % and if the pseudomarker is to the left of the typed
                   % marker replace the left marker by the pseudomarker
    lgt = pgeno; 
    ldt = ploc;
  end
end


% --------------------------------------------------------
function p = f2trans( x, y, theta )
% transition function for an intercross

  if(x==0)
    p = (y==0) * ((1-theta)^2) + (y==1) * ( 2*theta*(1-theta) )...
	+ (y==2) * (theta^2) + isnan(y);
  
  elseif(x==1)
    p = ( (y==0) + (y==2) )* (theta*(1-theta)) + ...
	(y==1) * (theta^2 + (1-theta)^2) + isnan(y);
   
  elseif(x==2)
    p = (y==2) * ((1-theta)^2) + (y==1) * ( 2*theta*(1-theta) )...
	+ (y==0) * (theta^2) + isnan(y);

  elseif(isnan(x))
    p = ( (y==0) + (y==2) ) * (1/4) + (y==1) * (1/2) + isnan(y);
    
  end
  
p = p';