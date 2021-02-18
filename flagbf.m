function flagout = flagbf( lod1, lod2, cv )
% FLAGBF Flag two-locus models based on one- and two-dimensional scans
%      based on Bayes factors
%
% FLAGBF(LOD1,LOD2,CV)  
% LOD1 = output of MAINSCAN or MAINESTIMATE
% LOD2 = output of PAIRSCAN
% CV = 1x3 vector of critical values  
%      [ overall int coat ]
%      overall = overall critical value determined from PAIRPERMUTE2;
%                this should be on the LOD (log base 10 ) scale
%      int = bayes factor cutoff for interactions; default 10
%      coat =  bayes factor cutoff for coat tail; default 10
%  
% See also: FLAGLOD, PERMUTEST2.  
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/09/24 20:55:08 $	
  
  % cv(1) = cv(1)*log(10); % default is to have cutoff on the log 10 scale
  
  nchrom = length( lod1 ); % number of chromosomes
  penalty = sqrt(lod1(1).n); % the Bayes factor per parameter penalty
  cross = lod1(1).cross; % the type of cross being analyzed
  
  % set the penalties for the different kinds of crosses
  if ( cross == 'bc' )
    penaltyint = penalty; % 4 df vs 3 df
    penaltyadd = penalty; % 3 df vs 2 df
  elseif ( cross == 'f2' )
    penaltyint = penalty^4; % 9 df vs 5 df
    penaltyadd = penalty^2; % 5 df vs 3 df
  else
    error( 'Unrecognized cross type' )
  end
  
  % print the initial message
  fprintf( '\n\tThe following loci were flagged:\n\n' );
  fprintf( '\tLoci\tLODmax\tBF\n' );

  for( i=1:nchrom )
    for( j=i:nchrom )
      tmp = lod2(i,j).lod;
      mmm = max( tmp(:) );

      % does the max exceed limits?
      if( mmm > cv(1) )

	% get the bayes factor for the interaction
	if( i==j )
	  bf = lod2(i,j).bfint/(lod2(i,j).bfadd*penaltyint);
	else 
	  bf = lod2(i,j).bf/(lod2(j,i).bf*penaltyint);
	end
	
	if( bf > cv(2) ) % if interaction important print message
	    lab = strcat( '\t', int2str(i),'x',int2str(j), '\t', ...
			  num2str(mmm), '\t', num2str(bf), '\n' );
	    fprintf( lab );
	end
	
	%else % else check for coat tails

	% get the bayes factors for additive model
	if( i==j )
	  bf = lod2(i,j).bfadd/(lod1(i).bf*penaltyadd);
	else 
	  bf1 = lod2(j,i).bf/(lod1(i).bf*penaltyadd);
	  bf2 = lod2(j,i).bf/(lod1(j).bf*penaltyadd);
	  bf = min( [ bf1 bf2 ] );
	end
	
	if( bf > cv(3) ) % if additive effect exists
	  lab = strcat( '\t', int2str(i),'+',int2str(j), '\t', ...
			num2str(mmm), '\t', num2str(bf), '\n' );
	  fprintf( lab );
	end
	
      end
	
    end
      
  end
  
