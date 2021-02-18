function flagout = flagbf1( lod1, cv )
% FLAGBF1 Flag one-locus models based on one-dimensional scans and Bayes 
%         factors
%
% FLAGBF1(LOD1,LOD2,CV)  
% LOD1 = output of MAINSCAN or MAINESTIMATE
% CV = critical value Bayes factor cutoff
%  
% See also: FLAGBF, PERMUTEST.  
  
  nchrom = length( lod1 ); % number of chromosomes
  penalty = sqrt(lod1(1).n); % the Bayes factor per parameter penalty
  cross = lod1(1).cross; % the type of cross being analyzed
  
  % set the penalties for the different kinds of crosses
  if ( cross == 'bc' )
    penalty = penalty; % 2 df vs 1 df
  elseif ( cross == 'f2' )
    penalty = penalty^2; % 3 df vs 1 df
  else
    error( 'Unrecognized cross type' )
  end
  
  % print the initial message
  fprintf( '\n\tThe following loci were flagged:\n\n' );
  fprintf( '\tLoci\t    BF\t\tmaxLOD\n' );

  for( i=1:nchrom )
    tmp = lod1(i).lod;
    mmm = max(tmp(:));
    % get the bayes factor
    bf = lod1(i).bf/penalty;
    if( bf > cv ) % if important print message
      fprintf( '\t%d\t%1.3e\t%4.3f\n',i,bf,mmm);
    end
  end




