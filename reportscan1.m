function reportscan1( lod1, cv )
% REPORTSCAN1 Report loci deemed interesting
%
% REPORTSCAN1(LOD1,CV)  
% LOD1 = output of MAINSCAN or MAINESTIMATE
% CV = LOD score cutoff determined from permutation tests or otherwise
%  
% See also: FLAGlod, PERMUTEST.  
 
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/09/27 22:25:54 $	
  
  nchrom = length( lod1 ); % number of chromosomes
  cross = lod1(1).cross; % the type of cross being analyzed
  
  % set the penalties for the different kinds of crosses
  if ( cross == 'bc' )
    df = 1; % 2 df vs 1 df
  elseif ( cross == 'f2' )
    df = 2; % 3 df vs 1 df
  else
    error( 'Unrecognized cross type' )
  end

  % print the initial message
  fprintf( '\n\tThe following loci were flagged:\n\n' );
  fprintf( '\tChrom\tcM\t LOD\t  pval\n' );

  for( i=1:nchrom )
    tmp = lod1(i).lod;
    [ mmm, idx ] = max( tmp(:) );
    % does the max exceed limits?
    if( mmm > cv )
      pq = 1-chi2cdf(2*mmm*log(10),df);
      cM = 100*lod1(i).mpos(idx);
      fprintf( '\t%d\t%1.0f\t%1.3f\t%f\n', i, cM, mmm, pq );
    end
  end  

