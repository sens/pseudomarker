function bf = getbf( lod )
% GETBF Get bayes factors from the output of a MAINSCAN
%
% GETBF(LOD)  
% LOD = output from MAINSCAN
% 
% Note that these Bayes factors are really the marginal distribution on
% the chromosomes of interest and are hence not normalized.
%  
%  See also GETBF2, MAINSCAN.


  nr = length( lod );
  
  bf = zeros( 1, nr );
  
  for( i=1:nr )
    bf(i) = lod(i).bf;
  end
  