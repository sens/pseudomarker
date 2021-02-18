function [ add, int ] = getbf2( lod )
% GETBF2 Get bayes factors from the output of a pairscan.
%
% [ADD,INT]=GETBF2(LOD)
% ADD = Bayes factors for the additive model.
% INT = Bayes factors for the model with interactions.  
% LOD = output from pairscan.
%
% Note that these Bayes factors are really the marginal distribution on
% the chromosomes of interest and are hence not normalized.
% 
% See also GETBF, PAIRSCAN.
  
  [ nr nc ] = size( lod );
  
  add = zeros( nr, nc );
  int = zeros( nr, nc );
  
  for( i=1:nr )
    for( j=i:nc )
      if( i==j )
	int(i,j) = lod(i,i).bfint;
	add(i,j) = lod(i,i).bfadd;	
      else
	int(i,j) = lod(i,j).bf;
	add(i,j) = lod(j,i).bf;		
      end
      
    end
  end