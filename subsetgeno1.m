function mdata = subsetgeno1( m, subset )
% SUBSETGENO1 Subsets a marker genotype structure array based on individuals
%
% M1 = SUBSETGENO1(M0,SUBSET)
% M0 = mother genotype structure array returned from READDATA or imputed
%      structure array returned by IMPUTE
% SUBSET = index of the individuals whose genotypes have to be extracted
%
% This function is useful for studying genotype patterns or for
% subsetting individuals.
% 
% See also: PLOTGENO, READDATA, SUBSETGENO, SUBSETGENO2, SUBSETGENO3.
%  

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/09/24 21:18:47 $	
  
  nchrom = length(m);
  mdata = m;

  if( isfield( m(1), 'geno' ) )
    imputed = 0;
  elseif ( isfield( m(1), 'igeno' ) )
    imputed = 1;
  else
    error( 'This is not a marker data structure.' );
  end
    
  for( i=1:nchrom )
    if( imputed==0 )
      mdata(i).geno = m(i).geno(subset,:);
    else
      mdata(i).igeno = m(i).igeno(subset,:,:);
    end
  end
