function newmk = subsetgeno3( m, npages )
% SUBSETGENO3 Subsets a marker genotype structure array by imputation number 
%
% M1 = SUBSETGENO1(M0,CHRID,NPAGE)
% M0 = mother genotype structure array returned from READDATA or imputed
%      structure array returned by IMPUTE
% NPAGES = imputation numbers
%
% This function is useful for studying genotype patterns or for
% subsetting individuals.
% 
% Example: subsetgeno2( m, [ 1 5 8 ] );
%          subsetgeno2( m, [ 1 5 8 ], [ 0.41 0.23 0.34 ] );
%          subsetgeno2( m, [ 1 5 8 ], [ 0.41 0.78; 0.32 0.45 0.63 0.92 ] );  
% See also: SUBSETGENO1.
%  
    
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/09/24 22:29:41 $	

% is this a marker or pseudomarker structure
  if( isfield( m(1), 'geno' ) )
    imputed = 0;
  elseif ( isfield( m(1), 'igeno' ) )
    imputed = 1;
  else
    error( 'This is not a marker data structure.' );
  end
  
  % number of chromosomes and number of individuals
  nchroms = length(m);
  if( imputed == 0 )
    n = size(m(1).geno,1);
  else
    n = size(m(1).igeno,1);
  end
  
  % subset the chromosomes according to chromosome id
  newmk = [];
  for( i=1:nchroms )  
    newmk = [ newmk m(i) ];
  end
  
  nchroms = length( newmk );

  for( i = 1:nchroms )
      if( imputed==1 )
	newmk(i).igeno = m(i).igeno(:,:,npages);
      end
  end

  




