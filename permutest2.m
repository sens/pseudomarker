function p = permutest2( y, x, z, pseudom, npermute, notmiss )
% PERMUTEST2 Perform a permutation test on PAIRSCAN maximum
%
% PERMUTEST2(Y,X,Z,FAKE,NPERMUTE)  
% PERMUTEST2(Y,X,Z,FAKE,NPERMUTE,SUBSET)   
% M = PERMUTEST2(Y,X,Z,FAKE,NPERMUTE)
% M = PERMUTEST2(Y,X,Z,FAKE,NPERMUTE,SUBSET) 
% Permutes the phenotypes and the covariates NPERMUTE times and calculates
% the maximum value of the BLOD on the natural logarithm scale as obtained
% from PAIRSCAN(Y,X,FAKE,SUBSET).  If the output arguement is not
% assigned, it will print out the 50th, 75th, 90th, 95th and 99th
% percentiles of both rows.
%  
% Y = phenotypes; CAN BE A MATRIX
% X = covariates, use [] if none  
% Z = interacting covariates, use [] if none   
% MDATA = qtl data structure  
% FAKE = pseudomarker data structure  
% NPERMUTE = number of permutations 
% SUBSET = not missing index; if missing, rows with no missing data will
%          be selected
%
% The output M is a matrix with two rows; the first row is the maximum of
% the BLOD score corresponding to the full pairwise model with interactions
% and the second row corresponds to the maximum BLOD score corresponding to
% the interaction term.  If no output argument is selected the 50th, 75th,
% 90th, 95th and 99th percentiles of the LOD scores corresponding to the
% full model are returned.
%  
% See also PAIRSCAN, PERMUTEST.

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2001/09/21 17:09:03 $	


  if( nargin == 5 )
    notmiss = phenoobs( [ y x z ] );
  end
 
  nchroms = length(pseudom);
  mmmfull = zeros(nchroms,nchroms);
  mmmint = zeros(nchroms,nchroms);
  p = zeros( 2, npermute );
  
  for( k=1:npermute )
    k
    ord = randperm(length(notmiss));
    pheno(notmiss,:) = y(notmiss(ord),:);
    
    %notmiss = notmiss(ord);
    
    if( ~isempty(x) )
      covar(notmiss,:) = x(notmiss(ord),:);
    else
      covar = [];
    end
    
    if( ~isempty(z) )
      covarint(notmiss,:) = z(notmiss(ord),:);
    else
      covarint = [];
    end

    lod = pairscan( pheno, covar, covarint, pseudom, notmiss );
    
    
    for( i=1:nchroms )
      for( j=i:nchroms )
	if( i==j )
 	  int = triu(lod(i,j).lod);
 	  add = tril(lod(i,j).lod)';
 	  mmmfull(i,j) = max(int(:));
 	  mmmint(i,j) = max(int(:)-add(:));
	else
 	  int = lod(i,j).lod;
 	  add = lod(j,i).lod';
 	  mmmfull(i,j) = max(int(:));
 	  mmmint(i,j) = max(int(:)-add(:));
	end	  
      end
    end
    
    p(1,k) = max(mmmfull(:));
    p(2,k) = max(mmmint(:));
  end


  if( nargout == 0 )
    levels = [ 50 75 90 95 99 ];
    [ levels' prctile( p', levels ) ]
  end
