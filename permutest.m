function p = permutest( y, x, z, pseudom, npermute, notmiss )
% PERMUTEST Perform a permutation test on MAINSCAN maximum
%
% PERMUTEST(Y,X,Z,FAKE,NPERMUTE)  
% PERMUTEST(Y,X,Z,FAKE,NPERMUTE,SUBSET)   
% M = PERMUTEST(Y,X,Z,FAKE,NPERMUTE)
% M = PERMUTEST(Y,X,Z,FAKE,NPERMUTE,SUBSET) 
% Permutes the phenotypes and the covariates NPERMUTE times and calculates
% the maximum value of the BLOD on the natural logarithm scale as obtained
% from MAINSCAN(Y,X,FAKE,SUBSET).  If the output arguement is not
% assigned, it will print out the 50th, 75th, 90th, 95th and 99th
% percentiles of both rows.
%  
% Y = phenotypes; CAN BE A MATRIX
% X = covariates, use [] if none  
% Z = interacting covariates, use [] if none  
% FAKE = pseudomarker data structure  
% NPERMUTE = number of permutations 
% SUBSET = not missing index; if missing, rows with no missing data will
%          be selected
%
% The output P is a matrix with two rows; the first row is the maximum of
% the BLOD score and the second is the maximum of the Bayes Factors.  The
% columns correspond to permutations.  If no outpur agrgument is selected
% the 50th, 75th, 90th, 95th and 99th percentiles of the LOD scores are
% returned. 
%  
% See also MAINSCAN, PERMUTEST2.
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/07/12 18:59:38 $	
  
  if( nargin == 5 )
    notmiss = phenoobs( [ y x z ] );
  end
 
  nchroms = length(pseudom);
  mmm = zeros(2,nchroms);  
  % try Bayes factor too
  p = zeros( 2, npermute );  

  pheno = y;
  covar = x;
  covarint = z;
  
  for( j=1:npermute )
    fprintf( '%d ', j );
    if( mod(j,20) == 0 )
      fprintf( '\n' );
    end
    
    ord = randperm(length(notmiss));
    pheno(notmiss,:) = y(notmiss(ord),:);
    
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

    
    lod = mainscan( pheno, covar, covarint, pseudom, notmiss );
    
    
    for( i=1:nchroms )
      mmm(1,i) = max(lod(i).lod);
      mmm(2,i) = max(lod(i).bf);      
    end
    
    p(1,j) = max(mmm(1,:));
    p(2,j) = max(mmm(2,:));
    
  end
  
  % return only the first row
  p = p(1,:);

  if( nargout == 0 )
    levels = [ 50 75 90 95 99 ];
    [ levels' prctile( p', levels ) ]
  end

