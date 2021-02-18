function twolod = pairscan( y, x, z, pseudom, notmiss )
% PAIRSCAN Plot a pairwise genome scan for both additive and interaction
%           effects.
%
% TWOLOD = PAIRSCAN(Y,X,Z,FAKE)  
% TWOLOD = PAIRSCAN(Y,X,Z,FAKE,OBSDX)  
% Y = phenotypes; can be a matrix
% X = additive covariates, use [] if none  
% Z = interacting covariates, use [] if none    
% M = marker data structure array
% FAKE = pseudomarker data structure array 
% OBSDX = observed data index
%  
% TWOLOD = two-dimensional lod score structure array which has the following
% fields 
% N = number of individuals
% CROSS = type of cross ('bc' or 'f2')
% CHRID1 = chromosome id corresponding to the row index
% CHRID2 = chromosome id corresponding to the col index
% CHROMLEN1 = chromosome length corresponding to the row index
% CHROMLEN2 = chromosome length corresponding to the col index
% MPOS1 = pseudomarker positions corresponding to the row index
% MPOS2 = pseudomarker posutions corresponding to the col index
% if row index is equal to the column index  
%     BFADD = *unscaled* Bayes factor for the additive model
%     BFINT = *unscaled* Bayes factor for the full model  
% else   
%     BF = *unscaled* Bayes factor, the marginal distribution of data;
%     if row index is greater than the col index, this corresponds to
%     full model, else the additive model
% LOD = LOD score vector, base 10

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.834 $ $Date: 2001/08/08 19:31:22 $	
  
  % what to do if the missing data index is missing
  if( nargin == 4 )
    notmiss = phenoobs( [ y x z ] );
  end
  
  y = y(notmiss,:);
  n = length(notmiss);
  
  if( ~isempty(x) )
    x = x(notmiss,:);
    rssx = rss(y,[ ones(n,1) x]);
    baselod = -(n/2)*log(rssx);
  else
    baselod = 0;
  end
  
  if( ~isempty(z) )
    z = z(notmiss,:);
  end
  
  % spacing between pseudomarkers
  spacing = pseudom(1).mpos(2) - pseudom(1).mpos(1);
  
  nchroms = length( pseudom ); % number of chromosomes
  
  
  % vector for the length of chromosomes
  chromlen = zeros( 1, nchroms ); 
  
  % position vector for plot
  pos = [];
  % number of pseudomarkers
  npseudo = [];
  chroms = zeros(1,nchroms);
  
  % get all the chromosome lengths and the number of pseudomarkers in
  % each chromosome
  for( i=1:nchroms )
    chroms(i) = pseudom(i).chrid;
    chromlen(i) = pseudom(i).chromlen;
    npseudo = [ npseudo length(pseudom(i).mpos) ] ;
  end

  cross = pseudom(1).cross;
  twolod = repmat( struct( 'n', n, 'cross', cross, ...
			    'chrid1', 0, 'chrid2', 0, 'bfadd', 0, ...
			    'bfint', 0, 'bf', 0, 'baselod', 0, 'lod', [], ...
			    'chromlen1', 0, 'chromlen2', 0, 'mpos1', [], ...
			    'mpos2', [] ), nchroms, nchroms );
  
  for( i=1:nchroms )
    for( j=i:nchroms )

      fprintf( '(%d,%d)', pseudom(i).chrid, pseudom(j).chrid );
      if(j==nchroms)
	fprintf('\n');
      end
      
      if(i==j)
	[thislod, bfi] = twoscan( y, x, z, pseudom(i).igeno(notmiss,:,:), ...
				 cross, '*' );
	[thatlod, bfa] = twoscan( y, x, z, pseudom(i).igeno(notmiss,:,:), ...
				 cross, 'a+b' );
	twolod( i, i ).baselod = baselod/log(10);
	twolod( i, i ).lod = (thislod+thatlod')/log(10);
	twolod( i, i ).chromlen1 = pseudom(i).chromlen;	
	twolod( i, i ).chromlen2 = pseudom(i).chromlen;		
	twolod( i, i ).bfadd = bfa;
	twolod( i, i ).bfint = bfi;
	twolod( i, i ).mpos1 = pseudom(i).mpos;
	twolod( i, i ).mpos2 = pseudom(i).mpos;	
	twolod( i, i ).chrid1 = pseudom(i).chrid;
	twolod( i, i ).chrid2 = pseudom(i).chrid;	
      else
	[thislod, bfi] = twoscan2( y, x, z, pseudom(i).igeno(notmiss,:,:), ...
				  pseudom(j).igeno(notmiss,:,:), cross, ...
				  '*' );
	[thatlod, bfa] = twoscan2( y, x, z, pseudom(i).igeno(notmiss,:,:), ...
				  pseudom(j).igeno(notmiss,:,:), cross, ...
				  'a+b' );
	twolod( i, j ).baselod = baselod/log(10);
	twolod( j, i ).baselod = baselod/log(10);
	twolod( i, j ).lod = thislod/log(10);
	twolod( j, i ).lod = thatlod'/log(10);	
	twolod( i, j ).chromlen1 = pseudom(i).chromlen;	
	twolod( i, j ).chromlen2 = pseudom(j).chromlen;		
	twolod( j, i ).chromlen1 = pseudom(j).chromlen;	
	twolod( j, i ).chromlen2 = pseudom(i).chromlen;		
	twolod( j, i ).bf = bfa;
	twolod( i, j ).bf = bfi;	
	twolod( i, j ).mpos1 = pseudom(i).mpos;
	twolod( i, j ).mpos2 = pseudom(j).mpos;	
	twolod( j, i ).mpos1 = pseudom(j).mpos;
	twolod( j, i ).mpos2 = pseudom(i).mpos;	
	twolod( i, j ).chrid1 = pseudom(i).chrid;
	twolod( i, j ).chrid2 = pseudom(j).chrid;	
	twolod( j, i ).chrid1 = pseudom(j).chrid;
	twolod( j, i ).chrid2 = pseudom(i).chrid;	
      end
      
      
    end

  end
  
