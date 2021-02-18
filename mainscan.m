function onelod = mainscan( y, x, z, pseudom, notmiss )
% MAINSCAN Plot and compute a one-dimensional genome scan.
%
% ONELOD = MAINSCAN(Y,X,Z,FAKE)    
% ONELOD = MAINSCAN(Y,X,Z,FAKE,OBSDX)  
% Y = phenotypes; can be a matrix
% X = additive covariates, use [] if none  
% Z = interacting covariates, use [] if none    
% FAKE = pseudomarker data structure array 
% OBSDX = observed data index; if absent only individuals with
%         phenotypes ( and covariates ) not missing are analyzed
%  
% ONELOD = one dimensional lod score structure array which has the following
% fields 
% N = number of individuals
% CROSS = type of cross ('bc' or 'f2')
% CHRID = chromosome id
% BF = *unscaled* Bayes factor, the marginal distribution of data
% LOD = LOD score vector, base 10
% MPOS = pseudomarker positions

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.837 $ $Date: 2001/08/08 19:27:53 $	
  
  if( nargin == 4 )
    notmiss = phenoobs( [ y x z ] );
  end
  

  n = length( notmiss );
  
  y = y(notmiss,:);
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
  
  [ n, ntr ] = size(y);
  nchroms = length( pseudom );
  chroms = zeros(1,nchroms);
 
  chromlen = zeros( 1, nchroms );
  pos = [];
  for( i=1:nchroms )
    chroms(i) = pseudom(i).chrid;
    chromlen(i) = pseudom(i).chromlen;
  end

  cross = pseudom(1).cross;
  onelod = repmat( struct( 'n', n, 'cross', cross, 'chrid', 0, ...
			   'bf', 0, 'baselod', 0, 'lod', [], ...
			   'mpos', [] ), 1, nchroms );
  
  for( i=1:nchroms )
    [thislod, bf] = onescan( y, x, z, ...
			     pseudom(i).igeno(notmiss,:,:), cross );
    onelod(i).chrid = chroms(i);
    onelod(i).bf = bf;
    onelod(i).baselod = baselod/log(10);    
    onelod(i).lod = thislod/log(10);
    onelod(i).mpos = pseudom(i).mpos;
    onelod(i).chromlen = pseudom(i).chromlen;
  end


