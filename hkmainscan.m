function onelod = hkmainscan( y, x, qtldata, notmiss, spacing, cross)
% EMMAINSCAN Plot and compute a one-dimensional genome scan using the
%            Haley-Knott method.
%            
% ONELOD = HKMAINSCAN(Y,X,M)      
% ONELOD = HKMAINSCAN(Y,X,M,OBSDX)    
% ONELOD = HKMAINSCAN(Y,X,M,OBSDX,SPACING)  
% ONELOD = HKMAINSCAN(Y,X,M,OBSDX,SPACING,CROSS)    
% Y = phenotypes; can be a matrix
% X = covariates, use [] if none  
% M = marker data structure array
% OBSDX = observed data index
% SPACING = spacing in Morgans of the LOD evaluation positions; if missing
%           0.02 Morgans.  The positions at which the LOD score is
%           evaluated is always between the first and last typed marker.
% CROSS = type of cross; 'bc' or 'f2'; if missing it is guessed
%  
% ONELOD = one dimensional lod score structure array which has the following
% fields 
%       n = number of individuals
%       cross = type of cross ('bc' or 'f2')
%       chrid = chromosome id
%       bf = *unscaled* Bayes factor, the marginal distribution of data
%       lod = LOD score vector
%       mpos = pseudomarker positions
%
% This data structure is compatible with the plotting functions
% associated with MAINSCAN.  
%  
% SEE ALSO: MAINSCAN, EMMAINSCAN, PLOTMAINSCAN, PLOTMAINLOCALIZE.
%


  % what to do if the observed data index is missing
  if( nargin < 4 )
    notmiss = phenoobs( [ y x ] ); 
  end
  
  % spacing between pseudomarkers
  if( nargin < 6 )
    spacing = 0.02;
  end

  % what to do if the missing data index is absent
  if( nargin < 6 )
    cross = guesscross( qtldata );
  end
  
  n = length( notmiss );
  
  y = y(notmiss,:);
  if( ~isempty(x) )
    x = x(notmiss,:);
    r = resid( y, [ ones(n,1) x ] );
    %    rssx = rss(y,[ ones(n,1) x]) / rss(y,ones(n,1))
    rssx = rss(y,[ ones(n,1) x]);
    baselod = -(n/2)*log(rssx);
    y = r;
  else
    baselod = 0;
  end
  
  [ n, ntr ] = size(y);
  nchroms = length( qtldata );
  chroms = zeros(1,nchroms);
 
  
  chromlen = zeros( 1, nchroms );
  pos = [];
  for( i=1:nchroms )
    chroms(i) = qtldata(i).chrid;
    chromlen(i) = qtldata(i).chromlen;
  end

  onelod = repmat( struct( 'n', n, 'cross', cross, 'chrid', 0, ...
			   'bf', 0, 'lod', [], 'mpos', [] ), 1, nchroms );
  
  for( i=1:nchroms )
    [thislod, bf] = hkonescan( y, [], qtldata(i).geno(notmiss,:), ...
	       qtldata(i).mpos, ...
	       qtldata(i).mpos(1):spacing:(qtldata(i).mpos(end)), cross );
    onelod(i).chrid = chroms(i);
    onelod(i).bf = bf;
    onelod(i).lod = thislod;
    onelod(i).mpos = qtldata(i).mpos(1):spacing:(qtldata(i).mpos(end));
    onelod(i).chromlen = qtldata(i).chromlen;
  end
