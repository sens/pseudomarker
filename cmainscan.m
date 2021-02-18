function onelod = cmainscan( y, x, pseudom, cqtl, notmiss )
% MAINSCAN Plot and compute a one-dimensional genome scan.
%
% ONELOD = MAINSCAN(Y,X,FAKE)    
% ONELOD = MAINSCAN(Y,X,FAKE,CQTL)    
% ONELOD = MAINSCAN(Y,X,FAKE,CQTL,OBSDX)  
% Y = phenotypes; can be a matrix
% X = covariates, use [] if none  
% FAKE = pseudomarker data structure array 
% CQTL = chromosome numbers and approximate positions to condition on;
%        [ chrid; mpos ]  
% OBSDX = observed data index; if absent only individuals with
%         phenotypes ( and covariates ) not missing are analyzed
%  
% ONELOD = one dimensional lod score structure array which has the following
% fields 
% N = number of individuals
% CROSS = type of cross ('bc' or 'f2')
% CHRID = chromosome id
% BF = *unscaled* Bayes factor, the marginal distribution of data
% LOD = LOD score vector
% MPOS = pseudomarker positions

  if( nargin == 3 )
    notmiss = phenoobs( [ y x ] );
  end
  

  n = length( notmiss );
  
  pheno = y(notmiss,:);
  if( ~isempty(x) )
    covar = x(notmiss,:);
    r = resid( pheno, [ ones(n,1) covar ] );
    %    rssx = rss(y,[ ones(n,1) x]) / rss(y,ones(n,1))
    rssx = rss(pheno,[ ones(n,1) covar]);
    baselod = -(n/2)*log(rssx);
    y = r;
  else
    baselod = 0;
  end
  
  [ n, ntr ] = size(y);
  nchroms = length( pseudom );
  chroms = zeros(1,nchroms);
 
  % spacing between pseudomarkers
  spacing = pseudom(1).mpos(2) - pseudom(1).mpos(1);
  
  chromlen = zeros( 1, nchroms );
  pos = [];
  for( i=1:nchroms )
    chroms(i) = pseudom(i).chrid;
    chromlen(i) = pseudom(i).chromlen;
  end

  % find our cross type
  cross = pseudom(1).cross;
  % make the lod data structure
  onelod = repmat( struct( 'n', n, 'cross', cross, 'chrid', 0, ...
			   'bf', 0, 'lod', [], 'mpos', [] ), 1, nchroms );
  
  % make the imputed data structure by getting pseudomarker genotypes on
  % the chromosomes and positions of interest; we are assuming that no
  % more than one QTL on one chromosome is being conditioned on
  cfake = getmarkergeno( pseudom( cqtl(1,:) ), cqtl(2,:) );
  
  % go through each chromosome 
  for( i=1:nchroms )

    % what are the ids of other chromosomes that need to be conditioned
    % on; this will give us a list of indexes
    otheridx = find( cqtl(1,:)~=i );
    
    % make the imputed data structure corresponding other pseudomarkers to
    % condition on
    cfake = getmarkergeno( pseudom( cqtl(1,otheridx) ), cqtl(2,otheridx) );
    
    % combine the other pseudomarkers with the chromosome of interest 
    fake = [ pseudom(i) cfake ];
    
    % make qtl data structure
    % the first row is a list of the chromosomes of interest 
    % the second row is the number of QTL on each; we assume only one on
    % each 
    qtl = [ i cqtl(1,otheridx); 1 ones(length(otheridx),1) ];
    [thislod, bf] = scan( y, [], fake, qtl, [], notmiss );
    qtla = [ cqtl(1,otheridx); ones(length(otheridx),1) ];
    [thisloda, bfa] = scan( y, [], fake(2:end), qtla, [], notmiss );
    onelod(i).chrid = chroms(i);
    onelod(i).bf = bf;
    % onelod(i).lod = thislod;
    onelod(i).lod = squeeze(thislod) - squeeze(thisloda);    
    onelod(i).mpos = pseudom(i).mpos;
    onelod(i).chromlen = pseudom(i).chromlen;
  end
