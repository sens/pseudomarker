function [lod,bf] = estimate( y, x, fake, qtl, varargin )
% ESTIMATE General function for getting estimated effects.
% 
% [LOD,BF] = ESTIMATE(Y,X,FAKE)    
% [LOD,BF] = ESTIMATE(Y,X,FAKE,QTL)  
% [LOD,BF] = ESTIMATE(Y,X,FAKE,QTL,...)  
%
% LOD = multidimensional array of LOD scores on the log base 10 scale
% BF = Bayes factor calculated by averaging the LOD scores
%  
% Y = vector of trait values  
% X = matrix of covariates, use [] if none
% FAKE = pseudomarker data structure; output of IMPUTE
% QTL = matrix of chromosome ID's and number of QTL in them
%       eg. [ 1 2 5; 2 1 1] which means that there are QTL on chromosomes
%       1 2 5 with two on chromosome 1; default is one QTL on each
%       chromosome
%  
% The other arguments appear in name-value pairs.  To specify which pair of
% QTL will be interacting use the field 'int2', with the value being a
% matrix with two columns.  The matrix [ 1 2 ] indicates that the first two
% QTL interact.  To specify which triplet of QTL will be interacting use the
% field 'int3', with the value being a matrix with two columns.  The matrix
% [ 1 2 3 ] indicates that the first three QTL interact.  The observed data
% index can be given as the field 'obsdx'.
%  
% See also: ONEESTIMATE, TWOESTIMATE, THREEESTIMATE, SCAN, WTAVERAGE.

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/07/12 18:42:16 $	

  % get the variable arguments
  nargin = length(varargin);
  if( length(varargin)>0 )
    nstep = 1;
    while( nstep <= nargin )
      argtype = varargin{nstep};
      switch argtype
       case 'int2'
	int2 = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'int3'
	int3 = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'obsdx'
	obsdx = varargin{ nstep+1 };
	nstep = nstep+2;
       otherwise
	error( 'Cannot recognize option.' )
	nstep = nargin;
      end
    end
  end
  
  % default obsdx
  if( ~exist(obsdx) )  
    obsdx = phenoobs( [ y x ] );
  end
  
  % default int2
  if( ~exist(int2) )  
    int2 = [];
    nint2 = 0;
  end
  
  % default int3
  if( ~exist(int3) )  
    int3 = [];
    nint3 = 0;
  end
  
  % default qtls
  if( ~exist(qtl) )
    qtl = zeros( 2, length(fake) );
    for( i=1:length(fake) )
      qtl( 1:2, i ) = [ fake(i).chrid; 1 ];
    end
  end

  % number of qtl
  nqtl = sum( qtl(2,:) );
  % type of cross
  cross = fake(1).cross;
  
  if( cross == 'bc' )
    % number of main effects
    nmain = nqtl;
    % number if two-way interactions
    nint2 = size(twoint,1);
    % number if three-way interactions    
    nint3 = size(threeint,1);  
    % number of effects to estimate
    neff = 1 + nmain + nint2 + nint3;
    % means initialized
    mmm = zeros( neff, 1 );
    % variances initialized
    vvv = zeros( neff, 1 );
  elseif( cross == 'f2' )
    % number of main effects
    nmain = 2 * nqtl;
    % number if two-way interactions
    nint2 = 4 * size(twoint,1);
    % number if three-way interactions    
    nint3 = 8 * size(threeint,1);  
    % number of effects to estimate
    neff = 1 + nmain + nint2 + nint3;
    % means initialized
    mmm = zeros( neff, 1 );
    % variances initialized
    vvv = zeros( neff, 1 );
  else
    error('Cannot recognize cross.');
  end
  
  % subset the phenotypes
  y = y(obsdx);
  [ n, ntr ] = size(y);  
  
  % subset the covariates if any
  if( ~isempty(x) )
    [ n, p ] = size(x(obsdx,:));
  end
  
  % number of pseudomarker imputations
  npages = size( fake(1).igeno, 3 );
  
  % begin device
  % this is a device to make the indexing of the fakes work smoothly
  maxchrid = max( qtl(1,:) );
  % make a blank imputed data structure
  newfake = repmat( struct( 'cross', cross, 'chrid', 0, 'igeno', [], ...
			    'mpos',  [], 'chromlen', 0 ), 1, maxchrid );
  % for those indexes that match a real chromosome, make it equal to the
  % appropriate fake element
  
  for( i=1:length(fake) )
    newfake( fake(i).chrid ) = fake(i);
  end
  
  % assign the new fake to the old one
  fake = newfake;
  % end device
  
  % number of markers and number of combinations possible
  [nposs,nmarkers] = ncombinations( qtl, fake );
  % chromosome ids
  chrid = zeros( 1, nqtl );
  % number of chromosomes
  nchroms = size(qtl,2);
  % selected chromosomes that appear in the model
  selchroms = qtl(1,:);
  
  
  % dimension of the lod array
  dimlod = zeros(1,nqtl);
  % the qtl numbering list
  cumnqtl = [ 0 cumsum( qtl(2,:) ) ];
  for( i=1:nchroms )
    dimlod( (cumnqtl(i)+1):cumnqtl(i+1) ) = nmarkers(i);
    chrid( (cumnqtl(i)+1):cumnqtl(i+1) ) = qtl(1,i);
  end
  
  
  % initialize lod array
  lod = zeros( [ npages dimlod ] );
  % base rss
  [rawss,df] = rss( y, ones(n,1) );     
  
  if( cross == 'bc' )
    modelmat = ones( n, nqtl+1 );
  elseif ( cross == 'f2' )
    modelmat = ones( n, 2*nqtl+1 );
  end
  
  % initialize the index of qtls here
  qtlindex = initqtlindex( qtl );
  % number of possibilities covered
  iposs=0;
  
  % the loop begins here
  while( iposs<nposs )
    
    for( ipages = 1:npages )
      
      % make the model matrices
      if( cross=='bc' )
	for( ichrom=1:nchroms )
	  % the qtl numbers of che qtl in the ichrom-th chromosome
	  iii = (cumnqtl(ichrom)+1) : cumnqtl(ichrom+1);
	  % make model matrix
	  modelmat( :, iii+1 ) =  bcmodel( ...
	      fake(selchroms(ichrom)).igeno( obsdx, qtlindex(iii), ipages ) );
	end
	
	% for two-way interactions
	iint2=0;
	while( iint2<nint2 )
	  iint2 = iint2+1;
	  qtl1 = int2( iint2, 1 );
	  qtl2 = int2( iint2, 2 );
	  mk1 = fake( chrid(qtl1) ).igeno( obsdx, qtlindex(qtl1), ipages );
	  mk2 = fake( chrid(qtl2) ).igeno( obsdx, qtlindex(qtl2), ipages );
	  modelmat( :, nqtl+1+iint2 ) = bcmodel(mk1).* bcmodel(mk2);
	end
	
	% for three-way interactions
	iint3=0;
	while( iint3<nint3 )
	  iint3 = iint3+1;
	  qtl1 = int3( iint3, 1 );
	  qtl2 = int3( iint3, 2 );
	  qtl3 = int3( iint3, 3 );
	  mk1 = fake( chrid(qtl1) ).igeno( obsdx, qtlindex(qtl1), ipages );
	  mk2 = fake( chrid(qtl2) ).igeno( obsdx, qtlindex(qtl2), ipages );
	  mk3 = fake( chrid(qtl3) ).igeno( obsdx, qtlindex(qtl3), ipages );
	  modelmat( :, 1+nqtl+nint2+iint3 ) = ...
	      bcmodel(mk1).* bcmodel(mk2).*bcmodel(mk3);
	end
	
	
      elseif( cross=='f2' )
	
	% main effects
	iqtl = 0;
	while( iqtl<nqtl )
	  iqtl = iqtl+1;
	  mk = fake( chrid( iqtl ) ).igeno( obsdx, qtlindex(iqtl), ipages );
	  effind = 1+2*(iqtl-1) + (1:2);
	  % modelmat( :, 1+effind ) =  f2model( mk );
	  modelmat( :, effind ) =  f2model( mk );	
	end
	
	% for two-way interactions
	iint2=0;
	while( iint2<nint2 )
	  iint2 = iint2+1;
	  qtl1 = int2( iint2, 1 );
	  qtl2 = int2( iint2, 2 );
	  mk1 = fake( chrid(qtl1) ).igeno( obsdx, qtlindex(qtl1), ipages );
	  mk2 = fake( chrid(qtl2) ).igeno( obsdx, qtlindex(qtl2), ipages );
	  effind = (1+2*nqtl)+4*(iint2-1) + (1:4);
	  modelmat( :, effind ) = prod_model( f2model(mk1), f2model(mk2) );
	end
	
	% for three-way interactions
	iint3=0;
	while( iint3<nint3 )
	  iint3 = iint3+1;
	  qtl1 = int3( iint3, 1 );
	  qtl2 = int3( iint3, 2 );
	  qtl3 = int3( iint3, 3 );
	  mk1 = fake( chrid(qtl1) ).igeno( obsdx, qtlindex(qtl1), ipages );
	  mk2 = fake( chrid(qtl2) ).igeno( obsdx, qtlindex(qtl2), ipages );
	  mk3 = fake( chrid(qtl3) ).igeno( obsdx, qtlindex(qtl3), ipages );	
	  effind = 1+2*nqtl+4*nint2+8*(nint3-1) + (1:8);	
	  modelmat( :, effind ) = prod_model( prod_model( ...
	      f2model(mk1), f2model(mk2) ), f2model(mk3) );
	end
	
      end

      [ mmm(ndx,:), vvv(ndx,:), rss ] = parameters( y, [ modelmat x ] );
      % [ rss, df ] = rss( y, [ modelmat x ] );
      ndx = mysub2ind( [ ipages qtlindex ], [ npages dimlod ] );
      lod( ndx ) = - (n/2) * log( rss/rawss ); 
      %qtlindex
    end
    
    iposs = iposs+1;
    
    %lod = mean( exp( lod - const )
    if( iposs<nposs)
      qtlindex = plusqtlindex( qtlindex, qtl, nmarkers );
    end
    
  end
  
  if( size(lod,1) > 1 )
    const = max( lod(:) );
    %lod = mean( exp( lod - const ) );
    %lod = lognormalmean( lod - const );
    zeroidx = ( lod > 0 ); % which entries are non-zero
    zeroidx = sum( zeroidx ); % number of lod entries are non-zero that are
			      % non zero
    idx = find( zeroidx > 0 ); % which lod entries are non-zero for at
                             % least one imputation
    lod = wtaverage( lod - const );  
    lod = log( lod ) + const;
    lod = squeeze( lod );
  else
    idx = ( lod > 0 ); % which entries are non-zero
  end
  
  %lod = squeeze( lod );
  %tmp = lod ( find( lod ~=0 ) );
  tmp = lod ( idx );
  bf = mean( exp( tmp ) );
  
  % convert to base 10
  lod = lod / log(10);
  
  
  
function [nposs,nmarkers] = ncombinations( qtl, fake )
% Function to compute the number of combinations to go through

  nmarkers = zeros(1,size(qtl,2));
  nposs = 1;
  for( i=1:size(qtl,2) )
    nmarkers(i) = length( fake(qtl(1,i)).mpos );
    for( j=1:qtl(2,i) )
      nposs = nposs * (nmarkers(i)-j+1)/j;
    end
  end
  
  
function init = initqtlindex( qtl )
% Function to initialize the QTL index function  

  init = zeros( 1, sum(qtl(2,:))' );
  qtlnum = 1;
  for( i=1:size(qtl,2) )
    for( j=1:qtl(2,i) )
      init( qtlnum ) = j;
      qtlnum = qtlnum+1;
    end     

  end

function ind = plusqtlindex( ind, qtl, nmarkers )
% Function to advance the QTL index function by one step

  cumqtlnum = cumsum(qtl(2,:));
  
  qtlnum = length( ind ); % qtl number
  stopornot = 0; % whether to stop the loop back or not
  whichchrom = size(qtl,2); % which chromosome is being considered
  whichone = qtl(2,whichchrom); % which qtl in a chromosome

  while( stopornot ~= 1 )

    if( ind(qtlnum) == nmarkers(whichchrom)-qtl(2,whichchrom)+whichone )
      % if the index is equal to the maximum possible value
      ind(qtlnum) = 0; % floor the index
      stopornot = -1; % indicate the a backtracking is to be done
      whichone = whichone - 1; % go to previous index

      if( whichone == 0 ) % if the last index in the chromosome has been
                          % reached 
	ind(qtlnum:(qtlnum+qtl(2,whichchrom)-1)) = 1:qtl(2,whichchrom);
	whichchrom = whichchrom - 1; % go to previous chromosome
	whichone = qtl(2,whichchrom); 

	if (whichchrom == 0 )
	  stopornot = 1
	end
      
      end
      qtlnum = qtlnum - 1; % go to previous index

      
    else

      if (stopornot == 0 )
	ind(qtlnum) = ind(qtlnum)+1;      
	stopornot = 1;
      elseif( stopornot == -1 )
	ind(qtlnum) = ind(qtlnum)+1;
	
	for( i=qtlnum+1:cumqtlnum(whichchrom) )
	  ind(i) = ind(i-1)+1;       
% $$$        
% $$$        for( ichrom=1:nchroms )
% $$$  	modelmat( :, (2*cumnqtl(ichrom)+2):(2*cumnqtl(ichrom+1)+1) ) = ...
% $$$  	    f2model( fake(ichrom).igeno(:,qtlindex,ipages) );
% $$$        end
      

	end
	
	stopornot = 1;
      end

    
    end
  
  end
  
      
  
function ndx = mysub2ind( sub, siz )
% $$$   n = length(sub);
% $$$   k = [1 cumprod(siz(1:end-1))];
% $$$   ndx = 1;
% $$$   for i = 1:n,
% $$$     ndx = ndx + (sub{i}-1)*k(i);
% $$$   end

k = [1 cumprod(siz(1:end-1))];
ndx = 1 + sum( (sub-1).*k );
