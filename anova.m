function table = anova( y, x, z, fake, qtl, model, obsdx )
%  General function for scanning.
% 
% [LOD,BF] = ANOVA(Y,X,Z,FAKE)    
% [LOD,BF] = SCAN(Y,X,Z,FAKE,QTL)  
% [LOD,BF] = SCAN(Y,X,Z,FAKE,QTL,MODEL)  
% [LOD,BF] = SCAN(Y,X,Z,FAKE,QTL,MODEL,OBSDX)
%
% LOD = multidimensional array of LOD scores on the log base 10 scale
% BF = Bayes factor calculated by averaging the LOD scores
%  
% Y = vector of trait values  
% X = matrix of additive covariates, use [] if none
% Z = matrix of interacting covariates, use [] if none  
% FAKE = pseudomarker data structure; output of IMPUTE
% QTL = matrix of chromosome ID's and number of QTL in them
%       eg. [ 1 2 5; 2 1 1] which means that there are QTL on chromosomes
%       1 2 5 with two on chromosome 1; default is one QTL on each
%       chromosome
% MODEL = model of the QTLs; structure with fields twoint and threeint.
%         eg. model= struct( 'twoint', [ 1 2 ], 'threeint', [ 1 2 3 ] );
%         default is empty
% OBSDX = observed data index; default is the individuals with no missing
%         phenotypes 
% See also: ONESCAN, TWOSCAN, THREESCAN, WTAVERAGE.

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2001/09/24 21:59:57 $	
  
if( nargin < 7 )  
  obsdx = phenoobs( [ y x z ] );
end

if( nargin < 6 )
  model = [];
end

if( nargin < 5 )
  qtl = zeros( 2, length(fake) );
  for( i=1:length(fake) )
    qtl( 1:2, i ) = [ fake(i).chrid; 1 ];
  end
  %qtl = [ 1:length(fake); ones(1,length(fake)) ];
end

cross = fake(1).cross;

y = y(obsdx);
[ n, ntr ] = size(y);  

if( ~isempty(x) )
  [ n, p ] = size(x(obsdx,:));
end

npages = size( fake(1).igeno, 3 );

% begin device
% this is a device to make the indexing of the fakes work smoothly
maxchrid = max( qtl(1,:) );
% make a blank imputed data structure
newfake = repmat( struct( 'cross', cross, 'chrid', 0, 'igeno', [], 'mpos', ...
		       [], 'chromlen', 0 ), 1, maxchrid );
% for those indexes that match a real chromosome, make it equal to the
% appropriate fake element

for( i=1:length(fake) )
  newfake( fake(i).chrid ) = fake(i);
end

% old code for doing the same thing as the above loop
%for( i=1:maxchrid )
%  for( j=1:length(fake) )
%    if( fake(j).chrid == i )
%      newfake(i) = fake(j);
%    end
%  end
%end

% assign the new fake to the old one
fake = newfake;
% end device

[nposs,nmarkers] = ncombinations( qtl, fake );
nqtl = sum( qtl(2,:) );
chrid = zeros( 1, nqtl );
nchroms = size(qtl,2);
selchroms = qtl(1,:);
if( ~isempty(model) )
  if( isfield( model, 'twoint' ) )
    ntwoint = size( model.twoint, 1 );
  else 
    ntwoint = 0;
  end
  if ( isfield( model, 'threeint' ) )
    nthreeint = size( model.threeint, 1 );
  else 
    nthreeint = 0;
  end
else
  ntwoint = 0;
  nthreeint = 0;
end




% dimension of the lod array
dimlod = zeros(1,nqtl);
% the qtl numbering list
cumnqtl = [ 0 cumsum( qtl(2,:) ) ];
for( i=1:nchroms )
  dimlod( (cumnqtl(i)+1):cumnqtl(i+1) ) = nmarkers(i);
  chrid( (cumnqtl(i)+1):cumnqtl(i+1) ) = qtl(1,i);
end



% $$$ for(i=2:nposs)
% $$$   qtlindex = plusqtlindex( qtlindex, qtl, nmarkers )
% $$$ end


% initialize lod array
lod = zeros( [ npages dimlod ] );
% base rss
[rawss,df] = rss( y, ones(n,1) );     

if( cross == 'bc' )
  modelmat = ones( n, nqtl+1 );
elseif ( cross == 'f2' )
  modelmat = ones( n, 2*nqtl+1 );
end
covintmodelmat = [];


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
      itwoint=0;
      while( itwoint<ntwoint )
	itwoint = itwoint+1;
	qtl1 = model.twoint( itwoint, 1 );
	qtl2 = model.twoint( itwoint, 2 );
	mk1 = fake( chrid(qtl1) ).igeno( obsdx, qtlindex(qtl1), ipages );
	mk2 = fake( chrid(qtl2) ).igeno( obsdx, qtlindex(qtl2), ipages );
	modelmat( :, nqtl+1+itwoint ) = bcmodel(mk1).* bcmodel(mk2);
      end

      % for three-way interactions
      ithreeint=0;
      while( ithreeint<nthreeint )
	ithreeint = ithreeint+1;
	qtl1 = model.threeint( ithreeint, 1 );
	qtl2 = model.threeint( ithreeint, 2 );
	qtl3 = model.threeint( ithreeint, 3 );
	mk1 = fake( chrid(qtl1) ).igeno( obsdx, qtlindex(qtl1), ipages );
	mk2 = fake( chrid(qtl2) ).igeno( obsdx, qtlindex(qtl2), ipages );
	mk3 = fake( chrid(qtl3) ).igeno( obsdx, qtlindex(qtl3), ipages );
	modelmat( :, 1+nqtl+ntwoint+ithreeint ) = ...
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
      itwoint=0;
      while( itwoint<ntwoint )
	itwoint = itwoint+1;
	qtl1 = model.twoint( itwoint, 1 );
	qtl2 = model.twoint( itwoint, 2 );
	mk1 = fake( chrid(qtl1) ).igeno( obsdx, qtlindex(qtl1), ipages );
	mk2 = fake( chrid(qtl2) ).igeno( obsdx, qtlindex(qtl2), ipages );
	effind = (1+2*nqtl)+4*(itwoint-1) + (1:4);
	modelmat( :, effind ) = prod_model( f2model(mk1), f2model(mk2) );
      end

      % for three-way interactions
      ithreeint=0;
      while( ithreeint<nthreeint )
	ithreeint = ithreeint+1;
	qtl1 = model.threeint( ithreeint, 1 );
	qtl2 = model.threeint( ithreeint, 2 );
	qtl3 = model.threeint( ithreeint, 3 );
	mk1 = fake( chrid(qtl1) ).igeno( obsdx, qtlindex(qtl1), ipages );
	mk2 = fake( chrid(qtl2) ).igeno( obsdx, qtlindex(qtl2), ipages );
	mk3 = fake( chrid(qtl3) ).igeno( obsdx, qtlindex(qtl3), ipages );	
	effind = 1+2*nqtl+4*ntwoint+8*(nthreeint-1) + (1:8);	
	modelmat( :, effind ) = prod_model( prod_model( ...
	    f2model(mk1), f2model(mk2) ), f2model(mk3) );
      end
       
    end
    
    if( ~isempty(z) )
      covintmodelmat = prod_model( modelmat( :, 2:end ), z );
    end
    
    [ rss, df ] = rss( y, [ modelmat x covintmodelmat ] );    
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
