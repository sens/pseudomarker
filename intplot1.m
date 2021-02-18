function intplot1( pheno, mk, mnames, labels )
% INTPLOT1 Plot interaction effects for one-QTL model.
%         
% INTPLOT1(PHENO,MK)
% INTPLOT1(PHENO,MK,MNAMES)    
% INTPLOT1(PHENO,MK,MNAMES,LABELS)
% PHENO = vector of phenotypes
% MK = vector of marker (or pseudomarker) genotypes  
% MNAMES = names of markers; default {'Marker1','Marker2'}
% LABELS = labels for genotypes; default {'A','B'} for backcross and
%          {'A','H','B'} for intercross
%
% If there are any missing phenotypes, those will be discarded. Whether the
% genotypes are markers or pseudomarkers is decided from the dimension of
% MK1 and MK2.  The rows correspond to individuals and columns to
% imputations.  The mean of each genotype combination is a weighted mean the
% mean over all imputations.  The variance is also computed in the same way.
% The weight that each imputation gets is proportional to RSS^(-n/2),
% where RSS is the residual sum of squares for that imputation.

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/09/24 23:14:48 $	

  % find dimensions of the markers/pseudomarkers
  [n,npages] = size(mk);
  
  % if marker genotypes, subset to get non-missing individuals
  if( npages==1 )
    % find which individuals have non-missing phenotypes and genotypes
    obsdx = phenoobs( [pheno mk] );
    % subset the phenotypes
    pheno = pheno(obsdx);
    mk = mk(obsdx);
  else
  % if pseudomarker genotypes, subset to get non-missing individuals    
    % find which individuals have non-missing phenotypes
    obsdx = phenoobs( pheno );
    % subset the phenotypes
    pheno = pheno(obsdx);
    mk = mk(obsdx,:);
  end

  % default marker names
  if( nargin<3 )
    mkname = 'Marker';
  else
    mkname = mnames;
  end
  
  % guess cross type and number of levels of the genotype factor
  if( max(mk) == 1 )
    cross = 'bc';
    nlevel = 2;
    
    % default labels for genotypes
    if(nargin<4)
      labels = {'A', 'B'};
    end
  
  elseif( max(mk) == 2 )
    cross = 'f2';
    nlevel = 3;
    % default labels
    if(nargin<4)
      labels = {'A', 'H', 'B'};    
    end
  end

  % make legend labels using the marker names
  % for( i=1:nlevel)
  %  legendlabels{i} = strcat( mk2name, '=', labels{i} );
  % end
  
  % number of genotype groups
  ngroup = nlevel;
  % group means initialized
  grmeans = zeros( ngroup, npages );
  % group vars initialized
  grvars = ones( ngroup, npages );
  wt = zeros( npages, 1 );  
  % make group ids corresponding to genotype combinations
  groupid = mk + 1;
  
  for( i = 1:npages )
    p = asfactor2( groupid(:,i) );
    % size(p)
    [b,r,v] = ols( pheno, p );
    % size(b)
    grmeans(:,i) = b;
    grvars(:,i) = (r/n)*diag(inv(v));  
    wt(i) = -(n/2)*log(r);
  end
  mx = max(wt(:));
  wt = exp(wt-mx);
  wt = wt/sum(wt(:));

  % mean is weighted group means
  grm = grmeans*wt;
  meanvar = grvars*wt;
  varmean = (grmeans - repmat(grm,1,npages)).*...
	    (grmeans - repmat(grm,1,npages)) *wt;
  grsd = sqrt( varmean + meanvar );
  
  % make markers for the intplots
    markers = 'o';
  
  % plus minus
  pm = [ -2 +2 ];

  % handle to plots that will be legended initialized
  h = zeros(1,nlevel);
  % make the plots now

  % for( j = 1:nlevel )

    % ids of the lines for a fixed level of marker 2
    id = 1:nlevel;
    % plot the connecting means with the marker and marker size
    % specified; handle is assigned
    h = plot( 0:(nlevel-1), grm(id), markers(1) );

    hold on;
    % plot connecting lines between means
    plot( 0:(nlevel-1), grm(id) );

    for( i = 1:nlevel )
      % for each fixed level of marker 1
      % get the group id
      id = i;
      % plot the error bars
      plot( repmat( i-1, 1, 2 ), grm(id) + pm*grsd(id) );
      % plot whiskers of error bars
      plot( i-1 + pm/50, repmat( grm(id) + 2*grsd(id), 1, 2 ) );
      plot( i-1 + pm/50, repmat( grm(id) - 2*grsd(id), 1, 2 ) );
    end
  
  % set the ticks on the x axis
  set( gca, 'xtick', 0:(nlevel-1) );
  xlim( [ -1/2 (nlevel -1/2) ] )
  % label ticks on x axis
  set( gca, 'xticklabel', labels );  
  % label the x axis
  xlabel( strcat( mkname, ' genotype' ) );    





