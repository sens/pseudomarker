function intplot2( pheno, mk1, mk2, mnames, labels )
% INTPLOT2 Plot interaction effects for two-QTL model.
%         
% INTPLOT2(PHENO,MK1,MK2)
% INTPLOT2(PHENO,MK1,MK2,MNAMES)    
% INTPLOT2(PHENO,MK1,MK2,MNAMES,LABELS)
% PHENO = vector of phenotypes
% MK1 = vector of first marker (or pseudomarker) genotypes  
% MK2 = vector of second marker (or pseudomarker) genotypes    
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
%	$Revision: 0.833 $ $Date: 2001/09/24 23:03:28 $	

  % find dimensions of the markers/pseudomarkers
  [n1,npages1] = size(mk1);
  [n2,npages2] = size(mk2);  
  
  % if they do not match, print error message
  if( (n1~=n2) | (npages1~=npages2) )
    error( 'Marker or pseudomarker dimensions do not match.' )
  end
  
  % set number of individuals
  n = n1;
  % set number of imputations
  npages = npages1;
  
  % if marker genotypes, subset to get non-missing individuals
  if( npages==1 )
    % find which individuals have non-missing phenotypes and genotypes
    obsdx = phenoobs( [pheno mk1 mk2] );
    % subset the phenotypes
    pheno = pheno(obsdx);
    mk1 = mk1(obsdx);
    mk2 = mk2(obsdx);  
  else
  % if pseudomarker genotypes, subset to get non-missing individuals    
    % find which individuals have non-missing phenotypes
    obsdx = phenoobs( pheno );
    % subset the phenotypes
    pheno = pheno(obsdx);
    mk1 = mk1(obsdx,:);
    mk2 = mk2(obsdx,:);  
  end

  % default marker names
  if( nargin<4 )
    mk1name = 'Marker1';
    mk2name = 'Marker2';
  else
    mk1name = mnames{1};
    mk2name = mnames{2};      
  end
  
  % guess cross type and number of levels of the genotype factor
  if( max(mk1) == 1 )
    cross = 'bc';
    nlevel = 2;
    
    % default labels for genotypes
    if(nargin<5)
      labels = {'A', 'B'};
    end
  
  elseif( max(mk1) == 2 )
    cross = 'f2';
    nlevel = 3;
    % default labels
    if(nargin<5)
      labels = {'A', 'H', 'B'};    
    end
  end

  % make legend labels using the marker names
  for( i=1:nlevel)
    legendlabels{i} = strcat( mk2name, '=', labels{i} );
  end
  
  % number of genotype groups
  ngroup = nlevel*nlevel;
  % group means initialized
  grmeans = zeros( ngroup, npages );
  % group vars initialized
  grvars = ones( ngroup, npages );
  wt = zeros( npages, 1 );  
  % make group ids corresponding to genotype combinations
  groupid = mk1 + mk2*nlevel + 1;
  % size(groupid)
  
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

  
  grm = grmeans*wt;
  meanvar = grvars*wt;
  varmean = (grmeans - repmat(grm,1,npages)).*...
	    (grmeans - repmat(grm,1,npages)) *wt;
  grsd = sqrt( varmean + meanvar );
  
  % make markers for the intplots
  if( cross == 'f2' )
    markers = [ 's' 'o' 'd' ];
    mksize = [ 6 8 9 ];
  elseif( cross =='bc' )
    markers = [ 'o' 'd' ];
    mksize = [ 6 9 ];
  end
  
  % plus minus
  pm = [ -2 +2 ];

  % handle to plots that will be legended initialized
  h = zeros(1,nlevel);
  % make the plots now
  for( j = 1:nlevel )

    % ids of the lines for a fixed level of marker 2
    id = (j-1)*nlevel + (1:nlevel);
    % plot the connecting means with the marker and marker size
    % specified; handle is assigned
    h(j) = plot( 0:(nlevel-1), grm(id), markers(j), 'markersize', ...
		 mksize(j) );

    hold on;
    % plot connecting lines between means
    plot( 0:(nlevel-1), grm(id) );

    for( i = 1:nlevel )
      % for each fixed level of marker 1
      % get the group id
      id = (i-1) + (j-1)*nlevel + 1;
      % plot the error bars
      plot( repmat( i-1, 1, 2 ), grm(id) + pm*grsd(id) );
      % plot whiskers of error bars
      plot( i-1 + pm/50, repmat( grm(id) + 2*grsd(id), 1, 2 ) );
      plot( i-1 + pm/50, repmat( grm(id) - 2*grsd(id), 1, 2 ) );
    end
  end
  
  % set the ticks on the x axis
  set( gca, 'xtick', 0:(nlevel-1) );
  xlim( [ -1/2 (nlevel-1/2) ] )
  % label ticks on x axis
  set( gca, 'xticklabel', labels );  
  % label the x axis
  xlabel( strcat( mk1name, ' genotype' ) );    
  % put legends
  legend( h, legendlabels, 0 )  
