function megalod = imp_mainscan( y, x, qtldata, spacing, npages, cross, ...
				chroms, notmiss )
% MAINSCAN2 Plot and compute a one-dimensional genome scan.
%
% y = phenotypes; can be a matrix
% x = covariates, use [] if none  
% qtldata = qtl data structure  
% notmiss = not missing index
% chroms = which chromosomes, use [] for all
%  
% onelod = lod score structure
%
  
  y = y(notmiss,:);
  if( ~isempty(x) )
    x = x(notmiss,:);
  end
  
  [ n, ntr ] = size(y);

  nchroms = length( qtldata );
  if( isempty(chroms) )
    chroms = 1:nchroms
  else
    nchroms = length(chroms);
  end
  
  fake = pseudo( qtldata, spacing, npages, cross, chroms );
  onelod = mainscan2( y, x, qtldata, fake, notmiss );
  plotmainscan( onelod );
  megalod = onelod;
  aaa = input( 'More iterations? (y/n)', 's' );
  i = 1;
  while( strcmp( aaa, 'y' ) == 1 )
    i = i + 1;
    fake = pseudo( qtldata, spacing, npages, cross, chroms );
    onelod = mainscan2( y, x, qtldata, fake, notmiss );
    plotmainscan( onelod, 'b' );
    
    for( j=1:length(onelod) )
      megalod(j).lod = (1-1/i)*megalod(j).lod + (1/i)*onelod(j).lod;
      megalod(j).bf  = (1-1/i)*megalod(j).bf  + (1/i)*onelod(j).bf;
    end
    plotmainscan( megalod, 'r' );
    aaa = input( 'More iterations? (y/n)', 's' );
  end
  
  