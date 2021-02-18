function obsdx = phenoobs( y )
% PHENOOBS Find the indexes of the rows for which no phenotype is missing
%
% OBSDX = PHENOOBS(Y)
% Y = matrix of phenptypes
% OBSDX = index of individuals who had no missing phenotype column  

  tr = size(y,2);
  if( tr==1 )
    obsdx = find( ~isnan( y ) );
  elseif( tr>1 )
    notmiss = (~isnan( y ));
    notmiss = prod( double(notmiss') )';
    obsdx = find( (notmiss == 1) );
  else
    error( 'Must have a trait vector.' );
  end