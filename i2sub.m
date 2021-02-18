function sub = i2sub( qtl, ind )
% I2SUB Convert global index to subindex
%

% chrid nqtl nmarkers
  
  nchroms = size( qtl,1 );
  nindex = length(ind);
  
  npositions = qtl(1,3);
  nqtl = qtl(1,2);

  tmp = zeros( npositions, nindex );
  nposs = zeros( npositions, 1 );
  
  for( i=nqtl:-1:1 )
    % get the number of positions
    for( j=1:npositions )
      nposs(j) = orderedpositions(j,i)
    end
    
function nposs = orderedpositions(nmark,ntuple)
  nposs = 1;
  for( i=1:ntuple )
    nposs = nposs*(nmark-i+1)
  end
  