function [g,mnames,mpos] = getmarkergeno( m, pos )
% GETMARKERGENO Get marker genotypes for particular position on a chromosome   
%
% G=GETMARKERGENO(M,MPOS)
% M = marker genotype structure  
% POS = approximate marker positions

  if( isfield( m(1), 'geno' ) )
    imputed = 0;
  elseif ( isfield( m(1), 'igeno' ) )
    imputed = 1;
  else
    error( 'This is not a marker data structure.' );
  end
  
  nchroms = length(m);
  if( imputed == 0 )
    n = size(m(1).geno,1);
  else
    n = size(m(1).igeno,1);
  end
  
  g = zeros(n,nchroms);
  mnames = [];

 
  for( i = 1:nchroms )
    [ mdist, idx ] = min( abs( m(i).mpos - pos(i) ) );
    if( imputed==0 )
      g(:,i) = m(i).geno(:,idx);
      mnames = [ mnames m(i).mnames(idx) ];
    else
      g(:,i) = m(i).igeno(:,idx,1);
      thischr = int2str( m(i).chrid );
      thismk = int2str( idx );
      mnames = [ mnames  strcat( thischr, '-', thismk )];
    end
    mpos(i) = m(i).mpos(idx);
  end

