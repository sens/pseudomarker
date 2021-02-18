function cross = guesscross( m )
% GUESSCROSS Guess cross type by looking at marker data structure
%
% CROSS = GUESSCROSS(M)
% CROSS = cross type; 'bc' or'f2'
% M = marker data structure
%   

  if( isfield( m(1), 'geno' ) )
    ggg = m(1).geno;
  elseif ( isfield( m(1), 'igeno' ) )
    ggg = m(1).igeno;
  else
    error( 'This is not a marker data structure.' );
  end
  
  maxgeno = max(ggg(:));
  
  if( maxgeno == 1 )
    cross = 'bc';
  elseif( maxgeno == 2 )
    cross = 'f2';
  else
    error('Cannot determine cross type.' );
  end
  