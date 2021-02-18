function vv = vvar( x, p )
% VVAR Estimate variance 

  sz = size( x );
  pp = repmat( p, 1, sz(2) );
  mm = sum( x.*pp );
  md = sum( (abs(x-repmat(mm,sz(1),1))).*pp );
  vv = md.^2 / (2/pi);