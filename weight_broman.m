function wt = weight_broman( y, x, basewt )
% WEIGHT_NORMAL Weight function for Karl Broman's data.
%
  ss = rss(y,x);
  wt = -(n/2) * ( log(ss) - log(baserss) );
