function wt = weight_normal( y, x, baserss )
% WEIGHT_NORMAL Weight function for normal distribution.
%
  ss = rss(y,x);
  wt = -(n/2) * ( log(ss) - log(baserss) );
