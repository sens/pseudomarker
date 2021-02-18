function wt = weight_binomial( y, x, basewt )
% WEIGHT_BINOMIAL Weight function for binomial distribution.
%
  
  wt = binowt( y, x ) - basewt;
  
