function wt = gausswt( y, X )
% NORMALWT Weight function for normal distribution.
%
% WT = GAUSSWT( Y, X )  
%
% Y = response vector
% X = design matrix


  [n,k] = size(y);
  [n,p] = size(X);
  ss = rss(y,X);
  wt = -(n/2) * log(ss);
