function p = marginal(df,rss,v,v0)
% MARGINAL Calculate marginal distribution.
%
% df = posterior degrees of freedom for variancs
% rss = posterior sum of squares
% v = posterior precision matrix
% v0 = prior precision matrix
  
  p =  -0.5*df*log(rss) - 0.5 * log( det(v) / det(v0) );