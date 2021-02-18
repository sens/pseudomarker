function [mmm,vvv] = unconditional( m, v, p )
% UNCOMDITIONAL Find unconditional mean and variance based on conditional
%               means and variances and the conditional probabilities.  
%
% m = vevtor of means
% v = vectot of variances
% p = vector of probabilities
%

  k = size(m,2);
  
  mmm = sum( m.* repmat( p, 1, k ) ); % unconditional mean
  vvv = sum( v.* repmat( p, 1, k ) ) + vvar( m, p ); % unconditional variance
  
