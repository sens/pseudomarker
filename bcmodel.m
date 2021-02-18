function u = bcmodel( gt )
% BCMODEL Model matrix for the backcross.
%
% BCMODEL(GT) gives the model matrix corresponding to backcross
% genotypes which are assumed to be coded as 0,1.  This coding gives
% the following effect estimate:
%      Effect = (level2-level1)/2,
%               the effect of a single allelic substitution
%  
  %u = 2*gt-1;
  u = gt+gt-1;    
