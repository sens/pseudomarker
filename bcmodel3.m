function u = bcmodel3( gt )
% BCMODEL3 Model matrix for the backcross.
%
% BCMODEL3(GT) gives the model matrix corresponding to backcross
% genotypes which are assumed to be coded as 0,1.  This coding gives
% the following effect estimate:
%      Effect = (level2-level1),
%               the effect of a single allelic substitution
%  
  u = gt-0.5;    
