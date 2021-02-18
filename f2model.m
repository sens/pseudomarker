function u = f2model( gt )
% F2MODEL Model matrix for the intercross.
%
% F2MODEL(GT) gives the model matrix corresponding to intercross
% genotypes which are assumed to be coded as 0,1,2.  This coding gives
% the following effect estimates:
%      Effect1 = (level2-level1)/2,
%                the additive effect of a single allelic substitution
%      Effect2 = 0.5 * ( (level2-level1) - (level1-level0) )
%                the dominance effect or the quadratic effect
%  
  v = gt-1;
  z = v.*v - 2/3;
  u = [ v z ];  

