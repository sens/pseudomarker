function u = bcmodel2( gt )
% BCMODEL2 Model matrix for the backcross.
%
% BCMODEL2(GT) gives the model matrix corresponding to backcross
% genotypes which are assumed to be coded as 0,1.  The estimated effects
% are the allelic effects.

  u = [ 1-gt gt ];
