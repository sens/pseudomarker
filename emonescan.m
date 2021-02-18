function [lod,bf] = emonescan( y, x, geno, mpos, qpos, cross )
% EMONESCAN Cmpute the approximate posterior distribution of QTL assuming
% there is only one of them. Uses the EM algorithm
% 
% y = vector of trait values  
% x = matrix of covariates, use [] if none
% geno = marker genotypes; can be missing
% mpos = marker positions
% qpos = positions where the LOD score has to be calculated  
% cross = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS
%         DOES NOT WORK FOR INTERCROSS YET
  
if( cross == 'bc' )
  [lod,bf] = emonescan_bc( y, x, geno, mpos, qpos );
elseif( cross == 'f2' )
  [lod,bf] = emonescan_f2( y, x, geno, mpos, qpos );
else
  error( 'Cross type not recognized.' );
end

