function [lod,bf] = hkonescan( y, x, geno, mpos, qpos, cross )
% HKONESCAN Cmpute the approximate posterior distribution of QTL assuming
% there is only one of them. Uses the Haley-Knott method.
% 
% y = vector of trait values  
% x = matrix of covariates, use [] if none
% geno = marker genotypes; can be missing
% mpos = marker positions
% qpos = positions where the LOD score has to be calculated  
% cross = cross type, 'bc' for BACKCROSS and 'f2' for INTERCROSS
%         DOES NOT WORK FOR INTERCROSS YET
  
if( cross == 'bc' )
  [lod,bf] = hkonescan_bc( y, x, geno, mpos, qpos );
elseif( cross == 'f2' )
  [lod,bf] = hkonescan_f2( y, x, geno, mpos, qpos );
else
  error( 'Cross type not recognized.' );
end

