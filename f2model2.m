function u = f2model2( gt )
% F2MODEL2 Model matrix for the intercross.
%
% F2MODEL2(GT) gives the model matrix corresponding to intercross
% genotypes which are assumed to be coded as 0,1,2.  This coding gives
% estimated effects equal to the allelic effects.
%  
  u = [ (gt==0) (gt==1) (gt==2) ];  
