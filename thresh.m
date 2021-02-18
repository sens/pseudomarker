function t = thresh( maxlod, p )
% THRESH Get thresholds from output of permutation tests on the LOD scale
%
% T=THRESH(MAXLOD)
% T=THRESH(MAXLOD,P)  
%  
% T = vector of thresholds on the LOD scale, i.e. log base 10.
% MAXLOD = output from either PERMUTEST or PERMUTEST2; it is a matrix of
%          the maximum of LOD scores (on the natural log scale)
% P = vector of percent levels of significance desired; by default it is 
%     [ 10 5 1 ]
% 
% See also: PERMUTEST, PERMUTEST2.  
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/09/20 22:55:54 $	
  

  if( nargin == 1 )
    p = [ 10 5 1 ];
  end
  
  t = prctile( maxlod', 100-p );
  