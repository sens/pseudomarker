function f=lod2f( l, n, df )
% LOD2F Convert log likelihood to proportion of variance
%                explained
%
% F=LOD2F(LOD,N,DF) returns the F value corresponding to a LOD score in
% when the numerator degrees of freedom is DF.  This is an approximate
% relationship. 

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2001/09/21 19:28:40 $	

  f = (n-df)*(exp(2*l*log(10)/n)-1)/df;
  