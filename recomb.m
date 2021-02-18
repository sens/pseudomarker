function y=recomb(x)
% RECOMB Calculate the recombination fraction.
%
% RECOMB(X) calculates the recombination fraction corresponding to a
% genetic distance of X Morgans.  Works for vectors too.  

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/07/10 15:49:30 $	

  y = (1-exp(-2*abs(x)))/2;