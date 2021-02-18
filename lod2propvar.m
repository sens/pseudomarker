function p=lod2propvar( l, n )
% LOD2PROPVAR Convert LOD score to proportion of variance explained
%
% P=LOGLIK2PROPVAR(L,N)
% Returns the proportion of variance explained corresponding to a
% loglikelihood of L when the sample size is N.  This is an approximate
% relationship.

p = 1-exp(-2*l*log(10)/n);     