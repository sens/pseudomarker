function p=loglik2propvar( l, n )
% LOGLIK2PROPVAR Convert log likelihood to proportion of variance
%                explained
%
% P=LOGLIK2PROPVAR(L,N)
% Returns the proportion of variance explained corresponding to a
% loglikelihood of L when the sample size is N.  This is an approximate
% relationship.

p = 1-exp(-2*l/n);     