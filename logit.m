function logodds=logit(p)
% LOGIT Calculate logit
% 
% L = LOGIT(X)
%
  logodds=log(p./(1-p));