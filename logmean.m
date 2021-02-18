function lmean = logmean( x )
% LOGMEAN Calculate mean on the log scale
%
% X is a vector or matrix  
%
  mmm = max(x);
  [nr,nc] = size(x);
  
  if(nc==1)
    z = x - mmm;
  else
    z = x - repmat(mmm,nr,1);
  end
  
  avg = mean(exp(z));
  lmean = mmm + log(avg);