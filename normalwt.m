function wt = normalwt( y, X, b0, v0, n0, rss0 )
% NORMALWT Weight function for normal distribution.
%
% WT = NORMALWT( Y, X )  
% WT = NORMALWT( Y, X, B0, V0, N0, RSS0 )
%
% Y = response vector
% X = design matrix
% B0 = prior mean for parameters
% V0 = prior precision for parameters (vector)
% N0 = prior df for variancs
% RSS0 = prior SS for variance
%
% The defaults are:
%    N0=0 RSS0=0 
%    B0=[mean(Y) ZEROS(1,P-1)] 
%    V0=ONES(1,p)  

  [n,k] = size(y);
  [n,p] = size(X);
  %ss = rss(y,x);
  %wt = -(n/2) * log(ss) - (p*k/2)*log(n);

  if( nargin==2 )
    n0=0;
    rss0=0;
    b0 = zeros(1,p);
    b0(1) = mean(y);
    v0 = ones(1,p);
  end
  
  [b,rss,v] = bls(y,X,b0,v0);
  wt = marginal(n+n0,rss+rss0,v,diag(v0));