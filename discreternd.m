function i = discreternd( w, n )
% DISCRETERND Simulates from a discrete distribution.
%
% DISCRETERND(W,N) simulates N samples from a discrtete distribution whose
% probability distribution is given by the vector W.  The vector W may be
% un-normalized.  The output is the index corresponding to the weight
% vector.
%
% See also DISCRETERND2.
  
  w = w/sum(w);
  w = cumsum(w);
  u = unifrnd(0,1,1,n);
  
  [U,W] = meshgrid(u,w);
  i = sum( (U-W) > 0 ) + 1;
