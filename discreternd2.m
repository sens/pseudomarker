function i = discreternd2( w )
% DISCRETERND2 Simulates from a discrete distribution.
%
% DISCRETERND2(W) simulates from a discrete distribution whose weights
% are given by the matrix W.  The columns correspond to different randdom
% variables.   
%
% See also DISCRETERND.  

  % 
  w = cumsum(w); 
  s = w(end,:); % last row
  [nrows, ncols]= size(w);
  A = repmat(s,nrows,1);
  W = w./A;
  u = unifrnd(0,1,1,ncols);
  U = repmat(u,nrows,1);

  i = sum( (U-W) > 0 ) + 1;
