function ww = bromanwt( y, x1, x2 )
% BROMANWT Calculate weights for Karl Broman's distribution on log scale.
%
% y = matrix of data whose first component is a 0/1 variable and the
%     second component is a lifetime.  The second component exists only
%     if the first component is 1.
%  
% x1 = grouping variable for the binomial part.
% x2 = model matrix for the normal part.
  
  
  w1 = binowt( y(:,1), x1 );
  w2 = gausswt( y(y(:,1)==1,2), x2(y(:,1)==1,:));
  ww = w1+w2;