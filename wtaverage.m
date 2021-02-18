function m = wtaverage( y )
% WTAVERAGE Average the weights; output and input on log scale.
% 
% Change this function to whatever you want the averaging function to
% be.  Examples would include LOGNORMALMEAN and MEAN on the exponential
% scale. 
%  
  
  m = lognormalmean( y );
