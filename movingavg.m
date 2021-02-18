function mm = movingavg( x, w )
% MOVINGAVG Calculate moving average of w observations

  n = size( x, 1 );
  
  mat = triu( ones(n) ) - triu( ones(n), w );
  
  if(w>1)
    mat( (n-w+2):end, : ) = [];
  end
    
  size(x)
  size(mat)
  mm = mat*x/w;