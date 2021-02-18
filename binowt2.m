function ww = binowt2( y, x )
% BINOWT2 Calculate binomial weights on the log scale.  

%  persistent BASE;
  
  u = unique(x);
  k = length(u);
  
  n = length(y);

%  if( exist( 'BASE' )==0 )
    BASE = zeros(n,k);
 % end
  
  for( i=1:k )
    BASE(:,i) = (x==u(i));
  end
  
  m = sum(BASE); % number in each group
  z = sum( prod_model( y, BASE ) ); % number of 1's in each group

  bbb = betaln( z+1, m-z+1 );
  ww = sum(bbb);