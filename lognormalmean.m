function m = lognormalmean( x )
% LOGNORMALMEAN Estimate mean of lognormal distribution.

  sz = size( x );

  if( length(sz) > 2 )
    xx = reshape( x, sz(1), prod( sz(2:end) ) );
    %    mm = mean( xx );
    %    vv = var( xx );
    [mm,vv] = trimmedmean( xx );
    mmm = exp( mm + 0.5 * vv );
%    [mmm,vv] = trimmedmean( exp(xx) );
    m = reshape( mmm, sz(2:end) );
  else
    %mm = mean( x );
    %vv = var( x );
    [mm,vv] = trimmedmean( x );
    m = exp( mm + 0.5 * vv );
  end

%  m=exp( median(x) );