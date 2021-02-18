function h=plotpairlocalize2( twolod, opt )
% PLOTPAIRLOCALIZE Plots the posterior distribution using a two-QTL model
% on a two different chromosomes.
%
% PLOTPAIRLOCALIZE2( TWOLOD, LEVELS )
% TWOLOD = lod output from the chromosomes of interest
% LEVELS = levels of the confidence regions to be shown; default is
%          [ 0.5 0.95 0.99 ]

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2002/01/11 19:40:53 $	
  
  thislod = twolod.lod*log(10);
  spacing1 = mean(diff(twolod.mpos1));
  spacing2 = mean(diff(twolod.mpos2));  
  maximum = max( thislod(:) );
  post = ( exp( thislod - maximum ) );
  const = sum( post(find(~isnan(post(:)))) );
  %post = 2 * post / ( spacing * spacing * const );
  post = post/const;
  M=max(post(:));
  npseudo1 = length( twolod.mpos1 );
  npseudo2 = length( twolod.mpos2 );

%  size( post )
%  npseudo

  if( nargin==1 )
    opt = [ 0.5 0.95 0.99 ];
  end    

  pos1 = [twolod.mpos1 twolod.mpos1(end)+ spacing1 ];
  pos2 = [twolod.mpos2 twolod.mpos2(end)+ spacing2 ];  
  %figure;
  
  [p,i]=sort(post(:));
  p(i) = cumsum(p);
  ppp = reshape( p, size( post ) );
  post = ppp;
  %post = log(ppp)/log(2);
  post = [ post zeros(npseudo1,1); zeros(1,npseudo2+1) ] ;
  
  map = zeros( 1000, 3 );
  nlev = length(opt);
  
  opt = [ 1 round(opt*1000) 1000 ];
  
  for( i=1:(nlev+1) )
    map( opt(i):opt(i+1), : ) = repmat( [ i-1 i-1 i-1 ]/nlev, ...
					opt(i+1)-opt(i)+1, 1 );	
  end
  map = map( 1000:-1:1,:);
  
  h=surf( 100*pos2, 100*pos1, post );
  shading flat;
  grid off;
  xlim( [ 0 100*max(pos2) ] );
  ylim( [ 0 100*max(pos1) ] );  
  view(2);
  
  colormap(map);
  


