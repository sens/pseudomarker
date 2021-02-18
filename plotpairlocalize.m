function h=plotpairlocalize( twolod, opt )
% PLOTPAIRLOCALIZE Plots the posterior distribution using a two-QTL model
% on a single chromosome.
%
% PLOTPAIRLOCALIZE( TWOLOD, LEVELS )
% TWOLOD = lod output from PAIRSCAN from the chromosome of interest
% LEVELS = levels of the confidence regions shaded; default is
%          [ 0.5 0.95 0.99 ]
%
% See also: PAIRSCAN.
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2002/01/11 19:41:27 $	

  thislod = twolod.lod*log(10);
  spacing = mean(diff(twolod.mpos1));
  maximum = max( thislod(:) );
  post = triu( exp( thislod - maximum ) );
  const = sum( post(find(~isnan(post(:)))) );
  post = post/const;
  M=max(post(:));
  npseudo = length( twolod.mpos1 );

  %  size( post )
  %  npseudo
  
  if( nargin==1 )
    opt = [ 0.5 0.95 0.99 ];
  end

  pos = [twolod.mpos1 twolod.mpos1(end)+ spacing ];
  [p,i]=sort(post(:));
  p(i) = cumsum(p);
  ppp = reshape( p, size( post ) );
  post = ppp;
  post = [ post zeros(npseudo,1); zeros(1,npseudo+1) ] ;
  
  map = zeros( 1000, 3 );
  nlev = length(opt);
  
  opt = [ 1 round(opt*1000) 1000 ];
  
  for( i=1:(nlev+1) )
    map( opt(i):opt(i+1), : ) = repmat( [ i-1 i-1 i-1 ]/nlev, ...
					opt(i+1)-opt(i)+1, 1 );	
  end
  map = map( 1000:-1:1,:);
  
  h=surf( 100*pos, 100*pos, post );
  shading flat;
  grid off;
  xlim( [ 0 100*max(pos) ] );
  ylim( [ 0 100*max(pos) ] );  
  view(2);
  
  colormap(map);
  

