function h=plotpairscan( twolod, type )
% PLOTPAIRSCAN Plots a two-dimensional genome scan  
%
% PLOTPAIRSCAN(TWOLOD) plots the results of a two-dimensional genome scan
% stored in TWOLOD.  Typically TWOLOD will be the output of a PAIRSCAN.
% The graph will be in the "proportion of variance" scale.  By using the
% PLOTPAIRSCAN(TWOLOD,TYPE) command, the scale of the plot can be
% controlled.  TYPE takes the following options:
%      'prop': proportion of variance
%      'lod': LOD scale, log base 10; default
%      'llik': log likelihood ratio scale, base e
%
% See also: PAIRSCAN, PLOTMAINSCAN.  
%
%   

  n = twolod(1).n;
  cross = twolod(1).cross;

  if( nargin==1 )
    type = 'lod';
  end
  
  
  switch cross
   case 'bc'
    dfmultiplier = 3;
   case 'f2'
    dfmultiplier = 2;
  end
  
    lod = [];
  
  [ a b ] = size( twolod );
  nchroms = a;
 
  % position vector for plot
  pos = [];
  % number of pseudomarkers
  npseudo = [];
  chroms = zeros(1,nchroms);
  spacing = zeros(1,nchroms);
  
  for( i=1:nchroms )
    chroms(i) = twolod(i,i).chrid1;
    chromlen(i) = twolod(i,i).chromlen1;
    npseudo = [ npseudo length(twolod(i,i).mpos1) ] ;
  end

    
  % set size of the lod matrix
  lod = zeros( sum(npseudo), sum(npseudo) );
  % how many pseudomarkers in to the left of each chromosome + 1
  npseudo = 1 + [ 0 cumsum(npseudo) ];
  
  % baselod
  basex = 0;
  
  for( i=1:nchroms )
    for( j=i:nchroms )
      % [i j ]
      if(i==j)
	% thislod = twolod( i, i ).lod;	
	switch type
	 case 'llik'
	  thislod = twolod( i, i ).lod - basex ;
	 case 'prop'
	  thislod = 1-exp( - (twolod( i, i ).lod - basex )*2/n);
	 case 'lod';
	  thislod = (twolod( i, i ).lod - basex )/log(10);
	end
     	lod( npseudo(i):(npseudo(i+1)-1), npseudo(i):(npseudo(i+1)-1) ) ...
	    = triu( thislod ) + dfmultiplier * ...
	    ( triu( thislod )' - tril( thislod ) ); 
	
      else
	switch type
	 case 'llik'
	  thislod = twolod( i, j ).lod-basex;
	  thatlod = twolod( j, i ).lod-basex;	
	 case 'prop'
	  thislod = 1-exp( -2 * (twolod( i, j ).lod-basex) / n);
	  thatlod = 1-exp( -2 * (twolod( j, i ).lod-basex) / n);	
	 case 'lod';
	  thislod = (twolod( i, j ).lod-basex)/log(10);
	  thatlod = (twolod( j, i ).lod-basex)/log(10);	
	end
	lod( npseudo(i):(npseudo(i+1)-1), npseudo(j):(npseudo(j+1)-1) ) ...
	    = thislod;
        lod( npseudo(j):(npseudo(j+1)-1), npseudo(i):(npseudo(i+1)-1) ) ...
	    = dfmultiplier*(thislod'-thatlod);
      end
      
	   
    end

    spacing(i) = twolod(i,i).mpos1(end) - twolod(i,i).mpos1(end-1);
    thispos = sum(chromlen(1:(i-1)))+ twolod(i,i).mpos1 + (i-1)*spacing(i); 
    pos = [ pos thispos ];
  end
  

  tmppos = [pos pos(end)+pos(2)-pos(1)];
  tmplod = [ lod zeros(npseudo(end)-1,1); zeros(1,npseudo(end)) ];
  axis square

  h=surf( tmppos, tmppos, tmplod );
  shading flat; 
  mmin = min( tmplod( find( tmplod~=0 ) ) );
  mmax = max( tmplod( find( tmplod~=0 ) ) );  
  caxis( [ mmin mmax ] );
  view(2);

  set(gca,'Xlim',[0 max(tmppos)]);  
  set(gca,'Ylim',[0 max(tmppos)]);
  set(gca,'TickDir','out');  
  set(gca,'TickLength',[ 0.01 0.01 ] );  
  labelpos = cumsum(chromlen+spacing)-spacing/2; 
  set(gca,'XTick',labelpos);
  set(gca,'XTickLabel',' ');  
  clabelpos = labelpos - chromlen/2;
  ml = max(labelpos);
  for( i=1:nchroms )
    text( 'position', [ clabelpos(i), -ml/40 ], ...
	  'string', int2str(chroms(i)), ...
	  'horizontalalignment', 'center' );
  end
  
  set(gca,'YAxisLocation','left');
  set(gca,'YTick',labelpos);
  set(gca,'YTickLabel','   ');  
  clabelpos = labelpos - chromlen/2;
  for( i=1:nchroms )
    text( 'position', [ -ml/25 clabelpos(i) ], ...
	  'string', int2str(chroms(i)) );
  end
    
  cl = colormap;
  colormap( cl( 7:end,: ) );


  h1=colorbar;

  
