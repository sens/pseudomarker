function plotmainestimate( onelod, type )
% PLOTMAINESTIMATE Plots the output from MAINESTIMATE
%
% PLOTMAINESTIMATE(LOD) plots the output LOD from MAINESTIMATE.
% There are two panels of this plot.  The upper panel plots the
% proportion of variance explained as a function of genome location.  The
% bottom panel plots the estimated effects in each chromosome along with
% the error bars corresponding to approximate 95% confidence intervals.
%  
% PLOTMAINESTIMATE(LOD,TYPE) plots the output LOD from MAINESTIMATE, but
% instead of plotting the proportion of variance explained, it plots the
% LOD score if TYPE='lod', the LOD score in natural logarithms if
% TYPE='loglik' and the proportion of variance explained if
% TYPE='propvar'.
%
% See also MAINESTIMATE.  

  if( nargin==1)
    type='propvar';
  end
  
  n = onelod(1).n;
  cross = onelod(1).cross;
%  figure   
  subplot( 2, 1, 1 );

  if( nargin==3 )
  %  perm = s
  end
  
  hold on;
  lod = [];
  nchroms = length( onelod );
  cumlen = 0;
  chromlen = zeros( 1, nchroms );
  chroms = zeros( 1, nchroms );
  for( i=1:nchroms )
    spacing = onelod(i).mpos(end) - onelod(i).mpos(end-1);
    thislod = onelod(i).lod;
    lod = [ lod thislod ];
    thispos = cumlen + onelod(i).mpos + (i-1) * spacing;
    cumlen = cumlen + onelod(i).chromlen;
    chromlen(i) = onelod(i).chromlen;
    chroms(i) = onelod(i).chrid;
    
    switch type
     case 'lod'
      plot( thispos, thislod/log(10) );           
     case 'propvar'
      plot( thispos, 1-exp(-2*thislod/n) );     
     case 'loglik'
      plot( thispos, thislod );     
    end

  end

  set(gca,'Ticklength',[0 0])
  labelpos = cumsum(chromlen+spacing);
  set(gca,'XTick',labelpos-chromlen/2);
  set(gca,'XTickLabel',chroms);  
  
  hold off
  
  
  subplot( 2, 1, 2 );
  if( cross=='bc' )
    effects = zeros(1,1);
    effvar = zeros(1,1);      
  elseif( cross=='f2' )
    effects = zeros(1,2);
    effvar = zeros(1,2);      
  end
  
  effectpos = 0;
  hold on
  for( i=1:nchroms )
    if( cross=='bc')
      effects=onelod(i).effects(2);
      effvar=onelod(i).effvar(2);    
    elseif( cross=='f2' )
      effects=onelod(i).effects(2:3);
      effvar=onelod(i).effvar(2:3);    
    end

    m = length(effects);
    effectpos = (2*i*(m*2))-1;
    for( j=1:m )
      plot(j+effectpos,effects(j),'.');
      plot([j j]+effectpos, [ effects(j) + 2*sqrt(effvar(j)) * [-1 1] ] );
    end
    
  end
  

  labelpos = 4*(1:nchroms)*m+(m-1)/2;
  set(gca,'XTick',labelpos);
  set(gca,'XTickLabel',chroms);  
  xlim( [0 4*m*nchroms+3] )
  refline(0,0);  
  
  