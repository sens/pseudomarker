function p = plotmissingprop( missprop )
% PLOTMISSINGPROP Plot the proportion of missing data
%
% PLOTMISSINGPROP(MISSPROP)  
% MISSPROP = output of the function MISSINGPROP
  
%  figure
  nchroms = length(missprop);
  mprop = missprop;

  hold on;
  ylim( [0 1] );
  prop=[];
  cumlen = 0;
  chromlen = zeros( 1, nchroms );
  chroms = zeros( 1, nchroms );
  for( i=1:nchroms )
    spacing = mprop(i).mpos(end) - mprop(i).mpos(end-1);
    thisprop = mprop(i).missprop;
    prop = [ prop thisprop ];
    thispos = cumlen + mprop(i).mpos + (i-1) * spacing;
    cumlen = cumlen + mprop(i).mpos(end);
    chromlen(i) = mprop(i).mpos(end);
    chroms(i) = mprop(i).chrid;
    
    plot( thispos, thisprop )
  end

  set(gca,'Ticklength',[0 0])
  labelpos = cumsum(chromlen+spacing);
  set(gca,'XTick',labelpos-chromlen/2);
  set(gca,'XTickLabel',chroms);  
  
  hold off

  
%  for(i=1:nchroms)
%    subplot(1,nchroms,i)
%    plot( mprop(i).mpos, mprop(i).missprop );
%    set(gca,'ylim',[0 1]);
%    set( gca, 'xtick', round(unique(mprop(i).mpos)*100)/100 );
%    set( gca, 'xticklabel', '' )
%  end