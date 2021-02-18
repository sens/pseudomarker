function p = plotsegdist( missprop )
% PLOTSEGDIST Plot the segregation distortion scan
%
% PLOTSEGDIST(FAKE)  
% FAKE = imputed data structure
  
  nchroms = length(missprop);
  mprop = missprop;

  hold on;
%  ylim( [-1 1] );
  prop=[];
  cumlen = 0;
  chromlen = zeros( 1, nchroms );
  chroms = zeros( 1, nchroms );
  for( i=1:nchroms )
    spacing = mprop(i).mpos(end) - mprop(i).mpos(end-1);
    thisprop = mprop(i).segdist;
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
  
  refline( 0, 0 );
  hold off

  
%  for(i=1:nchroms)
%    subplot(1,nchroms,i)
%    plot( mprop(i).mpos, mprop(i).missprop );
%    set(gca,'ylim',[0 1]);
%    set( gca, 'xtick', round(unique(mprop(i).mpos)*100)/100 );
%    set( gca, 'xticklabel', '' )
%  end