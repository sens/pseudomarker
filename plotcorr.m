function h=plotcorr( theta )
% PLOTCORR Plots correlations between pseudomarkers
%
% See also: PAIRSCAN.
%
%   

  n = theta(1).n;
  cross = theta(1).cross;

%  switch cross
%   case 'bc'
%    dfmultiplier = 3;
%   case 'f2'
%    dfmultiplier = 2;
%  end
  
    cor = [];
  
  [ a b ] = size( theta );
  nchroms = a;
 
  % position vector for plot
  pos = [];
  % number of pseudomarkers
  npseudo = [];
  chroms = zeros(1,nchroms);
  spacing = zeros(1,nchroms);
  
  for( i=1:nchroms )
    chroms(i) = theta(i,i).chrid1;
    chromlen(i) = theta(i,i).chromlen1;    
    npseudo = [ npseudo length(theta(i,i).mpos1) ] ;
  end

    
  % set size of the cor matrix
  cor = zeros( sum(npseudo), sum(npseudo) );
  % how many pseudomarkers in to the left of each chromosome + 1
  npseudo = 1 + [ 0 cumsum(npseudo) ];
  
  for( i=1:nchroms )
    for( j=i:nchroms )
      % [i j ]
      if(i==j)
	thiscor = theta( i, i ).cor;	
	cor( npseudo(i):(npseudo(i+1)-1), npseudo(i):(npseudo(i+1)-1) ) = ...
	    thiscor;
      else
	thiscor = theta( i, j ).cor;
	cor( npseudo(i):(npseudo(i+1)-1), npseudo(j):(npseudo(j+1)-1) ) ...
	    = thiscor;
        cor( npseudo(j):(npseudo(j+1)-1), npseudo(i):(npseudo(i+1)-1) ) ...
	    = thiscor';
	
      end
      
	   
    end

    spacing(i) = theta(i,i).mpos1(end) - theta(i,i).mpos1(end-1);
    thispos = sum(chromlen(1:(i-1)))+ theta(i,i).mpos1 + (i-1)*spacing(i); 
    pos = [ pos thispos ];
  end
  

  tmppos = [pos pos(end)+pos(2)-pos(1)];
  tmpcor = [ cor zeros(npseudo(end)-1,1); zeros(1,npseudo(end)) ];
  axis square

  h=surf( tmppos, tmppos, tmpcor );
  shading flat; 
  mmin = min( tmpcor( find( tmpcor~=0 ) ) );
  mmax = max( tmpcor( find( tmpcor~=0 ) ) );  
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

  
