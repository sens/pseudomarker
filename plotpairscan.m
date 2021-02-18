function h=plotpairscan( twolod, varargin )
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
%      'F': F statistic scale
% PLOTPAIRSCAN(TWOLOD,'type',TYPE,'onelod',ONELOD) 
% will plot the pairscan of type TYPE and will plot the lower
% triangle as LOD(Q1*Q2)-max(LOD(Q1),LOD(Q2)).  
% See also: PAIRSCAN, PLOTMAINSCAN.  

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.837 $ $Date: 2001/09/24 23:21:04 $	
   

  n = twolod(1).n;
  cross = twolod(1).cross;

  if( length(varargin) > 0 )
    if( length(varargin) == 1 )
      opts = struct('type',varargin{1});
    else
      opts = struct(varargin{:});
    end
  else 
    opts = [];
  end
  
    
  if( ~isfield(opts,'type') )
    type = 'lod';
  else 
    type = opts.type;
  end
  if( ~isfield(opts,'onelod') )
    onelod = [];
  else
    onelod = opts.onelod;
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
  nmarks = zeros(1,nchroms);  
  
  % are the pseudomarkers/markers equally spaced?
  rangediff = 0;
  for( i=1:nchroms )
    chroms(i) = twolod(i,i).chrid1;
    chromlen(i) = twolod(i,i).chromlen1;
    nmarks(i) = length(twolod(i,i).mpos1);
    npseudo = [ npseudo length(twolod(i,i).mpos1) ];
    rangediff = max( rangediff, range(diff(twolod(i,i).mpos1)) );
  end

    
  % set size of the lod matrix
  lod = zeros( sum(npseudo), sum(npseudo) );
  % how many pseudomarkers in to the left of each chromosome + 1
  npseudo = 1 + [ 0 cumsum(npseudo) ];
  
  % baselod
  basex = 0;
  
  switch cross
   case 'bc'
    dffull = 3;
    dfint = 1;
   case 'f2'
    dffull = 8;
    dfint = 4;
  end
  
  for( i=1:nchroms )
    for( j=i:nchroms )
      % [i j ]
      if(i==j)
	% thislod = twolod( i, i ).lod;	
	switch type
	 case 'llik'
	  thislod = (twolod( i, i ).lod - basex )*log(10);
	 case 'prop'
	  thislod = 1-exp( - (twolod( i, i ).lod - basex )*log(10)*2/n);
	 case 'lod';
	  thislod = (twolod( i, i ).lod - basex );
	 case 'F';
	  thislod = twolod( i, i ).lod - basex;
	 otherwise
	  error( 'Unknown type of plot.' )
	end
  
	% if we have F statistics we have to do something different
	if( strcmp( type, 'F' ) )
	  lod( npseudo(i):(npseudo(i+1)-1), npseudo(i):(npseudo(i+1)-1) ) ...
	      = lod2f( triu( thislod ), n, dffull ) + ...
	        lod2f( triu( thislod )' - tril( thislod ), n, dfint ); 
	else
	  lod( npseudo(i):(npseudo(i+1)-1), npseudo(i):(npseudo(i+1)-1) ) ...
	      = triu( thislod ) + dfmultiplier * ...
	      ( triu( thislod )' - tril( thislod ) ); 
	end
	
      else
	switch type
	 case 'llik'
	  thislod = log(10) * (twolod( i, j ).lod-basex);
	  thatlod = log(10) * (twolod( j, i ).lod-basex);	
	 case 'prop'
	  thislod = 1-exp( -2 * log(10) * (twolod( i, j ).lod-basex) / n);  
	  thatlod = 1-exp( -2 * log(10) * (twolod( j, i ).lod-basex) / n);
	 case 'lod';
	  thislod = (twolod( i, j ).lod-basex);
	  thatlod = (twolod( j, i ).lod-basex);	
	 case 'F';
	  thislod = twolod( i, j ).lod - basex;
	  thatlod = twolod( j, i ).lod - basex;	  
	end

	% if we have F statistics we have to do something different
	if( strcmp( type, 'F' ) )
	  lod( npseudo(i):(npseudo(i+1)-1), npseudo(j):(npseudo(j+1)-1) ) ...
	      = lod2f( thislod, n, dffull );
	  lod( npseudo(j):(npseudo(j+1)-1), npseudo(i):(npseudo(i+1)-1) ) ...
	      = lod2f( thislod'-thatlod, n, dfint );
	else
	  lod( npseudo(i):(npseudo(i+1)-1), npseudo(j):(npseudo(j+1)-1) ) ...
	      = thislod;
	  lod( npseudo(j):(npseudo(j+1)-1), npseudo(i):(npseudo(i+1)-1) ) ...
	      = dfmultiplier*(thislod'-thatlod);
	end
      end
      
	   
    end

    spacing(i) = ( twolod(i,i).mpos1(end) - twolod(i,i).mpos1(1) ) ...
	/ (nmarks(i)-1);
    if( rangediff < 1e-5 )
      thispos = sum(chromlen(1:(i-1)))+ twolod(i,i).mpos1 + ...
		(i-1)*spacing(i); 
      pos = [ pos thispos ];
    else
      spacing = 1;
      thispos = sum(nmarks(1:(i-1))) + (1:nmarks(i));
      pos = [ pos thispos ];
    end
  end

  tmppos = [pos pos(end)+pos(2)-pos(1)];
  %%%%%%
  %% NEW begin
  if( ~isempty(onelod) )
    mainlod = [];
    i=1;
    while( i<= nchroms )
      mainlod = [ mainlod onelod(i).lod ];
      i=i+1;
    end
    % lod = lod -diag(diag(lod)+alod');
    alod = repmat(mainlod,size(lod,1),1);
    blod = alod';
    mlod = zeros(size(lod,1),size(lod,1),2);
    mlod(:,:,1)=alod;
    mlod(:,:,2)=blod;  
    mlod = max(mlod,[],3);
    mlod = triu(mlod);
    lod = lod - mlod;
    lod = lod - diag(diag(lod)) + diag(mainlod);
  end
  %%% NEW end
  tmplod = [ lod zeros(npseudo(end)-1,1); zeros(1,npseudo(end)) ];
  axis square;

  h=surf( tmppos, tmppos, tmplod );
  shading flat; 
  mmin = min( tmplod( find( tmplod~=0 ) ) );
  mmax = max( tmplod( find( tmplod~=0 ) ) );  
  caxis( [ mmin mmax ] );
  view(2);

  set(gca,'Xlim',[min(tmppos) max(tmppos)]);  
  set(gca,'Ylim',[min(tmppos) max(tmppos)]);
  set(gca,'TickDir','out');  
  set(gca,'TickLength',[ 0.01 0.01 ] );  
  if( rangediff < 1e-5 )
    labelpos = cumsum(chromlen+spacing)-spacing/2;
  else 
    labelpos = cumsum(nmarks)+spacing/2;
    chromlen = (nmarks-1);
  end
  
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
  xlabel( 'Chromosome number' );
  ylabel( 'Chromosome numner' );

  hh=colorbar;
  hhh = get( hh, 'xlabel' );
  set( hh, 'xaxislocation', 'top' );
  switch type
   case 'lod'
    set( hhh, 'string', 'LOD score' );
    case 'prop'    
    set( hhh, 'string', 'Proportion of variance' );
   case 'llik'
    set( hhh, 'string', 'Log likelihood' );
   case 'F'
    set( hhh, 'string', 'F' );
  end
  
