function plotmainscan( onelod, type, varargin )
% PLOTMAINSCAN Plots a one-dimensional genome scan  
% 
% PLOTMAINSCAN(ONELOD) plots the results of a one-dimensional genome scan
% stored in ONELOD.  Typically ONELOD will be the output of a MAINSCAN.
% The graph will be in the "proportion of variance" scale.  By using the
% PLOTMAINSCAN(ONELOD,TYPE) command, the scale of the plot can be
% controlled.  TYPE takes the following options:
%      'prop': proportion of variance
%      'lod': LOD scale, log base 10; default
%      'llik': log likelihood ratio scale, base e; scale in which theyare
%              calculated 
%      'F': F statistic scale
% PLOTMAINSCAN(ONELOD,TYPE,...)
% One can also specify another option in a name-value pair.  The option
% is 'plotmarkerpos', which is either 'no' or a marker data structure
% from the output of READDATA or IMPORTDATA.
%
% See also: MAINSCAN, PLOTPAIRSCAN.  
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.838 $ $Date: 2002/01/08 23:32:05 $	

  cross = onelod(1).cross;
  
  if( nargin==1 )
    type='lod';
  end

  plotmarkerpos = argchk( varargin(:) );
  
  n = onelod(1).n;

  
  hold on;
  lod = [];
  nchroms = length( onelod );
  cumlen = 0;
  chromlen = zeros( 1, nchroms );
  chroms = zeros( 1, nchroms );
  for( i=1:nchroms )
    % spacing between pseudomarkers
    spacing = 0.1;
    % current chromosome lod
    thislod = onelod(i).lod;
    % overall lod is augmented
    % lod = [ lod thislod ];
    % positions of this 
    thispos = cumlen + onelod(i).mpos + (i-1)*spacing;
    if( isstruct(plotmarkerpos) )
      thismkpos = cumlen + plotmarkerpos(i).mpos + (i-1)*spacing;
      plot( 100*thismkpos, 0, 'k+' )
    end
    cumlen = cumlen + onelod(i).chromlen;
    chromlen(i) = onelod(i).chromlen;
    chroms(i) = onelod(i).chrid;

    
    
    switch type
     case 'lod'
      % plot( thispos, thislod/log(10), 'linewidth', 1 );
      plot( 100*thispos, thislod );      
     case 'prop'
      plot( 100*thispos, 1-exp(-2*thislod*log(10)/n) );     
     case 'llik'
      plot( 100*thispos, thislod*log(10) );     
     case 'F'
      if( cross == 'bc' )
	df = 1;
      elseif( cross == 'f2' )
	df = 2;
      end
      plot( 100*thispos, lod2f(thislod,n,df) );     
    end

  end

  set(gca,'Ticklength',[0 0])
  labelpos = 100*cumsum(chromlen+spacing);
  set(gca,'XTick',labelpos-100*chromlen/2);
  set(gca,'XTickLabel',chroms);  
  %xlim( [ -spacing labelpos(end) ] );
  xlim( [ 0 labelpos(end) ] );  
  hold off
 
  switch type
   case 'lod'
    ylabel( 'LOD score' );
   case 'prop'    
    ylabel( 'Proportion of variance' );
   case 'llik'
    ylabel( 'Log likelihood' );
   case 'F'
    ylabel( 'F' );
  end
  
  xlabel( 'Chromosome number' );


function plotmarkerpos = argchk( args )

  % default values
  plotmarkerpos = 'no';
  
  nvarargs = length( args );
  if( nvarargs>0 )
    nstep = 1;
    while( nstep <= nvarargs )
      argtype = args{nstep};
      switch argtype
       case 'plotmarkerpos'
	plotmarkerpos = args{nstep+1};
	nstep = nstep+2;
       otherwise 
	error('Argument unknown');
      end
    end
  end
  
