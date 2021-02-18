function plotmainlocalize( onelod, varargin )
% PLOTMAINLOCALIZE Plots posterior distribution of single QTL.  
%
% PLOTMAINLOCALISE(LOD)
% PLOTMAINLOCALISE(LOD,...)
% LOD = LOD structure returned from MAINSCAN; should be one element of
%       the structure of arrays returned  
% One can also specify another option in a name-value pair.  The option
% is 'plotmarkerpos', which is either 'no' or a marker data structure
% from the output of READDATA or IMPORTDATA.
%  
% See also: MAINSCAN.
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2002/01/07 21:24:19 $	
  
  plotmarkerpos = argchk (varargin(:) );
  
  lod = onelod.lod * log(10);
  maximum = max( lod );
  const = trapz( onelod.mpos, exp(lod-maximum) );
  hold on
  plot( 100*onelod.mpos, exp(lod-maximum)/const );
  if( isstruct(plotmarkerpos) )
    thismkpos = plotmarkerpos.mpos;
    plot( 100*thismkpos, 0, 'k+' )
  end
    
  
%  set(gca,'XTick', unique(qtldata(onelod.chrid).mpos) );
%  set(gca,'XTickLabel',round(unique(qtldata(onelod.chrid).mpos)*100)/100);  
%  set(gca,'XTick', onelod.mpos );
%  set(gca,'XTickLabel',onelod.mpos);  
  
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
  
