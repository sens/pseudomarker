function plotgeno( m, varargin )
% PLOTGENO Plot genotype pattern and missing data pattern
%       
% PLOTGENO(M)
% PLOTGENO(M,...)  
% M = marker data structure returned by READDATA
% Optional arguments can be passed as name-value pairs.
% sortby: vector phenotypes to sort the individuals by
% cap: caption string such as 'Genotype pattern'
% genolabel: cell array of genotype labels such as { 'B6', 'Het', Missing' }
%
% For backcross the two genotypes are red and turquoise; yellow ochre is
% missing.  For intercross the colors are red, dark brown, and turquoise for
% homozygous AA, heterozygous AB, and homozygous BB, respectively; yellow
% ochre is missing.

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/09/24 19:58:34 $	
  
  nchroms = length( m );
  cumlen = 0;
  nmarks = zeros( 1, nchroms );
  
  if( isfield( m(1), 'geno' ) )
    n = size(m(1).geno,1);
    imputed = 0;
  elseif ( isfield( m(1), 'igeno' ) )
    n = size(m(1).igeno,1);
    imputed = 1;
  else
    error( 'This is not a marker data structure.' );
  end

  cap = 'Genotype pattern';
  y = [];

  nvarargin = length(varargin);
  if( length(varargin)>0 )
    nstep = 1;
    while( nstep <= nvarargin )
      argtype = varargin{nstep};
      switch argtype
       case 'sortby'
	y = varargin{nstep+1};
	nstep = nstep+2;
       case 'cap'
	cap = varargin{nstep+1};
	nstep = nstep+2;
       case 'genolabel'
	genolabel = varargin{nstep+1};
	nstep = nstep+2;
       otherwise
	error( 'Cannot recognize option.' )
	nstep = nargin;
      end
    end
  end
  
  
  
  if( ~isempty(y) )
    [ysrt,ord]=sort(y);
    m = subsetgeno1(m,ord);
  end
  
  % find number of markers in each chromosome
  for( i=1:nchroms )
    nmarks(i) = length(m(i).mpos);
  end
  nmarks = [ 0 nmarks ];
  
  % total number of markers
  nnmarks = sum(nmarks);
  nmarks = cumsum( nmarks );
  midchrom = sum( nmarks < nmarks(end)/2 );
  
  totpos = (nnmarks + nchroms - 1)*1.33;
  
  cross = guesscross(m);
  
  %  cmap_bc = [ 1 0 0 ; 0 0 1; 0.8 0.8 0.8 ];
  cmap_bc = [ 0.85 0 0 ; 0 0.85 0.85; 0.9 0.7 0.1 ];
  %  cmap_f2 = [ 1 0 0 ; 0.9 1 0.3; 0 0 1; 0.8 0.8 0.8];
  cmap_f2 = [ 0.9 0 0 ; 0.45 0.3 0.3 ; 0 0.85 0.85; 0.9 0.7 0.1];
  
  if( cross == 'bc' )
    cmap = cmap_bc;
  elseif( cross == 'f2' )
    cmap = cmap_f2;
  end
  
  
  for( i=1:nchroms )
    subplot( 'position', [ 0.1 + (nmarks(i)+i-1)/totpos 0.15 ...
		    ( nmarks(i+1)-nmarks(i) ) / totpos 0.75 ] )
    
    if( imputed == 0 )
      ggg = m(i).geno;
    else
      ggg = m(i).igeno(:,:,1);
    end
    
    if( cross == 'bc' )
      ggg( isnan(ggg) ) = 2;
      if( exist('genolabel') == 0 )
	genolabel = { 'A', 'B', 'Missing' };
      end
    elseif( cross == 'f2' )
      ggg( isnan(ggg) ) = 3;
      if( exist('genolabel') == 0 )
	genolabel = { 'AA', 'AB', 'BB', 'Missing' };
      end
    end
    
%    image( ggg, 'cdatamapping', 'scaled' );
    image( ggg+1 );    
    if( i==1 )
      set(gca,  'xticklabel','', 'ticklength', [ 0 0 ] );
      xlabel( num2str(m(i).chrid) );
      ylabel( 'Individual number' );
      colormap( cmap );
    else

      set(gca,  'yticklabel','', 'xticklabel','', 'ticklength', [ 0 0 ] );
      xlabel( num2str(m(i).chrid) );
      colormap( cmap );
    end
    if( i==midchrom )
      text( (nmarks(i+1)-nmarks(i))/2,n*1.1, 'Chromosome number', ...
	    'horizontalalignment', 'center' );

      
      text( (nmarks(i+1)-nmarks(i))/2,-n*0.05, cap, ...
	    'horizontalalignment', 'center' );
    end
    
  end

    
  subplot( 'position', [ 0.9  0.15 ...
		    0.02 0.75 ] )
  
  if( cross == 'f2' )
    image( [ 1 2 3 4 ]' );
    set( gca, 'ytick', [ 1 2 3 4 ] );    
    set( gca, 'yticklabel', genolabel);    
    %set( gca, 'yticklabel', {'AA', 'AB', 'BB', 'Missing'});    
  else
    image( [ 1 2 3 ]' );
    set( gca, 'ytick', [ 1 2 3 ] );    
    %set( gca, 'yticklabel', {'A', 'B', 'Missing'});    
    set( gca, 'yticklabel', genolabel );    	
  end
    
  set( gca, 'yaxislocation', 'right' );
  set( gca, 'xtick', 1 );
  set( gca, 'xticklabel', 'Genotype' );
  set( gca, 'tickdir', 'out' )    
  set( gca, 'ticklength', [ 0 0 ] )
    

