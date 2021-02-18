function plotsummgeno( m, cap )
% PLOTSUMMGENO Plot summary of genotype pattern and missing data pattern
%       
% PLOTSUMMGENO(M)
% PLOTSUMMGENO(M,CAP)  
% M = marker data structure returned by READDATA
% CAP = string which is title of plot 
%
% For backcross the two genotypes are red and turquoise; yellow ochre is
% missing.  For intercross the colors are red, dark brown, and turquoise for
% homozygous AA, heterozygous AB, and homozygous BB, respectively; yellow
% ochre is missing.
%
% See also: PLOTGENO.
  
  nchroms = length( m );
  cumlen = 0;
  nmarks = zeros( 1, nchroms );
  n = size(m(1).geno,1);
  
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
    
    ggg = m(i).geno;

    if( cross == 'bc' )
      ggg( isnan(ggg) ) = 2;
    elseif( cross == 'f2' )
      ggg( isnan(ggg) ) = 3;
    end
    
    %    image( ggg, 'cdatamapping', 'scaled' );
    image( sort(ggg+1) );    
    
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

      if( nargin==1 )
	cap = 'Genotype pattern';
      end
      
      text( (nmarks(i+1)-nmarks(i))/2,-n*0.05, cap, ...
	    'horizontalalignment', 'center' );
    end
    
  end

    
  subplot( 'position', [ 0.9  0.15 ...
		    0.02 0.75 ] )
  
  if( cross == 'f2' )
    image( [ 1 2 3 4 ]' );
    set( gca, 'ytick', [ 1 2 3 4 ] );    
    set( gca, 'yticklabel', {'AA', 'AB', 'BB', 'Missing'});    
  else
    image( [ 1 2 3 ]' );
    set( gca, 'ytick', [ 1 2 3 ] );    
    set( gca, 'yticklabel', {'A', 'B', 'Missing'});    
  end
    
  set( gca, 'yaxislocation', 'right' );
  set( gca, 'xtick', 1 );
  set( gca, 'xticklabel', 'Genotype' );
  set( gca, 'tickdir', 'out' )    
  set( gca, 'ticklength', [ 0 0 ] )
    

