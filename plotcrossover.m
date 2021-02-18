function recombfrac = plotcrossover( m, varargin )
% PLOTCROSSOVER Plot and summarize crossover information on a chromosome
%       
% THETA = PLOTCROSSOVER(M)
% THETA = PLOTCROSSOVER(M,...)    
% THETA is a vector with two rows.  The first row is the estimated
% recombination fraction and the second row is the recombination fraction
% as specified by the supplied marker positions.
%  
%  M = marker data structure returned by READDATA
% The other options can be one or more of the following in the 
% following format: 
% Caption: 'cap', CAP
% Animal Ids: 'id', ID
% Subset of animals: 'subset', SUBSET
% If missing are to be filled in: 'filledin'
% Sort by phenotype: 'pheno', PHENO
% Do not plot animal Ids: 'noid'
% CAP is a string for the caption of the plot
% ID is a numeric vector or a vector of strings
% SUBSET is a numeric vector indicating the subset of animals to be
% displayed   
% PHENO is a numeric vector of phenotypes along which the animals have to
% be ordered
%  

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/09/28 00:17:38 $	

  if( isfield( m(1), 'geno' ) )
    n = size(m(1).geno,1);
    imputed = 0;
  elseif ( isfield( m(1), 'igeno' ) )
    n = size(m(1).igeno,1);
    imputed = 1;
  else
    error( 'This is not a marker data structure.' );
  end

   
  % default values
    noid = 1;
    cap = 'Crossover pattern';
    id = num2str( 1:n );
    subset = 1:n;
    y = 1:n;
    filledin = 0;
    
  % extra arguments to customize the plot; we have to let them know if
  % they want the ids of the individuals on the side, if they want to
  % plot a subset of the individuals
 
  nargin = length(varargin);
  if( length(varargin)>0 )
    nstep = 1;
    while( nstep <= nargin )
      argtype = varargin{nstep};
      switch argtype
       case 'subset'
	subset = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'pheno'
	y = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'id'
	id = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'cap'
	cap = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'filledin'
	filledin = 1;
	nstep = nstep+1;
       case 'noid'
	id = repmat(' ', 1,n );
	nstep = nstep+1;
       otherwise
	error( 'Cannot recognize option.' )
	nstep = nargin;
      end
    end
  end
  
  
  
  nchroms = length( m );
  cumlen = 0;
  nmarks = zeros( 1, nchroms );
  
 
  if( nargin>1 )
    y = y(subset);
    m = subsetgeno(m,subset);
    [ysrt,ord]=sort(y);
    m = subsetgeno(m,ord);
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
  cmap_bc = [ 0 0.85 0.85; 0.85 0 0; 1 1 1 ];
  %  cmap_f2 = [ 1 0 0 ; 0.9 1 0.3; 0 0 1; 0.8 0.8 0.8];
  cmap_f2 = [ 0 0.85 0.85; 0.45 0.3 0.3 ; 0.9 0 0 ; 1 1 1];
  
  if( cross == 'bc' )
    cmap = cmap_bc;
  elseif( cross == 'f2' )
    cmap = cmap_f2;
  end
  
  
  for( i=1:nchroms )
    subplot( 'position', [ 0.1 + (nmarks(i)+i-1)/totpos 0.15 ...
		    ( nmarks(i+1)-nmarks(i) ) / totpos 0.65 ] )
    
    if( imputed == 0 )
      ggg = abs(diff(m(i).geno,1,2));
    else
      ggg = abs(diff(m(i).igeno(:,:,1),1,2));
    end
    
    if( cross == 'bc' )
      ggg( isnan(ggg) ) = 2;
    elseif( cross == 'f2' )
      ggg( isnan(ggg) ) = 3;
    end
    
    mkspacing = diff(m(i).mpos); % marker spacings
    
    if( cross == 'bc' )
      ncrossover = sum(ggg) - 2*sum(ggg==2); % number of crossovers
      nnotmiss = sum( ggg~=2 ); % number of not missing cases
    elseif( cross == 'f2' )
      ncrossover = sum(ggg) - 3*sum(ggg==3); % number of crossovers
      nnotmiss = sum( ggg~=3 ); % number of not missing cases
    end

    % if data were "filled in" using Bev's strategy then estimate the
    % recombination fraction by dividing by the total number of animals;
    % else normalize by the number of not missing genotypes
    if( filledin == 0 )
      if( cross == 'bc' )
	theta = ncrossover./nnotmiss; % estimate recombination fraction
      elseif( cross == 'f2' )
	theta = ncrossover./(2*nnotmiss);
      else 
	error( 'Unknown cross type' )
      end

    else
      if( cross == 'bc' )
	theta = ncrossover/n;
      elseif( cross == 'f2' )
	theta = ncrossover/(2*n);
      else 
	error( 'Unknown cross type' )
      end
    end
    
    fprintf( '\tChromosome number %d:\n', m(i).chrid );
    recombfrac = [ theta; recomb(mkspacing) ];
    
    %    image( ggg, 'cdatamapping', 'scaled' );
    % plot the crossover pattern
    image( ggg+1 );
    
    % make the marker labels
    for( k=1:(size(ggg,2)+1) )
      h = text( k-0.5, 0, m(i).mnames(k) );
      set(h,'rotation',90)
    end
    
    
    if( i==1 )
      set( gca, 'yticklabel', '' );
      set(gca,  'xticklabel','', 'ticklength', [ 0 0 ] );
      xlabel( num2str(m(i).chrid) );
      %ylabel( 'Individual number' );
      colormap( cmap );
    
          % make the animal ids
    for( k=1:size(ggg,1) )
      if( isnumeric(id(k)) )
	idlabel = num2str(id(k));
      else
	idlabel = id(k);
      end
      
      h = text( 0, k, idlabel );
      set( h, 'fontsize', 10, 'horizontalalignment', 'right' );
    end

    else

      set(gca,  'yticklabel','', 'xticklabel','', 'ticklength', [ 0 0 ] );
      xlabel( num2str(m(i).chrid) );
      colormap( cmap );
    end
    if( i==midchrom )
      text( (nmarks(i+1)-nmarks(i))/2,n*1.1, 'Chromosome number', ...
	    'horizontalalignment', 'center' );

      
      text( (nmarks(i+1)-nmarks(i))/2,-n*0.25, cap, ...
	    'horizontalalignment', 'center' );
    end
    
  end

    
  subplot( 'position', [ 0.9  0.15 ...
		    0.02 0.65 ] )
  
  if( cross == 'f2' )
    image( [ 1 2 3 4 ]' );
    set( gca, 'ytick', [ 1 2 3 4 ] );    
    set( gca, 'yticklabel', {'0', '1', '2', 'Missing'});    
  else
    image( [ 1 2 3 ]' );
    set( gca, 'ytick', [ 1 2 3 ] );    
    set( gca, 'yticklabel', {'0', '1', 'Missing'});    
  end
    
  set( gca, 'yaxislocation', 'right' );
  set( gca, 'xtick', 1 );
  set( gca, 'xticklabel', 'Genotype' );
  set( gca, 'tickdir', 'out' )    
  set( gca, 'ticklength', [ 0 0 ] )
    

