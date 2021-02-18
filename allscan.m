function [lod1,lod2] = allscan(y,fake,pnames,varargin)
% ALLSCAN Makes scans of all phenotypes in a data

  % default values of arguments
  plot = 'no'; % other values are 'eps', 'jpeg', 'both'.
  pair = 'no';
  report = 'no';
  
  nargin = length(varargin);
  if( length(varargin)>0 )
    nstep = 1;
    while( nstep <= nargin )
      argtype = varargin{nstep};
      switch argtype
       case 'dir'
	dir = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'plot'
	plot = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'pair'
	pair = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'em'
	em = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'flagbf'
	flagbf = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'flaglod'
	flaglod = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'report'
	report = varargin{ nstep+1 };
	nstep = nstep+2;
       otherwise
	error( 'Cannot recognize option.' )
	nstep = nargin;
      end
    end
  end
  
% number of phenotypes  
npheno = size(y,2);

% mainscan lod structure for the first phenotype
lod1 = mainscan( y(:,1), [], [], fake );
% make an array of structures for one-dimensional LOD scores
lod1 = repmat( lod1, [ npheno 1 ] );

% if pairscans are desired
if( strcmp( pair, 'yes' ) )
  % make pairscan for first phenotype
  lod2 = pairscan( y(:,1), [], [], fake );
  % make an array of structures for two-dimensional LOD scores  
  lod2 = repmat( lod2, [ 1 1 npheno ] );
end

for( i=1:npheno )

  if( i>1 )
    lod1(i,:) = mainscan( y(:,i), [], [], fake );
    if( strcmp( pair, 'yes' ) )
      lod2(:,:,i) = pairscan( y(:,i), [], [], fake );
    end
  end
  

  % if reports are desired
  if( ~strcmp(report,'no'))
    % if only mainscan report is desired
    if( length(report) == 1 )
      reportscan(lod1(i,:),report);
    elseif( length(report) == 4 )
      % if pairscan reports are desired report from both mainscan and
      % pairscan 
      reportscan(lod1(i,:),report(1));      
      reportscan(lod1(i,:),lod2(:,:,i),report(2:end));
    end
  end
  
  % if plots are desired
  if( ~strcmp(plot,'no') )
    % mainscan
    figure('Position',[1 1 600 250],'paperpositionmode','auto')
    plotmainscan(lod1(i,:));
    ylabel( 'LOD' );
    xlabel( 'Chromosome number' );
    title( pnames(i) );  

    % if plots are desired
    if( strcmp(plot,'jpeg') | strcmp(plot,'both') )
      print( '-djpeg',  strcat(dir,'/',pnames{i},'1') );
    end
    if( strcmp(plot,'eps') | strcmp(plot,'both') )
      print( '-depsc',  strcat(dir,'/',pnames{i},'1') );
    end
    close
    
    % if pairscans are desired
    if( strcmp(pair,'yes') )
      figure;
      plotpairscan(lod2(:,:,i));
      ylabel( 'Chromosome number' );
      xlabel( 'Chromosome number' );
      colorbar
      title( pnames(i) );  
      
      % if plots are desired
      if( strcmp(plot,'jpeg') | strcmp(plot,'both') )
	print( '-djpeg',  strcat(dir,'/',pnames{i},'2') );
      end
      if( strcmp(plot,'eps') | strcmp(plot,'both') )
	print( '-depsc',  strcat(dir,'/',pnames{i},'2') );
      end
      close
      
    end
  end
  
end


