function butterfly( data, binwidth, labels )
% BUTTERFLY Make a butterfly plot for groups of data.
%
% BUTTERFLY( DATA, BINWIDTH )  
% BUTTERFLY( DATA, BINWIDTH, LABELS )
%  
% DATA = Cell array of data vectors or two-column matrix; if cell array,
%        each cell contains a vector of observations corresponding to a
%        group; if matrix, first column should be observations and the
%        second, the group id
% BINWIDTH = Bin width
% LABELS = Cell array of labels; if missing, index is used
%
% This way of plotting is suitable for comparing different sets of
% univariate data.  Vary the BINWIDTH to get the best look.
%  
% Example:
%    x1 = normrnd( 0, 1, 50, 1 );
%    x2 = normrnd( 2, 1, 100, 1 );  
%    butterfly( {x1,x2}, 0.2, {'Group1','Group2'} )
%    id = [ zeros(50,1); ones(100,1) ];
%    y = [ x1; x2 ];  
%    butterfly( [ y id ], 0.2, {'Group1','Group2'} )
%  
% See also HIST, HISTC, SUSHI.

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/07/09 16:43:15 $	

  % if the data is a matrix
  if( isnumeric(data) )
    % first col are the observations
    y = data(:,1);
    % second col are the group ids
    id = data(:,2);
    % get factor levels and number of levels
    levels = unique(id);
    nlevels = length(levels);
    % create cell array of dimension equal to number of levels
    data = cell(1,nlevels);
    % make cell array from the observations and the group ids
    for( i=1:nlevels )
      data{i} = y( find(id==levels(i)) );
    end
  end

  % number of groups
  ngroup = length( data );

  % if the labels are missing
    if( nargin < 3 )
      if( exist('levels')==1 )
	labels = num2str(levels);
      else
	labels = num2cell( 1:ngroup );
      end
    end
  
  
  % number, maximum and minimum values for each group
  big = zeros( 1, ngroup );
  small = zeros( 1, ngroup );
  
  % find the max and min for each group
  for( i=1:ngroup )
    big(i) = max( data{i} );
    small(i) = min( data{i} );
  end
  
  % the global max and min
  bbb = max( big(:) );
  sss = min( small(:) );

  % span of the data
  span = bbb - sss;

  % number of bins
  nbins = ceil( span / binwidth );
  span = nbins * binwidth;
  d = ( span - ( bbb - sss ) ) / 2;
  
  smallest = sss - d;
  biggest = bbb + d;
  h = (biggest-smallest)/nbins;
  
  edges = smallest:h:biggest;
  bincenters = (smallest+h/2):h:(biggest-h/2);
  
  counts = zeros( ngroup, nbins );
  size(counts);
  
  for( i=1:ngroup )
    a = histc( data{i}, edges );
    counts(i,:) = a(1:end-1)';
  end
  
  maxc = max(counts(:));
  
  groupdiff = maxc;
  
  y = repmat( bincenters, ngroup, 1 )';
  x = repmat( (1:ngroup)*maxc, nbins, 1 );
  
%  figure
  hold on
  xlim( [ 0 ngroup+1 ]);
  for( i=1:ngroup )
    for( j=1:nbins )
      y = [ bincenters(j), bincenters(j) ];
      x = i +  (counts(i,j)/(2*maxc)) * [ -1 1 ];
      plot( x,y);
    end
  end
  
  set(gca,'xtick',1:ngroup);

  if( ngroup == 1 )
    set(gca,'xticklabel',{labels});  
  else
    set(gca,'xticklabel',labels);  
  end
  