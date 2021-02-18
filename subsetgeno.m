function mdata = subsetgeno( m, varargin )
% SUBSETGENO Makes a subset of a marker genotype structure array
%
% M1 = SUBSETGENO(M0,...)
% M0 = mother genotype structure array returned from READDATA or imputed
%      structure array returned by IMPUTE
% Options should be passed in name-value pairs as follows:
% individuals = case number of the individuals
% chrid = chromosome ids of the chromosomes
% mpos = marker positions in CENTIMORGANS
%        if one-dimensional vector, nearest marker is selected
%        if two-dimensional range of markers in within the positions are
%        selected; each row is a chromosome and the number of rows must
%        match the number of chromosomes in chrid
% npages = imputation numbers to be selected
%  
% See also: SUBSETGENO1, SUBSETGENO2, SUBSETGENO3.
%  

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/09/24 22:30:27 $	
  
  % get the arguments
  [idx,chrid,mpos,npages]=argchk( m, varargin );
  % if the individuals argument is non-empty
  if( ~isempty(idx) )
    m = subsetgeno1(m,idx);
  end

  % if the chrid argument is non-empty
  if( ~isempty(chrid) )
    m = subsetgeno2(m,chrid);
  end

  % if the mpos argument is non-empty
  if( ~isempty(mpos) )
    mpos = mpos/100;
    % but if the chrid is empty make the chrid
    if( isempty(chrid) )
      chrid = zeros(1,length(m))
      for(i=1:length(m))
	chrid(i)=m(i).chrid;
      end
    end
    m = subsetgeno2(m,chrid,mpos);
  end
  
  % if the npages argument is non-empty
  if( ~isempty(npages) )
    m = subsetgeno3(m,npages);
  end
  
  mdata = m;
  
  
function [idx,chrid,mpos,npages]=argchk( m, args )

  if( isfield( m(1), 'geno' ) )
    imputed = 0;
  elseif ( isfield( m(1), 'igeno' ) )
    imputed = 1;
  else
    error( 'This is not a marker data structure.' );
  end

  % default values
  npages = [];
  mpos = [];
  chrid = [];
  idx = [];
  
  nvarargs = length( args );
  if( nvarargs>0 )
    nstep = 1;
    while( nstep <= nvarargs )
      argtype = args{nstep};
      switch argtype
       case 'individuals'
	idx = args{nstep+1};
	nstep = nstep+2;
       case 'chrid' 
	chrid = args{nstep+1};
	nstep = nstep+2;
       case 'mpos' 
	mpos = args{nstep+1};
	nstep = nstep+2;
       case 'npages' 
	npages = args{nstep+1};
	nstep = nstep+2;
       otherwise
	error( 'Argument type not recognized.' );
      end
    end
  end
  
  
