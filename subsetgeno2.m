function newmk = subsetgeno2( m, chrid, pos )
% SUBSETGENO2 Subsets a marker genotype structure array using marker positions 
%
% M1 = SUBSETGENO1(M0,CHRID,MPOS)
% M0 = mother genotype structure array returned from READDATA or imputed
%      structure array returned by IMPUTE
% CHRID = chromosome numbers to be subsetted
% MPOS = marker positions in Morgans; if a one-dimensional, then closest
%        marker is selected else all markers in the range are selected  
%
% This function is useful for studying genotype patterns or for
% subsetting individuals.
% 
% Example: subsetgeno2( m, [ 1 5 8 ] );
%          subsetgeno2( m, [ 1 5 8 ], [ 0.41 0.23 0.34 ] );
%          subsetgeno2( m, [ 1 5 8 ], [ 0.41 0.78; 0.32 0.45 0.63 0.92 ] );  
%
% See also: SUBSETGENO, SUBSETGENO1, SUBSETGENO3..
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2002/01/04 22:48:22 $	

  % is this a marker or pseudomarker structure
  if( isfield( m(1), 'geno' ) )
    imputed = 0;
  elseif ( isfield( m(1), 'igeno' ) )
    imputed = 1;
  else
    error( 'This is not a marker data structure.' );
  end
  
  % number of chromosomes and number of individuals
  nchroms = length(m);
  if( imputed == 0 )
    n = size(m(1).geno,1);
  else
    n = size(m(1).igeno,1);
  end
  
  % subset the chromosomes according to chromosome id
  newmk = [];
  % go through the chromosomes selected by the user
  for( i=1:length(chrid) )
    % indicator whether a chromosome was found or not
    found = 0;
    % initialize couter of marker data structure
    j = 1;
    while( found == 0 )
      % if counter exceeds the number of chromosomes, not found
      if( j > nchroms )
	found = -1;
      end
      % if chromosome found, assign and end search
      if( chrid(i) == m(j).chrid )
	newmk = [ newmk m(j) ];
	found = 1;
      end
      % increase counter
      j = j+1;
    end
    % if no matching chromosome found warn the reader
    if( found == -1 )
      fprintf( 'Chromosome %d not found.\n', chrid(i) );
    end
  end
  
  % if no matching chromosomes found, return error
  if( isempty(newmk) )
    error( 'No matching chromosomes.' )
  end
  

  if( nargin == 3 )
    m = newmk;
    % number of chromosomes selected
    nchroms = length( newmk );
    % numner of rows and columns in pos
    [ nr, nc ] = size( pos );

    % if pos is one-dimensional
    if( length( pos ) == nr*nc )
      if( nr*nc ~= nchroms )
	error( 'Chromosomes selected do not match number of positions.' ...
	       );
      end
      % go through the chromosomes
      for( i = 1:nchroms )
	% get the index of the nearest marker position
	[ mdist, idx ] = min( abs( m(i).mpos - pos(i) ) );
	if( imputed==0 )
	  newmk(i).geno = m(i).geno(:,idx);
	  newmk(i).mnames = m(i).mnames(idx);
	else
	  newmk(i).igeno = m(i).igeno(:,idx,:);
	  newmk(i).mnames = m(i).mnames(idx);
	end
	newmk(i).mpos = m(i).mpos(idx);
      end
    else
      % if positions ate not one-dimensional
      if( nr ~= nchroms )
	error( 'Chromosomes selected do not match number of positions.' ...
	       );
      end
      for( i = 1:nchroms )
	idx = find( (m(i).mpos >= pos(i,1) ) & ( m(i).mpos <= pos(i,2) ) );
	if( imputed==0 )
	  newmk(i).geno = m(i).geno(:,idx);
	  newmk(i).mnames = m(i).mnames(idx);
	else
	  newmk(i).igeno = m(i).igeno(:,idx,:);
	  newmk(i).mnames = m(i).mnames(idx);
	end
	newmk(i).mpos = m(i).mpos(idx);
      end
    end
  end
  




