function imputed = impute2( markerdata, spacing, npages, chroms, cross )
% IMPUTE2 Impute genotypes at pseudomarkers and typed locations.
%
% FAKE = IMPUTE(MARKERDATA)    
% FAKE = IMPUTE(MARKERDATA,SPACING)    
% FAKE = IMPUTE(MARKERDATA,SPACING,NPAGES)  
% FAKE = IMPUTE(MARKERDATA,SPACING,NPAGES,CHROMS)  
% FAKE = IMPUTE(MARKERDATA,SPACING,NPAGES,CHROMS,CROSS)
%  
% MARKERDATA = observed marker data as returned from the READDATA
%              function
% SPACING = if spacing between markers is greater than SPACING, create
%           pseudomarkers in between; default is 10 Morgans which
%           essectially means impute only at the typed markers
% NPAGES = number of pseudomarker datasets to be imputed; default 1 for
%          spacing of 100cM, 16 for 10cM, 64 for 5cM, and 256 for 2cM
% CHROMS = which chromosomes to impute; if omitted all chromosomes in
%          MARKERDATA are imputed
% CROSS = cross type; use 'bc' for backcross, 'f2' for intercross; if
%         absent the maximum value of the genotypes on the first
%         chromosome are used to decide if 'bc' or 'f2' are appropriate
%
% The output is an array of structures with the following fields:
% CROSS = type of cross
% CHRID = chromosome number
% IGENO = 3-dimensional array of imputed genotypes; the first index
%         (rows) is individuals, the second index (cols) is
%         pseudomarkers, and the third index (pages) is imputation number
% MPOS = pseudomarker positions
% CHROMLEN = chromosome length
%  
% See also READDATA, IMPUTE.  

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2002/01/09 22:49:58 $	

nchrom = length( markerdata );

if( nargin <= 1 )
  spacing = 10;
end

if( nargin <= 2 )
  if( spacing == 10 )
    % multiplier
    mult = 1/4;
  else
    % multiplier
    mult = floor( 0.1/spacing );
  end
  if( mult >=4 )
    % prevent too many imputations by making a default upper bound for
    % the multiplier; the user can override this, of course
    mult = 4;
  end
  % how many imputations to make by default
  npages = 16*mult^2;
end

if( nargin<=3 )
  chroms = 1:nchrom;
else
  nchrom = length(chroms);
end

if( nargin <= 4 )
  cross = guesscross( markerdata );
end



imputed = repmat( ...
    struct( 'cross', cross, 'chrid', 0, 'igeno', [], ...
	    'mpos', [], 'chromlen', 0 ), 1, nchrom );


for( i=1:length(chroms) )
  fprintf( '%d ', chroms(i) );
  mnames = markerdata(chroms(i)).mnames;
  [mpos,mnames] = makepos( markerdata(chroms(i)).mpos, spacing, mnames, ...
			   chroms(i) );
  igeno = jointfiller ( markerdata(chroms(i)).geno, ...
			markerdata(chroms(i)).mpos, ...
			mpos, npages, cross );
  imputed(i).chrid = markerdata(chroms(i)).chrid;
  imputed(i).igeno = igeno;
  imputed(i).mpos = mpos;
  
  imputed(i).mnames = mnames;

 
  imputed(i).chromlen = markerdata(chroms(i)).chromlen;
  clear mpos, igeno;
end
  fprintf( '\n' );


function [mp,mn] = makepos( mpos, spacing, mnames, chrid )
% MAKEPOS Make positions and names for imputations

  % mn = cell( length(mpos), 1 )
  % mp = zeros( length(mpos), 1 )  
  mp = mpos(1);
  mn = mnames(1);
  fakemk = 0;
  for( i=2:length(mpos) )
    d = mpos(i)-mpos(i-1);
    n = floor( d/spacing );
  
    if(n>0)
      delta = d/(n+1);
      mmm = mpos(i-1) + (1:n)*delta;
      mp = [ mp mmm mpos(i) ];
      tmpmnames = cell( 1, n );
      for( j=1:n )
	tmpmnames{j} = strcat( 'c', num2str(chrid), 'mk', num2str(j+fakemk) );
      end
      
      fakemk = fakemk + n;
      mn = [ mn tmpmnames mnames(i) ];
    else
      mp = [ mp mpos(i) ];
      mn = [ mn mnames(i) ];
    end
  end
  
  mp = sort(mp);
  mn = mn';
% $$$   mp = mpos;
% $$$   for( i=2:length(mpos) )
% $$$     d = mpos(i)-mpos(i-1);
% $$$     n = floor( d/spacing );
% $$$ 
% $$$     if(n>0)
% $$$       delta = d/(n+1);
% $$$       mmm = mpos(i-1) + (1:n)*delta;
% $$$       mp = [ mp mmm ];
% $$$     end
% $$$   end
% $$$ 
% $$$   mp = sort(mp);
  