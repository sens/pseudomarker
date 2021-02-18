function imputed = impute( markerdata, spacing, npages, chroms, cross )
% IMPUTE Function to compute an imputed data set from the observed markers.
%
% FAKE = IMPUTE(MARKERDATA)    
% FAKE = IMPUTE(MARKERDATA,SPACING)    
% FAKE = IMPUTE(MARKERDATA,SPACING,NPAGES)  
% FAKE = IMPUTE(MARKERDATA,SPACING,NPAGES,CHROMS)  
% FAKE = IMPUTE(MARKERDATA,SPACING,NPAGES,CHROMS,CROSS)
%  
% MARKERDATA = observed marker data as returned from the READDATA
%              function
% SPACING = spacing between pseudomarkers in Morgans; default 0.1 Morgans
% NPAGES = number of pseudomarker datasets to be imputed; default 16 for
%          spacing of 10cM, 64 for 5cM, and 256 for 2cM
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
% See also READDATA, IMPUTE2.  

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2002/01/04 22:41:04 $	

nchrom = length( markerdata );

if( nargin <= 1 )
  spacing = 0.1;
end

if( nargin <= 2 )
  mult = floor( 0.1/spacing );
  if( mult >=4 )
    mult = 4;
  end
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
  mpos = 0:spacing:(markerdata(chroms(i)).chromlen);
  igeno = jointfiller ( markerdata(chroms(i)).geno, ...
			markerdata(chroms(i)).mpos, ...
			mpos, npages, cross );
  imputed(i).chrid = markerdata(chroms(i)).chrid;
  imputed(i).igeno = igeno;
  imputed(i).mpos = mpos;
  
  % make marker names
  mnames = cell( length(mpos), 1 );
  for( j=1:length(mpos) )
    mnamesi{j} = strcat( 'c', num2str(i), 'mk', num2str(j) );
  end
  imputed(i).mnames = mnamesi;

  imputed(i).chromlen = markerdata(chroms(i)).chromlen;
  clear mpos, igeno;

end
  fprintf( '\n' );