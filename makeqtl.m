function qtl = makeqtl( fake, varargin )
% MAKEQTL Make a QTL data structure
%
% QTL = MAKEQTL(FAKE,...)
% FAKE = mother genotype structure array returned from READDATA or imputed
%      structure array returned by IMPUTE or IMPUTE2
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
% The output looks very similar to a marker or pseudomarker data
% structure, except that it has an added field which is "qtlid" in
% addition to "chrid".  The user can specify more than one QTL on the
% same chromosome, but just repeating a chromsome id in CHRID argument.
%
% Eg. qtl = makeqtl(fake,'chrid', [ 1 1 2 ], 'mpos', [ 30 70 25 ] );
%
% However, if the user decides to specify a range of pseudomarkers (or
% markers for a QTL), it is recommended that the ranges of the QTL not
% overlap.  The function does not check for this.  For example,
%
% qtl = makeqtl(fake,'chrid', [ 1 1 2 ], 'mpos', [ 25 40; 35 75; 20 60 ] );  
%  
% is allowed, but could lead to misleading inferences.
%  
% See also: SUBSETGENO, SUBSETGENO1, SUBSETGENO2, SUBSETGENO3, IMPUTE,
%           IMPUTE2. 
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.833 $ $Date: 2002/01/09 21:16:00 $	
  
qtl = subsetgeno( fake, varargin{:} );

nqtl = length( qtl );

for(i=1:nqtl)
  qtl(i).qtlid = i;
end
