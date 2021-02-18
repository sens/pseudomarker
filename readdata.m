function [ y, qtldata, pnames ] = readdata( varargin )
% READDATA Function to load data from a data directory
%
% [ Y, M, PNAMES ] = READDATA   
% [ Y, M, PNAMES ] = READDATA( DIR )   
% [ Y, M, PNAMES ] = READDATA( DIR, GMISS, PMISS ) 
% [ Y, M, PNAMES ] = READDATA( PHENO, GENO, CHRID, DIR, GMISS, PMISS )   
% DIR = directory name; default '.'
% GMISS = code for missing marker data; default 9
% PMISS = code for missing phenoytype data; default -999
%
% Reads data from the DIR directory.  It should contain the following
% files: 
% 1/ geno.dat: genotype data; individuals as rows, markers as columns
% 2/ pheno.dat: phenotype data; individuals as rows
% 3/ mnames.txt: names of the markers in geno.dat; rows are markers
% 4/ chrid.dat: chromosome id of all the markers; rows are markers
% 5/ markerpos.txt: file with the marker names and their positions in
%    CENTIMORGANS.  This file may contain marker names that are not
%    typed.
%
% The program will read all the files, order the markers in increasing map
% position and place them in the array of structures, M.  Each element of
% the array contains the marker data from a particular chromosome.  The
% program also sets the chromosome length by rounding up to the nearest
% 10cM.
%
% The resulting marker data structure array has the following fields:
% CHRID = chromosome number
% GENO = matrix of genotypes; can be missing
% MNAMES = vector marker names
% MPOS = vector of marker positions
% CHROMLEN = chromosome length which is estimated by rounding up to the
%            nearest 10cM from the rightmost marker

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/09/20 19:28:16 $	

if( nargin == 0 )
  dir = '.';
  gmiss = 9;
  pmiss = -999;
end

if( nargin == 1 )
  dir = varargin{1};
  gmiss = 9;
  pmiss = -999;
end

if( nargin == 3 )
  dir = varargin{1};
  gmiss = varargin{2};
  pmiss = varargin{3};
end

if (nargin==6)
  dir = varargin{4};
  gmiss = varargin{5};
  pmiss = varargin{6};
end

if( nargin <= 3 )  
  geno = load( strcat( dir, '/geno.dat' ) );
else 
  geno = varargin{2};
end

geno(find(geno == gmiss)) = NaN;
[ n, nmarkers ] = size(geno);

% read the phenotype data
if( nargin<=3 )
  y = load( strcat( dir, '/pheno.dat' ) );
else
  y = varargin{1};
end

y( find(y==pmiss) ) = NaN;
%id = load( strcat( dir, '/id.dat' ) );

% chromosome id
if( nargin<=3 )
  chrid = load( strcat( dir, '/chrid.dat' ) );
else
  chrid = varargin{3};
end

%m = length(chrid);

mnames = textread( strcat( dir, '/mnames.txt' ), ...
		    '%s', 'commentstyle','matlab' );
% read the marker names and marker positions
[ allmnames mpos ] = textread( strcat( dir, '/markerpos.txt'), ...
			      '%s %f', 'commentstyle','matlab' );

% adjust data to remove those components for whom the marker positions
% are not known
%geno = geno(:, find(~isnan(mpos)) );
%mnames = mnames( find(~isnan(mpos)) );
%chrid = chrid( find(~isnan(mpos)) );
%mpos( find(isnan(mpos)) ) = [];

newmnames = mnames;
newmpos = mpos;
mpos = [];
deletelist = [];

% i goes through the list of existing markers
for( i=1:nmarkers )
  % which marker name matches the ith marker genotyped
  www = strcmp( allmnames, mnames(i) );
  % if no marker in database matches a typed marker
  if( sum( www ) == 0 )
    % augment the delete list
    fprintf( strcat('Warning: Marker', ' ', mnames{i},...
		    ' position not found.\n' ) );
    deletelist = [ deletelist i ];
  else
    % marker position of the found marker
    mpos(i-length(deletelist)) = newmpos( www );
  end
end

% delete columns which did not have an entry in the marker positions
% database 
geno( :, deletelist ) = [];
chrid( deletelist ) = [];
newmnames( deletelist ) = [];
mnames = newmnames;

%geno = geno(:, find(~isnan(mpos)) );
%mnames = mnames( find(~isnan(mpos)) );
%chrid = chrid( find(~isnan(mpos)) );
%mpos( find(isnan(mpos)) ) = [];


%mpos(chrid==15 )
%mnames(chrid==15)
% ----
chroms = unique( chrid );
nchrom = length( chroms );

qtldata = repmat( struct( 'chrid', 0, 'geno', [], 'mnames', [], ...
			  'mpos', [], 'chromlen', 0 ), 1, nchrom );

for( i=1:nchrom )
  qtldata(i).chrid = chroms(i);
  qtldata(i).geno = geno( :, find(chrid == chroms(i)) );
  qtldata(i).mnames = mnames( find(chrid == chroms(i)) );
  qtldata(i).mpos = mpos( find(chrid == chroms(i)) )/100; % convert to Morgans
  tmp = qtldata(i).mpos;
  % makes up chromosome length
  qtldata(i).chromlen = ceil( max( tmp ) / 0.1 ) * 0.1;
  [ qtldata(i).mpos, ord ] = sort( tmp );
  tmp = qtldata(i).geno( :, ord );
  qtldata(i).geno = tmp;
  tmp = qtldata(i).mnames( ord );
  qtldata(i).mnames = tmp;
  
  clear tmp, ord;
end

fid = fopen( strcat( dir, '/pnames.txt' ) );

if( fid == -1 )
  ntr = size(y,2);
  pnames = cell(ntr,1);
  for( i=1:ntr )
    pnames(i) = {strcat( 'pheno', num2str(i) )};
  end
else
  fclose( fid );
  pnames = textread( strcat( dir, '/pnames.txt' ), '%s' );
end

