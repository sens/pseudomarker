function [ y, qtldata, pnames, cnames ] = importdata( filename, varargin )
% IMPORTDATA Function to load data from a text file in spreadsheet format
%
% [Y,M]=IMPORTDATA(FILE)
% [Y,M,PNAMES]=IMPORTDATA(FILE)
% [Y,M,PNAMES,CNAMES]=IMPORTDATA(FILE)  
% [Y,M]=IMPORTDATA(FILE,...)
% [Y,M,PNAMES]=IMPORTDATA(FILE,...)
% [Y,M,PNAMES,CNAMES]=IMPORTDATA(FILE,...)  
%
% Output arguments:
% Y = phenotype matrix
% M = marker data structure (see more about it below)
% PNAMES = phenotype names
% CNAMES = case names
%  
% Additional arguments that can be passed in NAME-VALUE pairs are:
% 1/ delimiter: ('tab') 'space' ',' 'comma' ';' 'semi' '|' 'bar'
% 2/ gcodes: If the genotypes are characters, then this cell array can
% specify what characters are supposed to be coded 0,1,missing (in case
% of backcross) and 0,1,2,missing (in the case if F2 intercross).
% Eg: {'a','b','-'}.
% 3/ gmiss: If genotypes are specified as numbers, this specifies the
% missing genotype code.  The default is 9.
% 4/ pmiss: Missing phenotype code.  The default is -999.  
%  
% This program expects data in a very rigid format.  It should be in
% spreadsheet format in a text file.  An example (tab delimited) file 
% would look like this:
% casenames	pheno1	pheno2	chr1mk1	chr1mk2	chr2mk1	chr2mk2
%                               1       1       2       2
%                               10      23      17      32  
% ind1		54	45	a	a	b	b
% ind2		59	25	a	b	a	a  
% ind3		74	55	b	b	b	a
%
% Row 1: Contains phenotype names and marker names.
% Row 2: Contains BLANKS (empty fields, IMPORTANT) for phenotypes
% and chromosome numbers for the markers.
% Row 3: Contains BLAKS for phenotypes and marker positions in
% CENTIMORGANS for the markers.
% Row 4+: Contains casenames (first column, OPTIONAL), phenotype values
% and genotypes (which can be characters or numbers);
%
% If the casenames are desired, the first column can be text and will give
% the case names.
% 
% The program will read the file, order the markers in increasing map
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
%
% SEE ALSO: READDATA.
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.835 $ $Date: 2002/01/12 02:18:45 $	

  % should casenames be returned
  if( nargout == 4 )
    casenames = 'yes';
  else 
    casenames = 'no';
  end
  
  % check the input arguments
  [delimiter,gcodes,gmiss,pmiss] = argchk( varargin );
  dlm = delimchk( delimiter );
  
  fid = fopen( filename, 'r' );
  s = fgetl( fid );

  % first line
  npheno = 0;
  nmarkers = 0;
  nfields = length( findstr(s,dlm) ) + 1;
  firstline = cell( nfields, 1 );
  for( i=1:nfields )
    [token,rem] = strtok(s,dlm);
    firstline{i} = token;
    s = rem;
  end


  % chromosome ids
  chrid = zeros( nfields, 1 );
  s = fgetl( fid );
  i = 1;
  rem = 'Something';
  while( ~isempty(rem) )
    [token,rem] = strtok(s,dlm);
    s = rem;
    num = str2num(token);
    if( ~isnan(num) )
      chrid(i) = num;
      i = i+1;
    end
  end

  chrid(i:end) = [];

  
  % marker positions
  mpos = zeros( 1, nfields );
  s = fgetl( fid );
  i = 1;
  rem = 'Something';
  while( ~isempty(rem) )
    [token,rem] = strtok(s,dlm);
    s = rem;
    num = str2num(token);
    if( ~isnan(num) )
      mpos(i) = num;
      i = i+1;
    end
  end
  
  mpos(i:end) = [];
  nmarkers = i-1;
  
  npheno = nfields - nmarkers;
  
  
  % marker names and phenotype names
  %  pnames = string( cell( npheno, 1 ) );
  %  mnames = string( cell( nmarkers, 1 ) );
  pnames = firstline( 1:npheno );
  mnames = firstline( npheno+1:end );

  % count number of individuals
  n = 0;
  s = fgetl( fid );
  while( s~=-1 )
    n=n+1;
    s = fgetl( fid );
  end
  
  frewind( fid); 
  s = fgetl(fid); 
  s = fgetl(fid);
  s = fgetl(fid);
    
  cnames = cell(n,1);
  pheno = zeros(n,npheno-1);
  geno = cell(n,nmarkers);
  
  for( i=1:n )
    s = fgetl(fid);
    [token,s] = strtok(s,dlm);
    cnames{i} = token;
    for( j=2:npheno )
      [token,s] = strtok(s,dlm);
      pheno( i, j-1 ) = str2num( token );
    end

    for( j=1:nmarkers )
      [token,s] = strtok(s,dlm);      
      geno{ i, j } = token;
    end
  end

  
  if( isempty(str2num(geno{1})) )
    % fprintf( 'String\n' )
    [nind,nmark]=size(geno);
    newgeno = zeros(nind,nmark);
    numgcodes = 0:(length(gcodes)-1);
    numgcodes(end) = gmiss;
    
    for( k = 1:length(numgcodes) )
      for( i=1:nind )
	for(j=1:nmark)
	  if(geno{i,j}==char(gcodes{k}) )
	    newgeno(i,j) = numgcodes(k);
	  end
	end
      end
    end  
  else
    [nind,nmark]=size(geno);
    newgeno = zeros(nind,nmark);
    for( i=1:nind )
      for( j=1:nmark )
	newgeno(i,j) = str2num(geno{i,j});
      end
    end
  end
  
  geno = newgeno;
  geno( find(geno==gmiss) ) = NaN;
  pheno( find(pheno==pmiss) ) = NaN;
  
  if( strcmp(casenames,'no') )
    y = [ str2num(char(cnames)) pheno ];
  else 
    y = pheno;
  end
  
  
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
 


function [delimiter,gcodes,gmiss,pmiss]=argchk( args )

  % default values
  delimiter = 'tab';
  gcodes = [];
  gmiss = 9;
  pmiss = -999;
  
  nvarargs = length( args );
  if( nvarargs>0 )
    nstep = 1;
    while( nstep <= nvarargs )
      argtype = args{nstep};
      switch argtype
       case 'delimiter'
	delimiter = args{nstep+1};
	nstep = nstep+2;
       case 'gcodes' 
	gcodes = args{nstep+1};
	nstep = nstep+2;
       case 'gmiss' 
	gmiss = args{nstep+1};
	nstep = nstep+2;
       case 'pmiss' 
	pmiss = args{nstep+1};
	nstep = nstep+2;
      end
    end
  end
  
	



function dlm = delimchk( delimiter )

  switch delimiter
   case 'tab'
    dlm = char(9);
   case 'space'
    dlm =  char(32);
   case 'comma'
    dlm = ',';
   case ','
    dlm = ',';
   case 'semi'
    dlm = ';'
   case ';'
    dlm = ';'
   case 'bar'
    dlm = ';'
   case '|'
    dlm = '|'
  end
  
