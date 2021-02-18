function [y,m,q] = makebc( varargin )
% MAKEBC Function to simulate data from multiple QTL model for backcross
%
% [Y,M]=MAKEBC(...)
% [Y,M,Q]=MAKEBC(...)  
% Y = phenotype vector
% M = marker data structure
% Q = QTL genotypes  
% Simulates data from a backcross.  The arguments should be in name and
% value pairs.  They are as follows (along with defaults).
% n: number of observations; default 100
% chromlen: vector of chromosome lengths; default 1  
% qpos: positions of the QTLs given as a matrix with two columns; the
%       first column is the chromosome number of the QTL, the second is
%       the position in Morgans of the QTL; default [ 1 0.5 ]
% mpos: marker position; can be scalar or a cell array of length equal to
%       the number f chromosomes; if scalar it assumes all markers are
%       equally spaced with the spacing equal to the scalar; if a cell
%       array, each cell is an array specifying the positions of the markers
%       on that chromosome; default 0.1
% mu: vector of QTL main effects
% mu2: matrix giving QTL two-way interaction effects; the first two columns
%      give the QTL numbers and the last one gives the effect size
% mu3: matrix giving QTL three-way interaction effects; the first three
%      columns give the QTL numbers and the last one gives the effect size
%
% See also: BCMODEL3, MAKEF2.
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.834 $ $Date: 2001/12/17 16:15:58 $	
  
% defaults
  n = 100; % number of individuals
  mu = 1; % the effect of the qtls
  qpos = [ 1 0.5 ]; % the first argument is the chromosome number and the
                    % second one is the position in *Morgans*
  chromlen = 1; % the chromosome lengths also declare the number of
                % chromosomes 
  mpos = 0.1; % can be a cell array of marker positions

  nargin = length(varargin);
  if( length(varargin)>0 )
    nstep = 1;
    while( nstep <= nargin )
      argtype = varargin{nstep};
      switch argtype
       case 'n'
	n = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'qpos'
	qpos = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'mu'
	mu = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'mu2'
	mu2 = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'mpos'
	mpos = varargin{ nstep+1 };
	nstep = nstep+2;
       case 'chromlen'
	chromlen = varargin{ nstep+1 };
	nstep = nstep+2;
       otherwise
	error( 'Cannot recognize option.' )
	nstep = nargin;
      end
    end
  end

  
  % number of chromosomes
  nchroms = length(chromlen);
  chroms = 1:nchroms;
  q = [];
  
  % marker positions
  if( iscell(mpos)==0 )
    spacing = mpos;
    mpos = cell(1,nchroms);
    for( i=1:nchroms )
      mpos{i} = 0:spacing:chromlen(i);
    end
  end
  
  % make marker data structure
  m = repmat( struct( 'cross', 'bc', 'chrid', 0, 'geno', [], ...
	      'mnames', [], 'mpos', [] ), 1, nchroms );
  
  % the error; QTL effects will be added
  y = normrnd(0,1,n,1);
  
  for( i=1:nchroms)
    cross = 'bc';
    qdx = find( qpos(:,1) == i );
    qposi = qpos(qdx,2);
    qposi = qposi(:);
    mposi = mpos{i};
    mposi = mposi(:);
    mui = mu(qdx);
    nqtli = length(qposi);

    % make marker names
    mnamesi = cell( length(mposi), 1 );
    for( j=1:length(mposi) )
      mnamesi{j} = strcat( 'c', num2str(i), 'mk', num2str(j) );
    end
          
    % list of all the positions, qtl and marker
    [ allpos, ord ] = sort( [ qposi; mposi ] );
    % make genotypes for backcross
    igeno = jointfiller ( repmat( NaN, n, 1 ), 0, allpos, 1, 'bc' );

    [ tmp, revord ] = sort( ord );
    
    if( nqtli>0 )
      qtlgeno = igeno(:,revord(1:nqtli) );

      y = y + bcmodel3(qtlgeno) * mui';
      q = [ q qtlgeno ];
    
      igeno(:,revord(1:nqtli))=[];
    end
    
    m(i).chrid = i;
    m(i).geno=igeno;
    m(i).mnames = mnamesi;
    m(i).mpos = mposi;
    m(i).chromlen = chromlen(i);
  end

  if(exist('mu2')>0)
    nint = size(mu2,1); % number of interactions
    for( i=1:nint )
      q1 = mu2(i,1);
      q2 = mu2(i,2);
      mu2i = mu2(i,3);
      y = y + prod_model( bcmodel3( q(:,q1) ), bcmodel3( q(:,q2) ) ) * mu2i;
    end
  end

  if(exist('mu3')>0)
    nint = size(mu3,1); % number of interactions
    for( i=1:nint )
      q1 = mu2(i,1);
      q2 = mu2(i,2);
      q3 = mu2(i,3);      
      mu3i = mu3(i,4);
      g1 = bcmodel3(q(:,q1));
      g2 = bcmodel3(q(:,q2));
      g3 = bcmodel3(q(:,q3));      
      y = y + prod_model( prod_model( g1, g2 ), g3 ) * mu3i;
    end
  end

