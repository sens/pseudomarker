function table = panova( y, x, z, qtl, model, obsdx )
% PANOVA Function to perform a Type III analysis of a QTL model
% 
% TABLE = PANOVA(Y,X,Z,QTL)    
% TABLE = PANOVA(Y,X,Z,QTL,MODEL)  
% TABLE = PANOVA(Y,X,Z,QTL,MODEL,OBSDX)
%
% Y = vector of trait values  
% X = matrix of additive covariates, use [] if none
% Z = matrix of interacting covariates, use [] if none  
% FAKE = pseudomarker data structure; output of IMPUTE
% QTL = qtl data structure; output of MAKEQTL
% MODEL = model of the QTLs; structure with fields twoint and threeint.
%         eg. model= struct( 'twoint', [ 1 2 ], 'threeint', [ 1 2 3 ] );
%         default is empty
% OBSDX = observed data index; default is the individuals with no missing
%         phenotypes 
%  
% The function will output an ANOVA table for doing Type III analysis for a
% given QTL model.  The QTLs have to be selected using the function MAKEQTL.
% The model is selected via the model function as in SCAN.  If missing, the
% model will assume that all QTL are additive.
% 
% If there are no interactions, the function will drop each of the main
% effect QTLs one by one.  If there are any interactions, the function will
% drop the interactions and any QTL not involved in the interactions.  The
% LOD score corresponding to each model is returned.
%  
% This LOD score is used for subsequent calculations.  The difference in the
% LOD scores between models is used to also calculate a proportion of
% variance explained by the model and the proportion of variance explained
% by the term dropped.  The p-values are based on a chi-squared
% approximation.
%  
% See also: ONESCAN, TWOSCAN, THREESCAN, SCAN, MAKEQTL, WTAVERAGE.

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.834 $ $Date: 2002/01/09 21:16:45 $	

  
  % if the model is missing, assume no interactions
  if( nargin < 5 )
    model = [];
  end
  
  % if the qtl argument is missing, assume one qtl on each chromosome 
  % if( nargin < 4 )
  %  qtl = zeros( 2, length(fake) );
  %  for( i=1:length(fake) )
  %    qtl( 1:2, i ) = [ fake(i).chrid; 1 ];
  %  end
  % end
  
  % make a copy of the qtl structure
  qtlcopy = qtl;
  % number of qtl
  nqtl = length(qtl);
  
  % remove the qtlid field in qtl by using tmpqtl
  % this is necessary so that we can use qtl in the scan function
  % tmpqtl = rmfield(qtl(1),'qtlid');
  for( i=1:nqtl)
    qtl(i).chrid = qtl(i).qtlid;
    tmpqtl(i) = rmfield(qtl(i),'qtlid');
  end
  qtl = tmpqtl;
  % cross type
  cross = qtl(1).cross;
  
  
  % get a list of models to compare; this is the hard part
  % algorithm:
  % if only additive effects, drop each term one by one
  % if any interactions, drop the interaction then drop the additive
  % effects one by one; if dropping a main effect, drop the interaction
  % also
  
  % make a list of models and qtls
  % if the twoint field in the model does not exist, then there are no

  % two-way interactions, otherwise the number of interactions is equal
  % to the number of rows in twoint
  
  % if the list of models is empty then the number of two-way
  % interactions is zero, else, if it is not empty, see if the twoint
  % field is there or not
  if( ~isempty(model) )
    if( isfield(model,'twoint') ) 
      twoint = model.twoint;
      ntwoint = size(twoint,1);
    else
      ntwoint = 0;
      twoint =  [];
    end
  else
    ntwoint = 0;
    twoint = [];
  end
  
  % the number of terms to consider
  % at this point we consider only dropping all interactions or all main
  % main effects one by one

   % we drop all terms in the model one by one
 
   % vector of qtl id's
   qtlid = zeros(1,nqtl);
   for( i=1:nqtl )
     qtlid(i) = qtlcopy(i).qtlid;
   end
   
   % which are the interacting qtl
   intqtlid = unique(twoint(:)');
   % which are the non-interacting qtl
   nonintqtlid = setdiff( qtlid, intqtlid );
   nnontwoint = length(nonintqtlid);
   
   % CHANGE 
   nterms = ntwoint + nnontwoint + 1;
   % Change nterms = ntwoint + nqtl + 1;

  lod = zeros( 1, nterms );
  bf = zeros( nterms, 1 );  
  % accordingly we set the degrees of freedom of each term
  if( cross == 'bc' )
    df = ones(nterms,1);
  elseif( cross == 'f2' )
    df = zeros(nterms,1);
    % CHANGE
    df(2:1+nnontwoint,1) = 2;
    df(2+nnontwoint:end,1) = 4;      
    % Change
    % df(2:1+nqtl,1) = 2;
    % df(2+nqtl:end,1) = 4;      
  end
  
  % first we note the full model; then we make models that drop each term
  % one by one
  
  % cell array that will carry the model information
  terms = cell(nterms,2);
  % names of the terms
  termnames = cell(nterms,1);
  termnames{1} = 'Total';
  % mother model
  % mother model has all the additive qtl
  terms{1,1} = 1:nqtl;
  % if there are any interactions, mother model has all the two-way
  % interactions 
  if(ntwoint>0)
    terms{1,2} = twoint;
  end
  
   % -- CHANGE begin ---
   if(nnontwoint>0)
     % make terms if there are no interactions
     for( j=1:nnontwoint )
       tmpqtl = 1:nnontwoint;
       termnames{j+1} = qtlcopy(nonintqtlid(j)).mnames{1};
       tmpqtl(j) = [];
       tmpqtlid = union( qtlid(nonintqtlid(tmpqtl)), intqtlid );
       terms{j+1,1} = tmpqtlid;
       tmptwoint = twoint;
       if( ~isempty(twoint) )
	 for( k = 1:length(tmpqtlid) )
	   tmptwoint(twoint==tmpqtlid(k)) = k;
	 end
       end
       terms{j+1,2} = tmptwoint;
     end
   end
   % make the terms if there are interactions
   if(ntwoint>0)
     for( j=1:ntwoint )
       tmptwoint = twoint;
       terms{1+nnontwoint+j,1} = qtlid;
       q1name = qtlcopy(tmptwoint(j,1)).mnames{1};
       q2name = qtlcopy(tmptwoint(j,2)).mnames{1};     
       termnames{j+1+nnontwoint} = strcat( q1name, '.', ' ', q2name );
       tmptwoint = twoint;
       tmptwoint(j,:) = [];
       terms{1+nnontwoint+j,2} = tmptwoint;
     end
   end
   
   % -- CHANGE end ----
   
% $$$    %%%% -- CHANGE BEGIN ----
% $$$    % drop main effects
% $$$    for( j=1:nqtl )
% $$$      tmpqtl = 1:nqtl;
% $$$      termnames{j+1} = qtlcopy(qtlid(j)).mnames{1};
% $$$      tmpqtl(j) = [];
% $$$      tmpqtlid = union( qtlid(qtlid(tmpqtl)), intqtlid );
% $$$      terms{j+1,1} = tmpqtlid;
% $$$      tmptwoint = twoint;
% $$$      if( ~isempty(twoint) )
% $$$        for( k = 1:length(tmpqtlid) )
% $$$  	tmptwoint(twoint==tmpqtlid(k)) = k;
% $$$        end
% $$$      end
% $$$      terms{j+1,2} = tmptwoint;
% $$$    end
% $$$  
% $$$    % make the terms if there are interactions
% $$$    if(ntwoint>0)
% $$$      for( j=1:ntwoint )
% $$$        tmptwoint = twoint;
% $$$        terms{1+nqtl+j,1} = qtlid;
% $$$        q1name = qtlcopy(tmptwoint(j,1)).mnames{1};
% $$$        q2name = qtlcopy(tmptwoint(j,2)).mnames{1};     
% $$$        termnames{j+1+nqtl} = strcat( q1name, '.', ' ', q2name );
% $$$        tmptwoint = twoint;
% $$$        tmptwoint(j,:) = [];
% $$$        terms{1+nqtl+j,2} = tmptwoint;
% $$$      end
% $$$    end
% $$$  
% $$$   %%%% -- CHANGE END ----

  % step through the terms
  if(nterms==2)
    % model information
    qtlmodel = [];
    qtlmodel(1,:) = terms{1,1};
    % make the model information from the terms
    qtlmodel(2,:) = ones(1,length(terms{1,1}));
    % make a structure from the information one has
    qtltwoint = struct( 'twoint', terms{1,2} );
    % run the scan function
    [ tmplod, bf ] = scan(y, x, z, qtl(terms{1,1}), qtlmodel, ...
			     qtltwoint );
    % lod(1) = log(bf(1))/log(10);
    % since the SCAN function now returns log of bf
    lod(1) = bf;    
    obsdx = phenoobs([y x z]);
    n = length(obsdx);
    [rawss0,df0] = rss( y(obsdx), ones(n,1) );
    if(~isempty(x))
      [rawss1,df1] = rss( y(obsdx), [ones(n,1) x(obsdx,:)] )  ;
      lod(2) = -(n/2)*log(rawss1/rawss0)/log(10);
    else
      lod(2) = 0;
    end
  else
  
    for( i=1:nterms )
      % model information
      qtlmodel = [];
      qtlmodel(1,:) = terms{i,1};
      % make the model information from the terms
      qtlmodel(2,:) = ones(1,length(terms{i,1}));
      % make a structure from the information one has
      qtltwoint = struct( 'twoint', terms{i,2} );
      % run the scan function
      [ tmplod, bf(i) ] = scan(y, x, z, qtl(terms{i,1}), qtlmodel, ...
			       qtltwoint );
      % lod(i) = log(bf(i))/log(10);
      % since the SCAN function now returns log of bf
      lod(i) = bf(i);      
    end
  end
    
  % create the observed data index
  obsdx = phenoobs( [ y x z ] );
  
  px = 0;
  pz = 0;
  % calculate df for the full model
  if( ~isempty(x) )
    [n,px] = size(x);
  end
  if( ~isempty(z) )
    [n,pz] = size(z);
  end
  
  % rank of the covariates
  if( px*pz > 0 )
    r = rank( [ x z ] );
  else
    r = px + pz;
  end

  % total number of observations
  n = length(obsdx);
  % total df
  totaldf = n-r-1;

  % make anova table
  table = zeros( nterms, 4 );
  
  % first column - df
  table(:,1) = df;
  table(1,1) = totaldf;
  % second column - lods
  table(:,2) = lod';
  % third column - difference of the lods
  table(:,3) = lod(1)-lod';
  % fourth column - model proportion variance explained
  table(:,4) = 100 * lod2propvar( lod', n );
  
  % fifth column - adjusted proportion variance
  table(:,5) = table(1,4) - table(:,4);
  % sixth column - pvalues
  table(1,6) = NaN;
  table(2:end,6) = 1-chi2cdf(2*log(10)*table(2:end,3),df(2:end));

  termnames = strjust( strvcat( 'Terms', strvcat( termnames ) ), 'right' );
  col1 = strjust( strvcat('df', strvcat( ...
      int2str(table(:,1) ) ) ), 'right' );
  col2 = strjust( strvcat('LOD', strvcat( num2str(table(:,2),'%0.2f' ) ) ), ...
		  'right' );
  col3 = strjust( strvcat('LOD-diff', ' ', strvcat( ...
      num2str(table(2:nterms,3),'%0.2f' ) ) ), 'right' );
  col4 = strjust( strvcat('Model %var', ...
			  strvcat( num2str(table(:,4),'%0.2f' ) ) ), 'right' );
  col5 = strjust( strvcat('Adj %var', ' ', strvcat( ...
      num2str(table(2:nterms,5),'%0.2f' ) ) ), 'right' );
  col6 = strjust( strvcat('p-value', ' ', strvcat( ...
      num2str(table(2:nterms,6),'%0.6f' ) ) ), 'right' );
  delim = repmat( '  ', nterms+1, 1);

  Anovatable = [ termnames delim col1 delim col2 delim col3 delim ...
		 col4 delim col5 delim col6 ] 





