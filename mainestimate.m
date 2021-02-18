function onelod = mainestimate( y, x, mdata, fake, notmiss )
% MAINESTIMATE Estimate effects from a one-dimensional scan of genome.
%
% MAINESTIMATE(Y,X,MDATA,FAKE,NOTMISS)  
%  
% Y = phenotypes; can be a matrix
% X = covariates, use [] if none  
% MDATA = observed marker data structure  
% FAKE = pseudomarker data structure  
% NOTMISS = not missing index
%  
% The function returns a structure array with the following components:
%    N = number of individuals
%    CROSS = cross type
%    CHRID = chromosome number
%    BF = un-normalized Bayes factor
%    LOD = log of the sum of the weights at a given location; the BLOD  
%    MPOS = pseudomarker positions
%    EFFECTS = posterior means of the effects at given pseudomarker
% 	       locations
%    EFFVAR = posterior variances of the effects   
%
% The coding of the effects is as follows:
%    For the backcross, the effects are those corresponding to
%    substituting the 1 allele for the 0 allele.
%    For intercross, the additive and dominance effects are estimated.
%    The additive effect is the average effect of a single allelic
%    substitution.  The dominance effect is the average of the difference
%    of the level2-level1 and level1-level0.
%
%  See also F2MODEL, ONEESTIMATE, PLOTMAINESTIMATE.
  
  y = y(notmiss,:);
  if( ~isempty(x) )
    x = x(notmiss,:);
  end
  
  [ n, ntr ] = size(y);
  nchroms = length( fake );
  chroms = zeros(1,nchroms);
 
  % spacing between fakearkers
  spacing = fake(1).mpos(2) - fake(1).mpos(1);
  
  chromlen = zeros( 1, nchroms );
  pos = [];
  for( i=1:nchroms )
    chroms(i) = fake(i).chrid;
    chromlen(i) = mdata(chroms(i)).chromlen;
  end

  cross = fake(1).cross;
  onelod = repmat( struct( 'n', n, 'cross', cross, ...
			   'chrid', 0, 'bf', 0, 'lod', [], ...
			   'mpos', [], 'effects', [], 'effvar', [] ), ...
		   1, nchroms );
  
  for( i=1:nchroms )
    [mmm,vvv,thislod,bf] = oneestimate( y, x, ...
			     fake(i).igeno(notmiss,:,:), cross );
    onelod(i).chrid = chroms(i);
    onelod(i).bf = bf;
    onelod(i).effects = mmm;    
    onelod(i).effvar = vvv;        
    onelod(i).lod = thislod;
    onelod(i).mpos = fake(i).mpos;
    onelod(i).chromlen = mdata(i).chromlen;
  end
