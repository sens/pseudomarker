function prop = missingprop( fake )
% MISSINGPROP Estimate the proportion of missing data
%
% MISSINGPROP(FAKE) calculates the proportion of missing information by
% calculating the average variance of the imputed genotypes in the structure
% FAKE and comparing it to the equilibrium variance of the genotypes.  The
% resulting output is a structure array with the following components
%   CROSS = cross type
%   CHRID = chromosome id  
%   MISSPROP = missing proportion
%   MPOS = pseudomarker positions  


  switch fake(1).cross
   case 'bc'
    basevar = 0.25;
    basemean = 0.5;
   case 'f2'
    basevar = 0.5;
    basemean = 1;
   otherwise
    error( 'Cross type not recognized.' )
  end

  nchroms = length(fake);
  prop = repmat( struct( 'cross', fake(1).cross, 'chrid', [], ...
			 'missprop', [], 'segdist', [], 'mpos', ...
			 [] ), 1, nchroms );
  
  
  for( i=1:nchroms )
    
    igeno = fake(i).igeno;
  
    sz = size( igeno );
    a = permute( igeno, [ 3 1 2 ] );
    b = reshape( a, [ sz(3) prod( sz(1:2) ) ] );
    c = var( b );
    d = reshape( c, sz(1), sz(2) );
    f = mean( b );
    g = reshape( f, sz(1), sz(2) );

    p = mean(d)/basevar;
    s = ( mean(g) - basemean );
    
    prop(i).chrid = fake(i).chrid;
    prop(i).missprop = p;
    prop(i).segdist = s;    
    prop(i).mpos = fake(i).mpos;
     
  end
  
