function qtl = panova( y, x, z, fake, loci, cv, model, obsdx )
% BACKWARD Perform model selection by dropping terms
% 

% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.834 $ $Date: 2002/01/09 21:16:45 $	

  
  % if the model is missing, assume no interactions
  if( nargin < 7 )
    model = [];
  end

p = 1;

while(p>cv)  
  qtl = makeqtl(fake,'chrid',loci(:,1),'mpos',loci(:,2));
  pa = panova(y,x,z,qtl);
  [p,i]=max(pa(:,6));
  
  loci = loci( [1:(i-2) i:end],:);
end


