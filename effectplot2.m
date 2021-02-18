function effectplot2( means, vars, labels, vscale )
% EFFECTPLOT2 Function to plot estimated effects and confidence intervals.
%  
% EFFECTPLOT2(MEANS,VARS,LABELS)  
% EFFECTPLOT2(MEANS,VARS,LABELS,VSCALE)    
% MEANS = vector of posterior means
% VARS = vector of posterior variances
% LABELS = labels for the effects
% VSCALE = scalar giving height at which a dotted horizontal line
%          is added; default height is 0
% See also EFFECTPLOT.
  
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.831 $ $Date: 2001/07/09 20:29:55 $	
  
  s = sqrt(vars);
  m = length(means);
  
%  figure
  hold on
  if( nargin~=4 )
    plot( [ 0 m+1 ], [0 0], ':')
  else
        plot( [ 0 m+1 ], [vscale vscale], ':')
  end
  
  for( i=1:m )
    plot(i,means(i),'.');
    plot([i i], [ means(i) + 2*s(i) * [-1 1] ] );
  end

  hold off
  set( gca, 'XTick', [ 1:m ]);
  set( gca, 'XTickLabel', labels );
