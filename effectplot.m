function effectplot( means, vars, labels, vscale )
% EFFECTPLOT Plots the effects and approximate 95% confidence intervals.
%  
% EFFECTPLOT(MEANS,VARS,LABELS )  
% MEANS = vector of posterior means
% VARS = vector of posterior variances
% LABELS = labels for the effects
%
% See also INTPLOT.
  
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
