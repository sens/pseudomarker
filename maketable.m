function [t,chi2,p]=maketable( c1, c2 )
% MAKETABLE Tabulate two columns of positive integers that may have
%           missing values

  c1 = c1+1;
  c2 = c2+1;
  
  a1 = find(isnan(c1));
  a2 = find(isnan(c2));  

  if( length(a1>0) )
    c1 = c1+1;
    c1(a1) = 1;
  end

  if( length(a2>0) )
    c2 = c2+1;
    c2(a2) = 1;
  end

  
  [t,chi2,p]=crosstab(c1,c2);