function wt = binowt( y, x )
% BINOWT Calculate binomial weights on the log scale.  

  % size of the phenotype matrix
  [n,k] = size(y);
  % size of the covariate matrix
  [n,p] = size(x);
  
  % if there is a covariate, set the prior to something sensible (?)
  if( nargin==2 )
    m = ones(n,1); %
    b0 = zeros(p,1);
    b0(1) = logit(mean(y));
    %v0 = 0.0001*ones(1,p);
    v0 = 0.0001*ones(1,p);    
  end

  % the following lines were the original lines
  % get the estimates, deviance and unscaled covariance from the Bayesian
  % logistic regression function
  % [b,dev,v] = blogitfit(y,x,m,b0,v0);
  %  wt = dev - 0.5*log( det(v) / prod(v0) );

  % the following line could have been used but was found to be
  % numerically unstable
  %[dev,b,v] = logitfit(y,x,m);
  
  [b,dev,v] = blogitfit(y,x,m,b0,v0);  
  wt = dev;
