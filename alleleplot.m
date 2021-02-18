function alleleplot( pheno, varargin )
% ALLELEPLOT Plot allelic effects or markers or pseudomarkers.
%         
% ALLELEPLOT(PHENO,MK)
% ALLELEPLOT(PHENO,MK,MNAME)    
% ALLELEPLOT(PHENO,MK,MNAME,LABELS)
%
% ALLELEPLOT(PHENO,MK1,MK2)
% ALLELEPLOT(PHENO,MK1,MK2,MNAMES)    
% ALLELEPLOT(PHENO,MK1,MK2,MNAMES,LABELS)
%
% PHENO = vector of phenotypes
% MK,MK1,MK2 = vector of marker (or pseudomarker) genotypes  
% MNAME,MNAMES = names of markers; default {'Marker'} or {'Marker1','Marker2'}
% LABELS = labels for genotypes; default {'A','B'} for backcross and
%          {'A','H','B'} for intercross
%
% If there are any missing phenotypes, those will be discarded. Whether
% the genotypes are markers or pseudomarkers is decided from the
% dimension of MK, MK1, or MK2.  The rows correspond to individuals and
% columns to imputations.  


% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.834 $ $Date: 2001/09/24 23:06:00 $	

  % number of arguments
  nargin = length( varargin );

  if( nargin==1 )
    mk = varargin{1};    
    intplot1( pheno, mk );
  elseif( isnumeric(varargin{2}) )
    mk1 = varargin{1};
    mk2 = varargin{2};
    if( nargin==2 )
      intplot2( pheno, mk1, mk2 );
    elseif( nargin==3 )      
      mnames = varargin{3};
      intplot2( pheno, mk1, mk2, mnames );
    elseif ( nargin==4 )
      mnames = varargin{3};      
      labels = varargin{4};      
      intplot2( pheno, mk1, mk2, mnames, labels );      
    end
  else
    mk = varargin{1};
    if( nargin==1 )
      intplot1( pheno, mk );
    elseif( nargin==2 )      
      mnames = varargin{2};
      intplot1( pheno, mk, mnames );      
    elseif ( nargin==3 )
      mnames = varargin{2};      
      labels = varargin{3};      
      intplot1( pheno, mk, mnames, labels );            
    end
  end
  
