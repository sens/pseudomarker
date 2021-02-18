function reportscan( varargin )
% REPORTSCAN Report loci deemed interesting
%
% REPORTSCAN(LOD1,CV) returns REPORTSCAN1(LOD1,CV)
% REPORTSCAN(LOD1,LOD2,CV) returns REPORTSCAN2(LOD1,LOD2,CV)   
%
% See also: REPORTSCAN1, REPORTSCAN2.
 
nargin = length( varargin );

if( nargin == 2 )
  reportscan1( varargin{1}, varargin{2} );
elseif( nargin == 3 )  
  reportscan2( varargin{1}, varargin{2}, varargin{3} );
end