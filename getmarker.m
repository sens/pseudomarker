function mm = getmarker( m, lod1, chroms )
% GETMARKER Get marker genotypes for maximum lod score on a particular
%           chromosome   
%
% MM=GETMARKER(M,LOD1,CHROM)
% MM = output marker genotype structure 
% M = marker genotype structure  
% LOD1 = LOD structure returned by MAINSCAN
% CHROMS = chromosome number


  if( nargin == 3 )
    nchroms = length(chroms);
  end

  if( nargin == 2 )
      nchroms = length(m);
      chroms = 1:nchroms;
  end


  
  n = size(m(1).geno,1);

  mm = repmat( struct( 'chrid', 0, 'geno', [], 'mnames', [], ...
			  'mpos', [], 'chromlen', 0 ), 1, nchroms );

  for( i = 1:nchroms )
    [ maxlod, idx ] = max(lod1(chroms(i)).lod);
    psmpos = lod1(chroms(i)).mpos(idx);
    [ mdist, idx ] = min( abs( m(chroms(i)).mpos - psmpos ) );
    % mnames = [ mnames m(chroms(i)).mnames(idx) ];
    mm(i).mnames = m(chroms(i)).mnames(idx);    
    % g(:,i) = m(chroms(i)).geno(:,idx);
    mm(i).geno = m(chroms(i)).geno(:,idx);    
    % mpos(i) = m(chroms(i)).mpos(idx);
    mm(i).mpos = m(chroms(i)).mpos(idx);
    mm(i).chromlen = m(chroms(i)).chromlen;
    mm(i).chrid = m(chroms(i)).chrid;    
  end


  