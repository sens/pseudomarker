function reportscan2( lod1, lod2, cv )
% REPORTSCAN2 Report pairs of loci deemed interesting
%
% REPORTSCAN2(LOD1,LOD2,CV)  
% LOD1 = output of MAINSCAN or MAINESTIMATE
% LOD2 = output of PAIRSCAN
% CV = 1x3 vector of critical values  
%      [ overall int coat ]
%      overall = overall critical value determined from PAIRPERMUTE2;
%                this should be on the LOD scale (log base 10) scale
%      int = chisq cutoff level for interaction
%      coat =   chisq cutoff level for coat tail
%  
% See also: FLAG, FLAG2, PERMUTEST2.  
 
  % cv(1) = cv(1)*log(10); % default is to have cutoff on the log 10 scale

  nchrom = length( lod1 ); % number of chromosomes
  cross = lod1(1).cross; % the type of cross being analyzed
  
  % set the penalties for the different kinds of crosses
  if ( cross == 'bc' )
    dfint = 1; % 4 df vs 3 df
    dfadd = 1; % 3 df vs 2 df
  elseif ( cross == 'f2' )
    dfint = 4; % 9 df vs 5 df
    dfadd = 2; % 5 df vs 3 df
  else
    error( 'Unrecognized cross type' )
  end

  % print the initial message
  fprintf( '\n\tThe following loci were flagged:\n\n' );
  fprintf( ...
      'Chrom\tcM1\tcM2\tLODfull\tLODint\tpval\tLODQ1\tpval\tLODQ2\tpval\n' );

  for( i=1:nchrom )
    for( j=i:nchrom )

      if( i==j)
	tmpint = triu( lod2(i,j).lod );
	[ mmmint, idx ] = max( tmpint(:) );
	tmpadd = tril( lod2(i,j).lod )';
	% mmmadd = max( tmpadd(:) );
	[s1,s2] = ind2sub( size(tmpint), idx );
	mmmadd = tmpadd(idx);
      else
	tmpint = lod2(i,j).lod;
	[ mmmint, idx ] = max( tmpint(:) );
	tmpadd = lod2(j,i).lod';
	% mmmadd = max( tmpadd(:) );
	[s1,s2] = ind2sub( size(tmpint), idx );
	mmmadd = tmpadd(idx);
      end	


      % does the max exceed limits?
      if( mmmint > cv(1) )
	
	% MI vs MA
	lodint = (mmmint - mmmadd);
	pint = 1-chi2cdf(2*lodint*log(10),dfint);
	% MA vs Q1
	lodq1 = mmmadd - lod1(i).lod(s1);
	pq1 = 1-chi2cdf(2*lodq1*log(10),dfadd);
	% MA vs Q2
	lodq2 = mmmadd - lod1(j).lod(s2);
	pq2 = 1-chi2cdf(2*lodq2*log(10),dfadd);
	
	% test for the interaction
	%lr = mmmint - mmmadd;
	if( (pint<cv(2)) | ( max(pq1,pq2)<cv(3) ) )
	  lab = strcat( int2str(lod1(i).chrid),':',...
			int2str(lod1(j).chrid), '\t', ...
			num2str(100*lod1(i).mpos(s1)), '\t', ...	    
			num2str(100*lod1(j).mpos(s2)), '\t', ...
 			num2str(mmmint,'%6.2f'), '\t', ...
			num2str(lodint,'%6.2f'), '\t', ...
			num2str(pint,'%0.4f'), '\t', ...
			num2str(lodq1,'%6.2f'), '\t', ...
			num2str(pq1,'%0.4f'), '\t', ...
			num2str(lodq2,'%6.2f'), '\t', ...
			num2str(pq2,'%0.4f'), '\n' );
	  fprintf( lab );
	end
      end
	
    end
    
  end
  

