function flaglod( lod1, lod2, cv )
% FLAGLOD Flag two-locus models based on one- and two-dimensional scans
%      based on likelihood ratios.
%
% FLAGLOD(LOD1,LOD2,CV)  
% LOD1 = output of MAINSCAN or MAINESTIMATE 
% LOD2 = output of PAIRSCAN
% CV = 1x3 vector of critical values  
%      [ overall int coat ]
%      overall = overall critical value determined from PAIRPERMUTE2;
%                this should be on the LOD scale (log base 10) scale
%      int = chisq cutoff level for interaction; usually 0.005
%      coat =   chisq cutoff level for coat tail; usually 0.005
%  
% See also: FLAGBF, PERMUTEST2.  
 
% Copyright 2000-2001: Saunak Sen
% Please cite: Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping", to appear in Genetics.  
%	$Revision: 0.832 $ $Date: 2001/09/24 20:54:47 $	
  
  % cv(1) = cv(1)*log(10); % default is to have cutoff on the log 10 scale

  nchrom = length( lod1 ); % number of chromosomes
  penalty = sqrt(lod1(1).n); % the Bayes factor per parameter penalty
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

  penaltyint = 0.5*chi2inv( 1-cv(2), dfint ); %
  penaltyadd = 0.5*chi2inv( 1-cv(3), dfadd ); % 

  % print the initial message
  fprintf( '\n\tThe following loci were flagged:\n\n' );
  fprintf( '\tLoci\tLODfull\tLOD\tpval\tcM1\tcM2\n' );

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
	
	% test for the interaction
	lr = (mmmint - mmmadd) * log(10);

	
	if( lr > penaltyint ) % if interaction important print message
	    lab = strcat( '\t', int2str(lod1(i).chrid),'x',...
			  int2str(lod1(j).chrid), '\t', ...
			  num2str(mmmint,'%6.2f'), '\t', ...
			  num2str(lr/log(10),'%6.2f'), '\t', ...
			  num2str(1-chi2cdf(2*lr,dfint),'%0.4f'), '\t', ...
			  num2str(100*lod1(i).mpos(s1)), '\t', ...	    
			  num2str(100*lod1(j).mpos(s2)), '\n' );
	    fprintf( lab );
	end
	
	%else % else check for coat tails

	% get the bayes factors for additive model
	%	if( i==j )
	%	  lr = mmmadd - max( lod1(i).lod );
	%	else 
	lr1 = mmmadd - lod1(i).lod(s1);
	lr2 = mmmadd - lod1(j).lod(s2);	    
	lr = log(10) * min( [ lr1 lr2 ] );
	%	end
	
	if( lr > penaltyadd ) % if additive effect exists
	  lab = strcat( '\t', int2str(lod1(i).chrid),'+',...
			int2str(lod1(j).chrid), '\t', ...
			num2str(mmmint,'%6.2f'), '\t', ...
			num2str(lr/log(10),'%6.2f'), '\t',...
			num2str(1-chi2cdf(2*lr,dfint),'%0.4f'), '\t', ...
			num2str(100*lod1(i).mpos(s1)), '\t', ...	    
			num2str(100*lod1(j).mpos(s2)), '\n' );
	  fprintf( lab );
	end
	
      end
      
    end
    
  end
