function theta = makecorr( fake, notmiss )
% MAKECORR Calculate correlations between pseudomarkers
%

  % number of individuals
  if( nargin == 1 )
    n = size( fake(1).igeno, 1 );
    notmiss = 1:n;
  else
    n = size(notmiss,1);
  end
  

  % spacing between fake markers
  spacing = fake(1).mpos(2) - fake(1).mpos(1);
  nchroms = length( fake ); % number of chromosomes
  % vector for the length of chromosomes
  chromlen = zeros( 1, nchroms ); 
  
  % position vector for plot
  pos = [];
  % number of fakearkers
  npseudo = [];
  chroms = zeros(1,nchroms);
  
  % get all the chromosome lengths and the number of fakearkers in
  % each chromosome
  for( i=1:nchroms )
    chroms(i) = fake(i).chrid;
    chromlen(i) = fake(i).chromlen;    
    npseudo = [ npseudo length(fake(i).mpos) ] ;
  end

  cross = fake(1).cross;
  theta = repmat( struct( 'n', n, 'cross', cross, ...
			  'chrid1', 0, 'chrid2', 0, 'cor', [], ...
			  'chromlen1', 0, 'chromlen2', 0,...
			  'mpos1', [], 'mpos2', [] ), nchroms, nchroms );
  
  for( i=1:nchroms )
    for( j=i:nchroms )

      %fprintf( '(%d,%d)', fake(i).chrid, fake(j).chrid );
      %if(j==nchroms)
	%fprintf('\n');
      %end
      
      %[ fake(i).chrid fake(j).chrid ]
      if(i==j)
	cor = corrcoef( fake(i).igeno(:,:,1) );
	theta( i, i ).cor = cor;
	theta( i, i ).chromlen1 = fake(i).chromlen;	
	theta( i, i ).chromlen2 = fake(i).chromlen;		
	theta( i, i ).mpos1 = fake(i).mpos;
	theta( i, i ).mpos2 = fake(i).mpos;	
	theta( i, i ).chrid1 = fake(i).chrid;
	theta( i, i ).chrid2 = fake(i).chrid;	
      else
	a = fake(i).igeno(:,:,1);
	b = fake(j).igeno(:,:,1);
	cor = corrcoef( [ a b ] );
	cor = cor( 1:size(a,2), size(a,2)+1:end );
	theta( i, j ).cor = cor;
	theta( j, i ).cor = cor';	
	theta( i, j ).chromlen1 = fake(i).chromlen;	
	theta( i, j ).chromlen2 = fake(j).chromlen;		
	theta( j, i ).chromlen1 = fake(j).chromlen;	
	theta( j, i ).chromlen2 = fake(i).chromlen;		
	theta( i, j ).mpos1 = fake(i).mpos;
	theta( i, j ).mpos2 = fake(j).mpos;	
	theta( j, i ).mpos1 = fake(j).mpos;
	theta( j, i ).mpos2 = fake(i).mpos;	
	theta( i, j ).chrid1 = fake(i).chrid;
	theta( i, j ).chrid2 = fake(j).chrid;	
	theta( j, i ).chrid1 = fake(j).chrid;
	theta( j, i ).chrid2 = fake(i).chrid;	
      end
      
      
    end

  end
  

