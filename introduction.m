%
%
%             An introduction to the PSEUDOMARKER programs
%
%                          Saunak Sen
%                   The Jackson Laboratory
%                         25 September 2001
%                 For Pseudomarker version 0.9
%
% This file is also a script file for Matlab, i.e. you can cut and paste
% this file into a Matlab window and it will run as a program.  As you can
% probably guess, everything that follows a % sign on a line is treated as a
% comment.
%
% The statistical theory behind these programs is discussed in the
% manuscript Sen and Churchill (2001) "A statistical framework for
% quantitative trait mapping," Genetics, 159:371-387.  For more
% information on the hypertension data please refer to Sugiyama
% et.al. (2001) "Concordance of Murine Quantitative Trait Loci for
% Salt-Induced Hypertension with Rat and Human Loci", Genomics, 71:70-77.
%
% The Pseudomarker programs are available at: 
%         www.jax.org/research/churchill/software/pseudomarker
% We also have an archive of QTL datasets that are available for the
% research community to download and analyze.  Point your browser to
%         www.jax.org/research/churchill/datasets/qtl/qtlarchive
% See also our web page: www.jax.org/research/churchill.
%
%
%                       HYPERTENSION DATA
%
% We will use the data collected by Dr. Bev Paigen's lab on salt-induced
% hypertension in mice to walk through an example data analysis.
%
% BACKGROUND ON THE DATA
%
% It is known that blood pressure in the C57BL/6J is raised by salt in the
% diet, whereas the A/J strain is not. Blood pressure measurements were
% obtained from 250 male mice from the backcross B6x(B6xA)F1. Initially only
% the extremes of the backcross population were typed.
%
% The original analysis of the data presented in Sugiyama et.al. (2001)
% reported the following QTL: two on chromosome 1, one each on chromosome 4,
% 5, 6, 15 and 11. There were three interacting loci: 1x5, 5x11, and 6x15.
%
% In our reanalysis we find strong evidence for five QTL, one each on
% chromosomes 1, 4, 6, 7, and 15. There are two pairs of interacting loci:
% 6x15 and 7x15. There was suggestive evidence for two QTL on chromosome and
% also on chromosome 4.


%                            FIRST STEPS
%
% Change directory to where the hypertension data has been stored.  The
% data must be in the pseudomarker format; for details please check our
% web page.
cd /home/ssen/qtl/data/hypertension
% Add to the search path so that Matlab can find the pseudomarker programs.
addpath /home/ssen/qtl/pseudo
% Check the version of Pseudomarker you are using.
pseudomarker
% For help on any function type "help functionname"; for example:
help impute


%                            READING DATA
%
% Read the data in first.
%
% Since the data is in the current directory
% "/home/ssen/qtl/data/hypertension" and has been formatted using the
% default convention, we can use this simple command to read in the data.
msgbox('Test', 'createmode','nonmodal')
[bp,mdata] = readdata;
% The phenotypes are in "bp" and the marker data structure is "mdata".
% The position of two markers was not known.  This is printed in a warning
% message.
%
% If the data is in a spreadsheet saved in "text" format, then it can
% also be read in using the "importdata" command.  
[bp,mdata,pnames]=importdata('hypertension-data.csv','delimiter','comma');
% In the above case, the file "hypertension-data.csv" was saved as a
% comma-delimited text file.  See the help for "importdata" for more
% details and options.

%                           BASIC PLOTS
%
% We do some exploratory data analysis to ensure data integrity; may help
% spot problems with the data or the data entry process.

% MAKE HISTOGRAM OF BLOOD PRESSURE.
%
figure; hist( bp )
title('Histogram of systolic blood pressure');
% To make histogram without colors use:
% bar(hist(bp,10),1,'w')


% An alternative way to do histograms is to use a butterfly plot; this is
% specially useful for comparing the cross data with the parentals.
figure; butterfly( {bp}, 1, 'Blood pressure' );
title('Butterfly plot of blood presuure')
% CHECK THE MISSING DATA PATTERN.
%
% What did the raw genotypes look like?
figure; plotgeno( mdata, 'genolabel', {'B6','Het','Missing'} )
% Notice that almost half the genotypes look missing.  The first 100 or
% so individuals look more completely genotyped than the rest.
% 
% Let us plot the genotype data again.  This time we will sort the
% individuals according to their blood pressure.
figure; plotgeno( mdata, 'sortby', bp, 'genolabel', ...
		  {'B6','Het','Missing'}, 'cap', ...
		  'Genotype pattern: mice sorted by blood pressure')
% We can now see the effects of the selective genotyping strategy that
% was used for the cross.  Notice how the top 50 and botton 50 mice are
% relatively completely genptyped, while there are gaps for the
% intermediate mice.
%
% If you stare at the figure hard enough you will be able to identify some
% QTLs; most clear is the one on chromosome 6; See also chromosomes 1 and 4
% Notice how the hets on chromosome 6 are associated with higher blood
% pressure.  On chromosomes 1 and 4, the B6 homozygous mice have higher
% blood pressure.

%                       IMPUTATION STEP

% CREATE PSEUDOMARKER DATA STRUCTURE
fake = impute( mdata );

% The default is to use 16 imputations spaced 10cM apart. The cross type is
% automatically determined.  If desired, we can use different number of
% imputations and pseudomarker density.


%                       EXPLORATORY SCANS

% PERFORM ONE-DIMENSIONAL SCAN
% 
% The function "mainscan" takes the phenotype vector (bp), covariates ([],
% none in our case), the imputed data structure (fake) to perform the
% one-dimensional genome scan.  You can also specify which subset of animals
% to analyze; by default all animals with no missing phenotypoes and
% covariates are analyzed.
lod1 = mainscan( bp, [], [], fake );

% PERFORM TWO-DIMENSIONAL SCAN
%
% This works the same way as the one-dimensional scan.
lod2 = pairscan( bp, [], [], fake );

% PLOT THE SCANS
% 
% Mainscan results plotted on the LOD scale by default; we will plot it on
% the proportion of variance scale.
figure;plotmainscan(lod1,'prop')
title( 'Mainscan with proportion of variance explained' );

% Pairscan results plotted on a proportion of variance scale; the lower
% triangle plots the proportion of variance explained by a two-QTL model
% while the upper triangle plots the proportion of variance explained by the
% interaction term over the additive term (inflated by a factor of 3)
figure;plotpairscan(lod2,'prop')
title( 'Pairscan with proportion of variance explained' );

% You can also plot the scan results on the LOD scale by using
figure;plotmainscan(lod1)
title( 'Mainscan with LOD score' );
figure;plotpairscan(lod2)
title( 'Pairscan with LOD score' );



%                       MODEL SELECTION
%
% PERMUTATION TESTS
%
% Brace yourself! This may take some time...  
%
% For illustration purposes we will perform the permutation tests on
% pseudomarkers with only one imputation.  This is not the recommended
% procedure, but will speed up the process of getting permutations
% thresholds coniderably.
%
% Make pseudomarkers with only one imputation.
fakesmall = impute(mdata,0.1,1);
% Perform the permutation tests for the mainscans.
maxlod1=permutest(bp,[],[],fakesmall,1000);
% Perform the permutation tests for the paircans.
maxlod2=permutest2(bp,[],[],fakesmall,40);

% Find out the threshold for the mainscan cutoff on the LOD scale.  This
% command will give us the 10%, 5% and 1% thresholds.
thresh1 = thresh( maxlod1 )
% Using the threshold from the mainscan we can flag the main loci.
reportscan(lod1,3)

% Do the same for the pairscans.
thresh2 = thresh( maxlod2 )
% For the pairscan thresholds, the first column will give the threholds
% for the full model as opposed to the null model.  The second column
% will give the thresholds for the interaction test.


% SELECTING A SET OF TWO-LOCUS MODELS
%
% The flagging routines help us select a small number of two-locus
% models. 
%
% Flagging based on the Bayes factor:
% The first two arguments are the mainscan and pairscan structure
% the last argument is a vector with three numbers [ cv int coat ]
% cv is the overall cutoff from the pairscan permutations (5.1) the number
%       5.1 is based on 100 pairscan permutations performed by the author;
%       in this example we have done only 10
%   int is the Bayes factor cutoff for interactions (0.5)
%   coat is the Bayes factor cutoff for coat-tail effects (0.5)
% Please play around with the last two numbers; they indicate how
% stringent you would like to be
flagbf( lod1, lod2, [ 5.1 0.5 0.5 ] )

% Flagging baesed on the likelihood ratio tests:
% This is similar to the Bayes factor version, except that the secondary
% testing is done not on the basis of Bayes factors, but likelihood ratio
% tests.
flaglod( lod1, lod2, [ 5.1 0.001 0.001 ] )

% If we prefer to flag loci using likelihood ratio tests, we can also use
% the following command.
reportscan(lod1,lod2, [ 5.1 0.001 0.001 ] )


%                    ANALYSES WITH SELECTED LOCI
%
% After a few loci have been selected, we may want to pull them out and
% investigate them further.  We select the pseudomarkers closest to the loci
% flagged.
selectfake = subsetgeno( fake, 'chrid', [ 1 4 6 15 ], ...
			 'mpos', [ 70 30 70 20 ] );

% PLOTTING ALLELIC EFFECTS
% Next we plot the allelic effects of selected loci using the function
% "alleleplot".
figure;
subplot(2,2,1)
alleleplot( bp, selectfake(1).igeno, 'Chr1QTL', {'B6','Het'} )
subplot(2,2,2)
alleleplot( bp, selectfake(2).igeno, 'Chr4QTL', {'B6','Het'}  )
subplot(2,2,3)
alleleplot( bp, selectfake(4).igeno, 'Chr15QTL', {'B6','Het'} )
subplot(2,2,4)
alleleplot( bp, selectfake(3).igeno, selectfake(4).igeno, ...
	    {'Chr6QTL','Chr15QTL'}, {'B6','Het'}  )
% Note that all the QTL on chromosomes 1, 4 and 15 have strong main
% effects.  The interaction between loci on chromosomes 6 and 15 can be
% seen since the locus on chromosome 15 has an effect only when the locus
% on chromosome 6 is heterozygous.

% TYPE III ANALYSIS
% We can also do a "Type III" analysis using the general "scan" function.
% First we will fit the largest model with main effects on chromosomes 1,
% 4, 6 and 15, and a 6x15 interactions.  Then we will drop each of these
% and see how much the LOD score drops as a result.

% Fit full model:
[lod0,bf0] = scan( bp, [], [], selectfake, [ 1 4 6 15; 1 1 1 1 ], ...
		   struct( 'twoint', [ 3 4 ] ) );
% Drop QTL on chr1:
[lod1,bf1] = scan( bp, [], [], selectfake, [ 4 6 15; 1 1 1 ], ...
		   struct( 'twoint', [ 2 3 ] ) );
% Drop QTL on chr4
[lod4,bf4] = scan( bp, [], [], selectfake, [ 1 6 15; 1 1 1 ], ...
		   struct( 'twoint', [ 2 3 ] ) );
% Drop 6x15 interaction
[lod6x15,bf6x15] = scan( bp, [], [], selectfake, [ 1 4 6 15; 1 1 1 1 ] );
% Drop drop QTLs on chr6 and chr15
[lod6n15,bf6n15] = scan( bp, [], [], selectfake, [ 1 4; 1 1 ] );

% Compare models to the largest model:
comparison_lods = [ lod0-lod1 lod0-lod4 lod0-lod6x15 lod0-lod6n15 ]
% Calculate nominal p-values associated with those terms
comparison_pvalues = [ 1-chi2cdf(2*log(10)*(lod0-lod1),1)...
  1-chi2cdf(2*log(10)*(lod0-lod4),1)...
  1-chi2cdf(2*log(10)*(lod0-lod6x15),1)...
  1-chi2cdf(2*log(10)*(lod0-lod6n15),1) ]
% We can see that all terms in the model do contribute and we keep them.
% So our selected model has QTL on chromosomes 1, 4, 6, and 15.  The last
% two loci interact.



% SAVING RESULTS 
%
% At this point you may want to save the work we have done so far.
% The command
% save HYPERTENSION
% will save the variables in the current workspace in a file called
% HYPERTENSION.mat in our working directory.




%                      FINE MAPPING
% 
% PREPARATORY STEPS
%
% On to fine mapping now.  We will concentrate on the chromosomes that are
% of interest to us. First we make a pseudomarker data structure at a
% higher resolution (5cM) and larger number of imputations (64) than the
% default.  We also concentrate on chromosomes 1, 4, 6 and 15 which look
% most promising on the basis of the one-way and two-way scans.
fakebig = impute( mdata, 0.05, 64, [ 1 4 6 15 ] );
% Mainscan at higher resolution.
lodbig1 = mainscan( bp, [], [], fakebig );
% Pairscan at higher resolution.
lodbig2 = pairscan( bp, [], [], fakebig );
% Plot the scan results.
figure; plotmainscan(lodbig1);
title( 'Closeup Mainscan with LOD score' )
figure; plotpairscan(lodbig2);
title( 'Closeup Pairscan with LOD score' )

% LOCALIZING ONE QTL
%
% Posterior distribution of the location of QTL on chromosome 1.
figure; plotmainlocalize( lodbig1(1) )
title('Posterior distribution of single QTL on chromosome 1')
% You may try any others you want.

% LOCALIZING TWO QTL
%
% Are there two QTL on chromosome 1?  If so where are they?  The following
% command will make 50%, 95% and 99% confidence regions assuming there are
% two QTL on chromosome 1.
figure; plotpairlocalize( lodbig2(1,1), [ 0.5 0.95 0.99 ] );
title('Posterior distribution of two QTL on chromosome 1')
% The confidence region seems quite diffuse.  This indicates that that there
% is not enough evidence in favor of two QTL on chromosome 1, although
% the evidence is suggestive.


%                           ESTIMATING EFFECTS
%
% ESTIMATE SCANS
%
% We can perform a one-dimensional "estimate scan"; this is a genome scan
% with the estimates of the QTL effects on each chromosome.
lod3 = mainestimate( bp, [], mdata, fake, 1:250 );
figure; plotmainestimate( lod3 );
subplot(2,1,1)
title( 'Estimate scan with proportion of variance explained' )

% The procedure for estimating the effects of the QTL are somewhat tedious;
% this will change in future versions of the code.  The functions for
% estimation are called "oneestimate", "twoestimate", "twoestimate2",
% "threeestimate", "threeestimate2" and "threeestimate3".  Their use is
% illustrated below.

% SINGLE QTL ON A CHROMOSOME
%
% Estimate effect of QTL on chromosome 4
igeno = fake(4).igeno;
[mm4,vv4] = oneestimate( bp, [], igeno, 'bc' )
% Estimate effect of QTL on chromosome 1
igeno = fake(1).igeno;
[mm1,vv1] = oneestimate( bp, [], igeno, 'bc' ) 
% The first output argument has the estimates of the effects, the second
% element has the variances for the estimates the effects are coded in
% factorial notation which means that mm1(1) is always the grand mean and
% mm1(2) is the half of the effect of an allelic substitution.

% EFFECT PLOTS
%
% Plot the main effects of chromosome 1 4 6 15
% Estimate effects first.
igeno = fake(6).igeno;
[mm6,vv6] = oneestimate( bp, [], igeno, 'bc' )
igeno = fake(15).igeno;
[mm15,vv15] = oneestimate( bp, [], igeno, 'bc' )
% Effect plot
figure;effectplot( [ mm1(2) mm4(2) mm6(2) mm15(2) ], ...
		   [ vv1(2) vv4(2) vv6(2) vv15(2) ], ...
		   {'chr1', 'chr4', 'chr6', 'chr15'} );
title('Estimated main effects of QTL on chromosomes 1, 4, 6, and 15')

% TWO QTL ON A CHROMOSOME OR DIFFERENT CHROMOSOMES 
%
% Let us take the case of QTL on chromosome 6 and 15; they seem to interact.
igeno6 = fake(6).igeno;
igeno15 = fake(15).igeno;
% Estimate the additive effects for QTL on different chromosomes.
[mm6a15,vv6a15] = twoestimate2( bp, [], igeno6, igeno15, 'bc', 'a+b' )
% Estimate the interaction effects as well.
[mm6x15,vv6x15] = twoestimate2( bp, [], igeno6, igeno15, 'bc', 'a*b' )
% Plot the main effects as well as the interaction effect.
figure;effectplot( [ mm1(2) mm4(2) mm6(2) mm15(2) mm6x15(3)], ...
		   [ vv1(2) vv4(2) vv6(2) vv15(2) vv6x15(3)], ...
		   {'chr1', 'chr4', 'chr6', 'chr15' 'chr6x15'} );
title('Estimated effects of QTL in final model')

% To estimate the effects of two QTL on the same chromosome, such as
% chromosome 1, use the following:
igeno1 = fake(1).igeno;
[mm1a1,vv1a1] = twoestimate( bp, [], igeno1, 'bc', 'a+b' )
figure;effectplot( mm1a1(2:3), vv1a1(2:3), {'chr1a', 'chr1b'} );
title('Estimated effects of two QTL on chromosome 1')
% This figure is another piece of evidence that there is not sufficient
% evidence for two QTL on chromosome 1.






