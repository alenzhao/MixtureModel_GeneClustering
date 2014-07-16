***************************************************
* GITTAR: GIbbs sampler for finding Transcription *
* factor TARget genes				  *
* by Li Shen, UCSD				  *
* 09-04-2007 Ver.1.0				  *
***************************************************

I.   Overview:

GITTAR employs a Gibbs sampler to estimate each gene's probability 
of being transcription factor (TF) target using both sequence and 
ChIP-chip log-ratios.

II.  Input:

The input to GITTAR consists of at least three files: 
1. A file contains the binding consensus motif of the TF. IUPACs are 
allowed to use here.
2. A file contains each gene's ChIP-chip ratio in TAB delimited format
3. A file contains each gene's promoter sequence in TAB delimited format

e.g. "./gittar -m YY1.motif -e YY1.expr -s YY1.seq". This says that YY1's 
motif, ChIP-chip ratios and promoter sequences are stored in the above 
three files.

The consensus motif is a short segment of DNA sequence describing the binding
specificity of a transcription facotr's binding site. e.g. CCATNTT describes
the binding motif of transcription factor YY1. GITTAR makes an assumption that
the DNA sequence describing the binding specificity of a TF should be around 21
bps. Therefore the core motif should be about 6-8 bps in length. GITTAR will
then concatenate a flanking region of length 7 to both sides of the core motif
and calculate the PSFM for the extended motif from data. Although all these numbers can
be specified as parameter values in GITTAR, it is suggested to follow the above
guidelines. Using too long or too short core motifs is not recommended.
Likewise, too long flanking regions are not recommended to use.

Note: The default input gene IDs to GITTAR are Unigene IDs. Using different 
types of IDs won't influence GITTAR running. But it will make gene annotation 
non-workable. You must follow the gene annotation format in "Hs.data.annot" to 
make it work!!!

GITTAR assumes the target genes' have higher log-ratios than background mean. 
If you want to apply GITTAR to experiments other than ChIP-chip, like RNAi
where the deletion of TF causes it's target genes repressed, you have to do
some work to reverse the sign of your input log-ratios.

III. Output:

The output of GITTAR are two files describing both the results and parameter
settings:
1. All parameter settings and each gene's name, probability and other 
information
2. The position specific frequency matrix (PSFM) of both target and 
non-target genes, including the ratios of target/non-target matrices

You can specify the names of the output files, or GITTAR will automatically 
generate them for you.

IV.  Additonal parameters:

There are additional parameters which may affect the results of GITTAR. 
If you don't specify them, GITTAR will use the default values which are usually 
suitable for most cases.

The optional parameters are explained here:

[-mis] mismatches allowed to consensus motif(Default=1). This parameter is 
useful to find candidate genes.

[-a] annotation file. See "Hs.data.annot" for format.

[-o] output(gene_and_other_info PSFM). e.g. "-o hyy1.info hyy1.log". This 
specify "hyy1.info" for all genes' information and "hyy1.log" for PSFM.

[-v] verbose mode. Switch on to make GITTAR noisy: it will output various 
information during GITTAR running. Default=off.

[-h] hyperprio parameters(hyper_samples=50 mu=2.0 sigma^2=1.0
flanking_length=7). e.g. "-h 50 2.0 1.0 7". This says that the 
number of hypothetical samples equal 50; hypothetical mean equal 2.0;
hypthetical sigma square equal 1.0; flanking region equal 7; You don't have 
to specify all of them. You may only set the first a few ones. 
e.g. "-h 30" sets the number of hypothetical samples to 30 and use the default 
values for the rest parameters.

[-n] DO NOT use prio information(overide hyperprio parameters). Turn on this 
switch will shut off prio information and will make GITTAR run in
maximum-likelihood fashion. GITTAR was designed in Bayesian statistics.

[-g] Gibbs sampler parameters(starting_point=10 burn-in=100 window=100). e.g. 
"-g 10 100 100". This will set GITTAR to run from 10 random starting points 
with a burn-in period of 100 iterations and a window size of 100 iterations. 
The final results of GITTAR are based on the average of the samples in the 10 
windows of 100 size each.

[-misc] Miscellaneous parameters(stds_of_intial_target=4 skip_core_region=0
random_background=0 threshold_for_target=0.5 stds_enfore=1). e.g. "-misc 4 0 0 
0.5 1". This sets the initial target genes to be 4 standard deviations above 
mean; core motif region is considered in calculating sequence probability; DO 
NOT choose a random binding site for genes with ratios below mean; threshold
for identifying a gene as target; number of standard deviation to force the
means of target and background to separate.
 


