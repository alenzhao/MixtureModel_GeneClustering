#include "mysub.h"
#include <sstream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include "CmdLine.h"

using namespace std;

unsigned BURN_IN = 100;	// iterations required for burn-in.
unsigned WINDOW = 100;	// number of iterations used to calculate probabilities.
unsigned STA_PNT = 10;	// number of random starting points.
unsigned ENFORN = 1;	// number of standard deviations that should be enforced between target mean and background mean.
unsigned INIV = 4;		// Initial labels based a number of standard deviations.
double THRLD = 0.5;		// threshold for target probability.
int SKIPCORE = 0;		// skip core motif region when calculate sequence likelihood.
int RANDBG = 0;			// randomly pick a binding site for low ratio genes. 
bool VERBOSE = false;	// verbose mode. output execution procedure. 

OrfMotif::OrfMotif()
{}

OrfMotif::OrfMotif(string name, VLoc locOrg, VLoc locRev, string seq, double val, double prob, 
				   int label, double avg_prob, bool target, unsigned loc, char orien)
{
	this->name = name;
	this->locOrg = locOrg;
	this->locRev = locRev;
	this->seq = seq;
	this->val = val;
	this->prob = prob;
	this->label = label;
	this->avg_prob = avg_prob;
	this->target = target;
	this->loc = loc;
	this->orien = orien;
}

OrfMotif::OrfMotif(const OrfMotif& r)
:name(r.name), locOrg(r.locOrg), locRev(r.locRev), seq(r.seq), val(r.val), prob(r.prob), label(r.label), 
avg_prob(r.avg_prob), target(r.target), loc(r.loc), orien(r.orien)
{}

OrfMotif::~OrfMotif()
{}	


int main(int argc, char* argv[])
{
        // Read in parameters from command line.
        CCmdLine cmdLine;

        if(cmdLine.SplitLine(argc, argv) < 3)
        {
                cerr << "Usage: ./gittar -m consensus_motif -e log_ratio_file -s promoter_sequence" << endl;
		cerr << "Additional parameters:" << endl;
		cerr << "[-mis] mismatches to consensus motif(Default=1)" << endl;
		cerr << "[-a] annotation file" << endl;
		cerr << "[-o] output(gene_and_other_info PSFM)" << endl;
		cerr << "[-v] verbose mode" << endl;
		cerr << "[-h] hyperprio parameters(hyper_samples=50 mu=2.0 sigma^2=1.0 flanking_length=7" << endl;
		cerr << "[-n] DO NOT use prio information(overide hyperprio parameters)" << endl;
		cerr << "[-g] Gibbs sampler parameters(starting_point=10 burn-in=100 window=100)" << endl;
		cerr << "[-misc] Miscellaneous parameters(stds_of_intial_target=4 skip_core_region=0 random_background=0 threshold_for_target=0.5 stds_enfore=1)" << endl;

                return 1;
        }

        string m, e, s;
        try
        {
                m = cmdLine.GetArgument("-m", 0);
                e = cmdLine.GetArgument("-e", 0);
                s = cmdLine.GetArgument("-s", 0);
        }
        catch(int)
        {
                cerr << "Wrong arguments!" << endl;
                return 1;
        }
	string mis = cmdLine.GetSafeArgument("-mis", 0, "1");	// mismatches.
	unsigned misallow = atoi(mis.data());
	string a = cmdLine.GetSafeArgument("-a", 0, "");	// gene annotation file.
	string o1 = cmdLine.GetSafeArgument("-o", 0, "");	// output gene and other information.
	string o2 = cmdLine.GetSafeArgument("-o", 1, "");	// output PSFM.
	
	// Hyper-prio parameters.
	HypePrio hypePrio;
	string v = cmdLine.GetSafeArgument("-h", 0, "50");
	hypePrio.V = atoi(v.data());
	string mu = cmdLine.GetSafeArgument("-h", 1, "2.0");
	hypePrio.MU = atof(mu.data());
	string sigma2 = cmdLine.GetSafeArgument("-h", 2, "1.0");
	hypePrio.SIGMA_2 = atof(sigma2.data());
	string flanking = cmdLine.GetSafeArgument("-h", 3, "7");
	hypePrio.FLANKING = atoi(flanking.data());

	// DO NOT use prior information? 
	if(cmdLine.HasSwitch("-n"))
		hypePrio.V = 1;

	// Gibbs-sampler parameters.
	string sta_pnt = cmdLine.GetSafeArgument("-g", 0, "10");
	STA_PNT = atoi(sta_pnt.data());
	string burn_in = cmdLine.GetSafeArgument("-g", 1, "100");
	BURN_IN = atoi(burn_in.data());
	string window = cmdLine.GetSafeArgument("-g", 2, "100");
	WINDOW = atoi(window.data());
 
	// Miscellaneous parameters.
	string iniv = cmdLine.GetSafeArgument("-misc", 0, "4");
	INIV = atoi(iniv.data());
	string skipcore	= cmdLine.GetSafeArgument("-misc", 1, "0");
	SKIPCORE = atoi(skipcore.data());
	string randbg = cmdLine.GetSafeArgument("-misc", 2, "0");
	RANDBG = atoi(randbg.data());
	VERBOSE = cmdLine.HasSwitch("-v");	// verbose mode switch.


	// Load data into memory.
	string motifSeq; // string of motif sequence.
	Expr expr; // hash table of expression values.
	Seq seq; // hash table of sequences.
	GeneAnno geneAnno; // hash table of gene annotation info.

	if(VERBOSE)
		cout << "Loading...motif...expression...sequence..." << endl;
	if(load_data(m, e, s, a, motifSeq, expr, seq, geneAnno)) // load data.
	{
		cerr << "Error occurs during data loading!" << endl;
		return 1;
	}
	if(VERBOSE)
		cout << "Data load complete!" << endl;

	// Normalize expression levels.
	double bkm = mean(expr);
	double bkstd = stnd(expr, bkm);
	if(VERBOSE)
		cout << "All gene mean: " << bkm << ", standard deviation: " << bkstd << endl;

	// Find all genes containing the motif with mismatch.
	VGene vGene;
	motifgene(motifSeq, hypePrio.FLANKING, expr, seq, misallow, vGene);
	if(VERBOSE)
		cout << "Total genes identified: " << vGene.size() << endl;
		
	// Calculate pseudo-count for DNA bases among all sequences.
	double dnaPor[4];
	sudocnt(seq, dnaPor);
	if(VERBOSE)
	{
		cout << "DNA bases portion in all sequences." << endl;
		cout << "A: " << dnaPor[0] << endl;
		cout << "C: " << dnaPor[1] << endl;
		cout << "G: " << dnaPor[2] << endl;
		cout << "T: " << dnaPor[3] << endl;
	}
				
	// Gibbs sampler starts from here...
	double tPSFM[4][40], nPSFM[4][40]; // target and non-target PSFM.
	double mTarget, s2Target;
	SST eLen = motifSeq.length() + 2*hypePrio.FLANKING;
	setzero(tPSFM, 4, (unsigned)eLen);
	setzero(nPSFM, 4, (unsigned)eLen);
	mTarget = s2Target = 0.0;
	for(unsigned i = 0; i < STA_PNT; i++)
	{
		// Generate candidate genes' initial labels.
		genlbl(vGene, bkm, bkstd);

		// Use Gibbs Sampler to calculate the estimated labels.
		if(VERBOSE)
			cout << "Gibbs sampler starting..." << endl;
		double mTarget_, s2Target_;
		double tPSFM_[4][40], nPSFM_[4][40];
		/*	
			Here go the parameters information: 
			vGene: Information about all genes that contain motif.
			dnaPor: DNA pseudo-count.
			seq: Promoter sequences for all genes.
			mVal: background expression mean.
			sVal: background expression std.
		*/
		int iter;
		iter = gsamp(vGene, mTarget_, s2Target_, tPSFM_, nPSFM_, misallow, dnaPor, motifSeq, seq, bkm, bkstd, hypePrio);
		if(VERBOSE)
			cout << "Gibbs sampler stops with " << iter << " iterations." << endl;
		mTarget += mTarget_;
		s2Target += s2Target_;
		add2mat(tPSFM, tPSFM, tPSFM_, (unsigned)eLen);
		add2mat(nPSFM, nPSFM, nPSFM_, (unsigned)eLen);
		if(VERBOSE)
		{
			cout << "Target genes in last iteration: " << ntar(vGene) << endl;
			cout << "Non-target genes in last iteration: " << nnon(vGene) << endl;
		}
		sum_iter(i, vGene);
	}
	// Summarize the results from different starting points.
	avg_iter(vGene);
	unsigned nTar = ntar(vGene, true);
	unsigned nNon = nnon(vGene, true);
	mTarget /= STA_PNT;
	s2Target /= STA_PNT;
	normat(tPSFM, (unsigned)eLen, STA_PNT);
	normat(nPSFM, (unsigned)eLen, STA_PNT);
		
	// Write results into files.
	if(VERBOSE)
		cout << "Writing results into files..." << endl;
	time_t time_tag = time(NULL);
	if(o2 == "")
	{
		ostringstream logstream;
		logstream << time_tag << "_" << motifSeq << "_mis" << misallow << ".log";
		o2 = logstream.str();
	}

	// PSFMs and profile ratios.	
	if(write_log(o2, (unsigned)eLen, tPSFM, nPSFM) != 0)
		cerr << "Writing matrices information failed!" << endl;
	else if(VERBOSE)
		cout << "PSFMs and profile ratios write complete!" << endl;

	// Re-calculate binding site according to the averaged PSFM.
	map<string, string> mbindloc;	// A map to store all binding sites of each gene.
	for(unsigned i = 0; i < vGene.size(); i++)
		mbindloc[vGene[i].name] = bindsite(vGene[i], vGene[i].target, seq[vGene[i].name], tPSFM, nPSFM, hypePrio.FLANKING);

	// Write genes information.
	// Sort genes according to their probabilities and expression values.
	sort(vGene.begin(), vGene.end(), cmp);
	if(o1 == "")
	{
		ostringstream infostream;
		infostream << time_tag << "_" << motifSeq << "_mis" << misallow << ".info";
		o1 = infostream.str();
	}
	if(write_info(o1, vGene, geneAnno, mTarget, sqrt(s2Target), bkm, bkstd, 
		(unsigned)seq.size(), nTar, nNon, hypePrio, motifSeq, mbindloc) != 0)
		cerr << "Writing genes information failed!" << endl;
	else if(VERBOSE)
	{
		cout << "Genes information write complete!" << endl;
		cout << "Finished!" << endl;
	}
	
	return 0;

}




