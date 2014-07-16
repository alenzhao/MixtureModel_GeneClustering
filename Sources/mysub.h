/*	GITTAR: mysub.h

	Header file for "mysub.cpp", "mathsub.cpp" and "Gibbs.cpp".
*/


#include <cctype>
#include <string>
#include <map>
#include <vector>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <limits>
#include <sstream>

using namespace std;

extern unsigned BURN_IN;	// iterations required for burn-in.
extern unsigned WINDOW;		// number of iterations used to calculate probabilities.
extern unsigned STA_PNT;	// number of random starting points.
extern unsigned ENFORN;		// number of standard deviations that should be enforced between target mean and background mean.
extern unsigned INIV;		// Initial labels based a number of standard deviations.
extern double THRLD;		// threshold for target probability.
extern int SKIPCORE;		// skip core motif region when calculate sequence likelihood.
extern int RANDBG;		// randomly pick a binding site for low ratio genes. 
extern bool VERBOSE;		// switch for verbose mode. Default is off.

typedef string::size_type SST;
typedef map<char, double> DNAPor; // hash table for nucleotide base to pseudo-count.
typedef map<string, double> Expr; // hash table for expression ratios.
typedef map<string, string> Seq; // hash table for sequences.
typedef struct {
	unsigned loc;
	unsigned mis;
} MatchSite;
typedef vector<MatchSite> VLoc;	// vector to store binding locations.
typedef struct {
	string geneName;
	string geneDesc;
	string refSeq;
	string locusLink;
} GeneInfo;	// struct to store information about each gene.
typedef map<string, GeneInfo> GeneAnno; // hash table for gene information/annotation.

typedef struct {
	double n, s_2, avg, SIGMA_2, MU, V;
} ExprPost;

typedef struct {
	unsigned FLANKING; // length of flanking region.
	double SIGMA_2, MU; // shape, inverse-scale, mean.
	unsigned V; // number of hyperthetical samples.
} HypePrio;	// struct for hype-prior info.

class OrfMotif // class for each candidate gene's information.
{
public:
	OrfMotif();
	OrfMotif(string name, VLoc locOrg, VLoc locRev, string seq, double val, 
		double prob, int label, double avg_prob, bool target, unsigned loc, char orien);
	OrfMotif(const OrfMotif& r);
	~OrfMotif();
	string name;
	VLoc locOrg, locRev;
	double val;
	// These parameters are for each step in each iteration.
	string seq;
	double prob;
	int label;
	// These parameters are for the final results.
	double avg_prob;
	bool target;
	unsigned loc;
	char orien;
};
typedef vector<OrfMotif> VGene; // vector for ORF that contain motif.

typedef struct {
	unsigned TP;
	unsigned FP;
	unsigned TN;
	unsigned FN;
} Chars;	// struct for statistics of classifier performance: true positive, false negative, etc.

string& str2upper(string& str);
string& str2lower(string& str);
VLoc locmotif(const string motif, const unsigned FLANKING, const string seq, const unsigned int misallow);
string emotif(const int startLoc, const SST mLen, const unsigned FLANKING, const string& seq);
string copyrev(const string seq);
int gsamp(VGene& vGene, double& mTarget, double& s2Target, double tPSFM[][40], double nPSFM[][40], 
		  const unsigned mismatch, const double dnaPor[4], const string motif, Seq& seq, 
		  const double bkm, const double bkstd, const HypePrio hypePrio, const bool sim = false);
void sumup_seq(VGene& vGene, double tPSFM_sum[][40], double nPSFM_sum[][40], const double dnaPor[4]);
void cal_psfm(VGene& vGene, double tPSFM[][40], double nPSFM[][40], const double dnaPor[4], unsigned& nTar, unsigned& nNon);
void sub_gene(OrfMotif& orfMotif, unsigned& nTar, unsigned& nNon, double tPSFM_sum[][40], 
			  double nPSFM_sum[][40], double& tValSum, double& tValSumSq, const double dnaPor[4]);
void add_gene(OrfMotif& orfMotif, unsigned& nTar, unsigned& nNon, double tPSFM_sum[][40], 
			  double nPSFM_sum[][40], double& tValSum, double& tValSumSq, const double dnaPor[4]);
double mean(const VGene& v);
double mean(bool target, const VGene& v);
double mean(const Expr& expr);
double stnd(const VGene& v, const double m);
double stnd(bool target, const VGene& v, const double m);
double stnd(const Expr& expr, const double m);
double sum_val(bool target, const VGene& v);
double sumsq_val(bool target, const VGene& v);
inline double normpdf(const double x, const double u, const double s);
void setzero(double matrix[][40], const unsigned row, const unsigned col);
void chksign(VGene &v);
char basepair(char base);
bool cmp(OrfMotif a, OrfMotif b);
double seqlk(OrfMotif& gene, double& pTargetSeq, double& pNontSeq, const string seq, 
				  const double tPSFM[][40], const double nPSFM[][40], 
				  const unsigned FLANKING, const string& motif, bool hiratio, bool loratio);
double seqlk(OrfMotif& gene, double psfm[][40], const unsigned FLANKING);
string bindsite(OrfMotif& gene, const bool target, const string seq, const double tPSFM[][40], const double nPSFM[][40], const unsigned FLANKING);
unsigned dnaidx(char base);
void write_psfm(ofstream& h, const unsigned len, const double psfm[][40], const double info[40] = 0);
void write_gene(ofstream& h, const VGene& vGene, const bool target, GeneAnno& geneAnno, map<string, string>& mbindloc);
void normat(double matrix[][40], const unsigned len, const double sum);
void normat(double tarmat[][40], const double orgmat[][40], const unsigned len, const double sum);
void entropy(const double matrix[][40], const unsigned len, double info[40]);
void prof_ratio(const unsigned eLen, const double tMat[][40], const double nMat[][40], double rMat[][40]);
int load_data(const string& fMotif, const string& fExpr, const string& fSeq, const string& fGeneAnno, string& motifSeq, Expr& expr, Seq& seq, GeneAnno& geneAnno);
void motifgene(const string& motifSeq, const unsigned FLANKING, Expr& expr, const Seq& seq, const unsigned mis, VGene& vGene);
void motifgene(Expr& expr, const Seq& seq, VGene& vGene);
void genlbl(VGene& vGene, const double bkm, const double bkstd);
void sudocnt(const Seq& seq, double portion[4]);
int write_log(const string logname, const unsigned len, double tPSFM[][40], double nPSFM[][40]);
void write_log(ofstream& log_h, const unsigned len, double tPSFM[][40], double nPSFM[][40]);
int write_info(const string infoname, const VGene& vGene, GeneAnno& geneAnno, const double tm, const double tstd, 
			   const double bkm, const double bks, const unsigned nTotal, const unsigned nTar, const unsigned nNon, 
			   const HypePrio& hypePrio, const string motif, map<string, string>& mbindloc);
unsigned ntar(const VGene& v, bool final = false);
unsigned nnon(const VGene& v, bool final = false);
double calc_lambda(const VGene& v);
inline void post_expr(const unsigned v, const double mu, const double sigma_2, 
					  const unsigned n, const double avg, const double s_2, 
					  double& shift, double& scale, unsigned& dof);
inline double prob_expr(const double x, const double alpha, const double beta, const double mu);
void prior_seq(const unsigned mismatch, const double dnaPor[4], const string motif, 
			   const unsigned eLen, const unsigned FLANKING, const unsigned V, double U[][40]);
void prior_seq(const double dnaPor[4], const string motif, const unsigned eLen, const unsigned FLANKING, 
			   const unsigned V, double U[][40]);
unsigned chose(const unsigned m, const unsigned n);
unsigned factorial(const unsigned n);
void add2mat(double tar_mat[][40], const double mat1[][40], const double mat2[][40], const unsigned len);
void avg_prob(VGene& vGene, const vector<double>& vProb);
void cal_chars(VGene& vGene, Chars& charsInfo);
inline double ran_freq(const double psfm[][40], unsigned col);
unsigned numofN(const string motif);
int Bernoulli(const double theta);
int Multinomial(const double log_p[255], const unsigned n);
inline void safeguard(double& mTarget, double& s2Target, const double bkm, const double bkstd, const double val);
void sum_iter(const unsigned i, VGene& vGene);
void avg_iter(VGene& vGene);
bool iupac(const char nuc, const char pac);
double matscore(double& pt, double& pn, const string& seq, const double tPSFM[][40], const double nPSFM[][40], unsigned FLANKING, 
				const string& motif, bool hiratio);


/********* Stat routines from others ***************/
double gamma ( double x );
double student_pdf ( double x, double a, double b, double c );
double d_huge ( void );
double d_epsilon ( void );
double normal_pdf ( double x, double a, double b );
/***************************************************/
















