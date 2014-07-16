#include "mysub.h"
// All math subroutines.
// Run Gibbs sampling procedure. 
int gsamp(VGene& vGene, double& mTarget, double& s2Target, double tPSFM[][40], double nPSFM[][40], 
		  const unsigned mismatch, const double dnaPor[4], const string motif, Seq& seq, 
		  const double bkm, const double bkstd, const HypePrio hypePrio, const bool sim)
{
	SST eLen = vGene[0].seq.length(); // Extended motif length.
	// Calculation of parameters for all genes.
	// In each iteration, parameters are updated based on these values.
	double tPSFM_sum[4][40], nPSFM_sum[4][40];	// Sum of weights of each base on each position.
	setzero(tPSFM_sum, 4, (unsigned)eLen); 
	setzero(nPSFM_sum, 4, (unsigned)eLen); // Set all elements of PSFM to zero.
	sumup_seq(vGene, tPSFM_sum, nPSFM_sum, dnaPor);
	
	unsigned nTar = ntar(vGene);
	unsigned nNon = nnon(vGene);
	double lambda0; // lambda value for the whole population.
	double lambda;	// lambda value of the current run.

	double tValSum = sum_val(true, vGene); // sum for targets.
	double tValSumSq = sumsq_val(true, vGene); // squared sum for targets.
	double nValSum = sum_val(false, vGene); // sum for non-targets.
	double nValSumSq = sumsq_val(false, vGene); // squared sum for non-targets.
	mTarget = 0.0;
	s2Target = 0.0; // mean and variance of target, averaged through iterations after convergence.
	setzero(tPSFM, 4, (unsigned)eLen);
	setzero(nPSFM, 4, (unsigned)eLen);

	double U0[4][40], U1[4][40];
	// Calculate values for hyper-parameters of prior Dirichlet distribution.
	prior_seq(mismatch, dnaPor, motif, (unsigned)eLen, hypePrio.FLANKING, hypePrio.V, U0);
	prior_seq(dnaPor, motif, (unsigned)eLen, hypePrio.FLANKING, hypePrio.V, U1);

	vector<double> vProb; // vector to store sum of weights after burn-in. 
	vProb.assign(vGene.size(), 0.0);

	unsigned iter = 0;
	srand((unsigned)time(NULL));	// random seed for Bernoulli.
	while(true) // judge convergence.
	{
		iter++; // Number of iteration.
		if(nTar == 0) // Check whether convergence goes to edge.
		{
			cerr << "Target gene number becomes zero! Process stops." << endl;
			break;
		}
		else if(nNon == 0)
		{
			cerr << "Non-target gene number becomes zero! Process stops." << endl;
			break;
		}

		double loglike = 0.0; // sum of predictive log-likelihoods of all genes.
		unsigned nCal = 0;	// number of genes that are involved in calculation.
		for(unsigned i = 0; i < vGene.size(); i++) // Estimate all gene's label.
		{
						
			sub_gene(vGene[i], nTar, nNon, tPSFM_sum, nPSFM_sum, tValSum, tValSumSq, dnaPor); // subtract gene i.
			double tPSFM_[4][40], nPSFM_[4][40]; // PSFM (in current iteration) for target or non-target.
			// Calculate posterior for sequence parameters.
			add2mat(tPSFM_, tPSFM_sum, U1, (unsigned)eLen);
			add2mat(nPSFM_, nPSFM_sum, U0, (unsigned)eLen);
			double u1_sum = nTar + 1 + hypePrio.V;
			double u0_sum = nNon + 1 + hypePrio.V;
			// Store normalized PSFM into current PSFM.
			normat(tPSFM_, (unsigned)eLen, u1_sum);
			normat(nPSFM_, (unsigned)eLen, u0_sum);

			nCal++;
			// mean and std of expression values for target genes.
			double mTarget_ = tValSum/nTar; // mean of current step.
			double s2Target_ = (tValSumSq + (mTarget_*mTarget_)*nTar - 2*mTarget_*tValSum)/(nTar-1); // variance of current step.
			if(mTarget_ < bkm + ENFORN*sqrt(s2Target_)) // Add enforcement to separate target and background means.
				mTarget_ = bkm + ENFORN*sqrt(s2Target_);
			//safeguard(mTarget_, s2Target_, hypePrio.positive, bkm, bkstd, vGene[i].val);	// safe-guarding mean and variance, avoid local maximum.

			lambda = ((double)nTar+0.5)/(double)vGene.size(); // mixing parameter.
			// Predictive posterior parameters.
			double shift, scale;
			unsigned dof;
			post_expr(hypePrio.V, hypePrio.MU, hypePrio.SIGMA_2, nTar, mTarget_, s2Target_, shift, scale, dof);

			// Calculate posterior probability based on expression data.
			double pTargetExp, pSafe, pNontExp;
			// Target probability.
			if(dof < 30)
			{
				pTargetExp = log(student_pdf(vGene[i].val, shift, scale, dof));
				pSafe = log(student_pdf(vGene[i].val, shift, bkstd, dof)); // safe-guard variance.
			}
			else
			{
				pTargetExp = log(normal_pdf(vGene[i].val, shift, scale));
				pSafe = log(normal_pdf(vGene[i].val, shift, bkstd)); // safe-guard variance.
			}
			// Non-target probability.
			if(nNon - 1 < 30)
				pNontExp = log(student_pdf(vGene[i].val, bkm, bkstd, nNon-1));
			else // if DOF is larger than 30, use normal distribution to approx.
				pNontExp = log(normal_pdf(vGene[i].val, bkm, bkstd));

			// Safe-guard technique from Ron!
			if(vGene[i].val < bkm && pNontExp - pTargetExp < pNontExp - pSafe)
				pTargetExp = pSafe;
			else if(vGene[i].val > shift && pNontExp - pTargetExp > pNontExp - pSafe)
				pTargetExp = pSafe;

			// Calculate likelihood (target - nontar) based on sequence data.
			double pTargetSeq, pNontSeq;
			if(sim)
			{
				pTargetSeq = seqlk(vGene[i], tPSFM_, hypePrio.FLANKING);
				pNontSeq = seqlk(vGene[i], nPSFM_, hypePrio.FLANKING);
			}
			else
				seqlk(vGene[i], pTargetSeq, pNontSeq, seq[vGene[i].name], tPSFM_, nPSFM_, 
					hypePrio.FLANKING, motif, vGene[i].val > mTarget_, vGene[i].val < bkm);

			// Calculate the probability being target or non-target for gene i.
			double pNont, pTarget;
			pTarget = pTargetExp + pTargetSeq + log(lambda);
			pNont = pNontExp + pNontSeq + log(1-lambda);

			vGene[i].prob = 1/(1 + exp(pNont - pTarget)); // Update the probability being target.
			vGene[i].label = Bernoulli(vGene[i].prob);

			loglike += pTarget + log(1+exp(pNont - pTarget)); // sum up all genes' predictive log-likelihoods.
			add_gene(vGene[i], nTar, nNon, tPSFM_sum, nPSFM_sum, tValSum, tValSumSq, dnaPor); // Add gene i back.

			// Add probs into vector after burn-in.
			if(iter > BURN_IN)
				vProb[i] += vGene[i].prob;
		}
		// Sum up statistics of parameters after burn-in.
		if(iter > BURN_IN)
		{
			// Expression mean and variance.
			double mTarget_ = tValSum/nTar;
			double s2Target_ = (tValSumSq + (mTarget_*mTarget_)*nTar - 2*mTarget_*tValSum)/(nTar-1);
			double shift, scale;
			unsigned dof;
			post_expr(hypePrio.V, hypePrio.MU, hypePrio.SIGMA_2, nTar, mTarget_, s2Target_, shift, scale, dof);
			mTarget += shift;
			s2Target += scale*scale;
			// Sequence frequency matrices.
			double tPSFM_[4][40], nPSFM_[4][40];
			setzero(tPSFM_, 4, (unsigned)eLen);
			setzero(nPSFM_, 4, (unsigned)eLen);
			add2mat(tPSFM_, tPSFM_sum, U1, (unsigned)eLen);
			add2mat(nPSFM_, nPSFM_sum, U0, (unsigned)eLen);
			double u1_sum = nTar + 1 + hypePrio.V;
			double u0_sum = nNon + 1 + hypePrio.V;
			normat(tPSFM_, (unsigned)eLen, u1_sum);
			normat(nPSFM_, (unsigned)eLen, u0_sum);
			add2mat(tPSFM, tPSFM, tPSFM_, (unsigned)eLen);
			add2mat(nPSFM, nPSFM, nPSFM_, (unsigned)eLen);
		}

		lambda0 = (double)nTar/(double)vGene.size();
		if(VERBOSE)
		{
			cout << "Iter: ";
			cout.width(3); cout.setf(ios::left);
			cout << iter << " lambda: ";
			cout.width(20); cout.setf(ios::left); cout.setf(ios::showpoint);
			cout << lambda0 << " log-likelihood: ";
			cout.width(20); cout.setf(ios::left); cout.setf(ios::showpoint);
			cout << eLen + 3 + loglike/nCal << endl;
		}

		if(iter >= BURN_IN + WINDOW)
			break;
	}
	
	avg_prob(vGene, vProb);
	mTarget /= WINDOW;
	s2Target /= WINDOW;
	normat(tPSFM, (unsigned)eLen, WINDOW);
	normat(nPSFM, (unsigned)eLen, WINDOW);

//	chksign(vGene); // check whether necessary to inverse signs.
		
	return iter;
}

// Calculate the counts of DNA bases at each position.
void sumup_seq(VGene& vGene, double tPSFM_sum[][40], double nPSFM_sum[][40], const double dnaPor[4])
{
	SST eLen = vGene[0].seq.length(); // Extended motif length.
	for(unsigned i = 0; i < vGene.size(); i++) // Calculate PSFM and probability distribution.
	{
		for(unsigned j = 0; j < eLen; j++) // Statistics of nucleotides.
		{
			// 1. k-th base letter from j-th gene sequence,
			// 2. mapping base letter to numeric index, 
			// 3. updating corresponding PSFM.
			if(vGene[i].seq[j] == 'N' || islower(vGene[i].seq[j]))
			{
				for(unsigned k = 0; k < 4; k++)
				{
					tPSFM_sum[k][j] += vGene[i].label*dnaPor[k];
					nPSFM_sum[k][j] += (1 - vGene[i].label)*dnaPor[k];
				}
			}
			else
			{
				tPSFM_sum[dnaidx(vGene[i].seq[j])][j] += vGene[i].label; // Target PSFM.
				nPSFM_sum[dnaidx(vGene[i].seq[j])][j] += (1 - vGene[i].label); // Non-target PSFM.
			}
		}
	}
	for(unsigned i = 0; i < 4; i++) // Add pseudo-count.
	{
		for(unsigned j = 0; j < eLen; j++)
		{
			tPSFM_sum[i][j] += dnaPor[i];
			nPSFM_sum[i][j] += dnaPor[i];
		}
	}
}

// Subtract gene i from all statistics.
void sub_gene(OrfMotif& orfMotif, unsigned& nTar, unsigned& nNon, double tPSFM_sum[][40], 
			  double nPSFM_sum[][40], double& tValSum, double& tValSumSq, const double dnaPor[4])
{
	if(orfMotif.label == 1) // update number of target and non-targets.
		nTar--;
	else
		nNon--;
	SST eLen = orfMotif.seq.length();
	// Subtract gene i from PSFM.
	for(unsigned j = 0; j < eLen; j++)
	{
		if(orfMotif.seq[j] == 'N' || islower(orfMotif.seq[j]))
		{
			for(unsigned k = 0; k < 4; k++)
			{
				tPSFM_sum[k][j] -= orfMotif.label*dnaPor[k];
				nPSFM_sum[k][j] -= (1 - orfMotif.label)*dnaPor[k];
			}
		}
		else
		{
			tPSFM_sum[dnaidx(orfMotif.seq[j])][j] -= orfMotif.label; // Target PSFM.
			nPSFM_sum[dnaidx(orfMotif.seq[j])][j] -= (1 - orfMotif.label); // Non-target PSFM.
		}
	}

	// Subtract gene i from expression sums.
	tValSum -= orfMotif.val*orfMotif.label;
	tValSumSq -= (orfMotif.val*orfMotif.val)*orfMotif.label;
}

// Add gene i back to all statistics.
void add_gene(OrfMotif& orfMotif, unsigned& nTar, unsigned& nNon, double tPSFM_sum[][40], 
			  double nPSFM_sum[][40], double& tValSum, double& tValSumSq, const double dnaPor[4])
{
	if(orfMotif.label == 1) // update number of target and non-targets.
		nTar++;
	else
		nNon++;
	SST eLen = orfMotif.seq.length();
	// Add gene i back PSFM.
	for(unsigned j = 0; j < eLen; j++)
	{
		if(orfMotif.seq[j] == 'N' || islower(orfMotif.seq[j]))
		{
			for(unsigned k = 0; k < 4; k++)
			{
				tPSFM_sum[k][j] += orfMotif.label*dnaPor[k];
				nPSFM_sum[k][j] += (1 - orfMotif.label)*dnaPor[k];
			}
		}
		else
		{
			tPSFM_sum[dnaidx(orfMotif.seq[j])][j] += orfMotif.label; // Target PSFM.
			nPSFM_sum[dnaidx(orfMotif.seq[j])][j] += (1 - orfMotif.label); // Non-target PSFM.
		}
	}

	// Add gene i back expression sums.
	tValSum += orfMotif.val*orfMotif.label;
	tValSumSq += (orfMotif.val*orfMotif.val)*orfMotif.label;
}

int Bernoulli(const double theta) // Sample from Bernoulli distribution.
{
	if((double)rand()/RAND_MAX > theta)
		return 0;
	else
		return 1;
}

int Multinomial(const double log_p[255], const unsigned n) // Sample from Multinomial distribution.
{
//	double p;

	return 0;
}

inline void safeguard(double& mTarget, double& s2Target, const double bkm, const double bkstd, const double val)
{
	if(mTarget < bkm + ENFORN*bkstd)	// Force means to separate.
		mTarget = bkm + ENFORN*bkstd;

	// Avoid large variance to mislead.
	if((val < bkm && s2Target > bkstd*bkstd) || (val > mTarget && s2Target < bkstd*bkstd))
		s2Target = bkstd*bkstd;
}

double mean(bool target, const VGene& v) // Calculate target or non-target mean.
{
	double sum = 0.0;
	unsigned num = 0;
	for(unsigned i = 0; i < v.size(); i++)
	{
		if(target && v[i].target)
		{
			sum += v[i].val;
			num++;
		}
		else if(!target && !v[i].target)
		{
			sum += v[i].val;
			num++;
		}
	}

	return sum/(double)num;
}

double mean(const Expr& expr) // Calculate mean for all genes.
{
	double sum = 0.0;
	for(Expr::const_iterator ci = expr.begin(); ci != expr.end(); ci++)
		sum += ci->second;
	
	return sum/(double)expr.size();
}

double mean(const VGene& vGene) // Calculate mean for all candidate genes.
{
	double sum = 0.0;
	for(unsigned i = 0; i < vGene.size(); i++)
		sum += vGene[i].val;
	
	return sum/(double)vGene.size();
}
	
// Calculate standard deviation.
double stnd(bool target, const VGene& v, const double m)
{
	double sq = 0.0;
	unsigned num = 0;
	for(unsigned int j = 0; j < v.size(); j++)
	{
		if(target && v[j].target)
		{
			sq += pow(v[j].val - m, 2);
			num++;
		}
		else if(!target && !v[j].target)
		{
			sq += pow(v[j].val - m, 2);
			num++;
		}
	}		
	sq /= (double)(num-1);
	sq = sqrt(sq);

	return sq;
}

double stnd(const Expr& expr, const double m) // Calculate standard deviation for all genes.
{
	double sq = 0.0;
	for(Expr::const_iterator ci = expr.begin(); ci != expr.end(); ci++)
		sq += pow((ci->second - m),2);
	sq /= (double)(expr.size() - 1);
	sq = sqrt(sq);

	return sq;
}

double stnd(const VGene& vGene, const double m) // Calculate standard deviation for all candidiate genes.
{
	double sq = 0.0;
	for(unsigned i = 0; i < vGene.size(); i++)
		sq += pow(vGene[i].val - m, 2);
	sq /= (double)(vGene.size() - 1);
	sq = sqrt(sq);

	return sq;
}

// Calculate normal (1-side) probabilistic density function.
inline double normpdf(const double x, const double u, const double s) 
{
	if(x == u)
		return -log(s);
	else
		return -pow((x-u),2)/(2*pow(s,2)) - log(s);
}

void setzero(double matrix[][40], const unsigned row, const unsigned col) // Set all elements of matrix to zero.
{
	for(unsigned i = 0; i < row; i++)
		for(unsigned j = 0; j < col; j++)
			matrix[i][j] = 0.0;
}

double matscore(double& pt, double& pn, const string& seq, const double tPSFM[][40], const double nPSFM[][40], 
				unsigned FLANKING, const string& motif, bool hiratio)
{
	pt = pn = 0.0;
	for(unsigned i = 0; i < seq.length(); i++)
	{
		if(SKIPCORE != 0 && i >= FLANKING && i < FLANKING + motif.length()) // skip core motif region if specified.
			continue;
		if(seq[i] == 'N' || islower(seq[i]))
		{
			pt += ran_freq(tPSFM, i);
			pn += ran_freq(nPSFM, i);
		}
		else
		{
			pt += log(tPSFM[dnaidx(seq[i])][i]);
			pn += log(nPSFM[dnaidx(seq[i])][i]);
		}
	}

	return pt - pn;
}

// Calculate the maximum likelihood for a gene being target given sequence and PSFM.
double seqlk(OrfMotif& gene, double& pTargetSeq, double& pNontSeq, const string seq, 
				  const double tPSFM[][40], const double nPSFM[][40], 
				  const unsigned FLANKING, const string& motif, bool hiratio, bool loratio)
{
	SST mLen = gene.seq.length() - 2*FLANKING;
	double pSeq0;
	if(!RANDBG || !loratio)
	{
		unsigned mis0;
		bool tag = true;	// true if there is no binding site at original sequence.
		// Consider original sequence first.
		for(unsigned i = 0; i < gene.locOrg.size(); i++)
		{
			string emSeq = emotif(gene.locOrg[i].loc, mLen, FLANKING, seq);
			double pTargetSeq_, pNontSeq_;
			double pSeq = matscore(pTargetSeq_, pNontSeq_, emSeq, tPSFM, nPSFM, FLANKING, motif, hiratio);
			// Find the sequence segment with maximum likelihood.
			//if(i==0 || gene.locOrg.at(i).mis < mis0 || pSeq > pSeq0 && gene.locOrg.at(i).mis <= mis0)
			if(i==0 || pSeq > pSeq0)
			{
				tag = false;
				pTargetSeq = pTargetSeq_;
				pNontSeq = pNontSeq_;
				pSeq0 = pSeq;
				mis0 = gene.locOrg.at(i).mis;
				gene.seq = emSeq;
				gene.loc = gene.locOrg.at(i).loc;
				gene.orien = 'F';
			}
		}
		// Consider reverse sequence then.
		string revseq = copyrev(seq);
		for(unsigned i = 0; i < gene.locRev.size(); i++)
		{
			string emSeq = emotif(gene.locRev[i].loc, mLen, FLANKING, revseq);
			double pTargetSeq_, pNontSeq_;
			double pSeq = matscore(pTargetSeq_, pNontSeq_, emSeq, tPSFM, nPSFM, FLANKING, motif, hiratio);
			// Find the sequence segment with maximum likelihood.
			//if(tag || gene.locRev.at(i).mis < mis0 || pSeq > pSeq0 && gene.locRev.at(i).mis <= mis0)
			if(tag || pSeq > pSeq0)
			{
				pTargetSeq = pTargetSeq_;
				pNontSeq = pNontSeq_;
				pSeq0 = pSeq;
				mis0 = gene.locRev.at(i).mis;
				gene.seq = emSeq;
				gene.loc = gene.locRev.at(i).loc;
				gene.orien = 'R';
			}
		}
	}
	else // If genes have low ratios, randomly pick a binding site.
	{
		unsigned nBind = (unsigned)gene.locOrg.size() + (unsigned)gene.locRev.size();
		unsigned idx = (unsigned)((double)(rand()-1)/RAND_MAX*nBind);
		if(idx < gene.locOrg.size()) // Select binding site from forward strand.
		{
			gene.seq = emotif(gene.locOrg[idx].loc, mLen, FLANKING, seq);
			gene.loc = gene.locOrg.at(idx).loc;
			gene.orien = 'F';
		}
		else // Select binding site from reverse strand.
		{
			string revseq = copyrev(seq);
			idx -= (unsigned)gene.locOrg.size();
			gene.seq = emotif(gene.locRev[idx].loc, mLen, FLANKING, revseq);
			gene.loc = gene.locRev.at(idx).loc;
			gene.orien = 'R';
		}
		pSeq0 = matscore(pTargetSeq, pNontSeq, gene.seq, tPSFM, nPSFM, 0, "", false);
	}
	
	return pSeq0;
}

// Determine one gene's binding site based on ratio between target and non-target PSFMs.
string bindsite(OrfMotif& gene, const bool target, const string seq, 
			  const double tPSFM[][40], const double nPSFM[][40], const unsigned FLANKING)
{
	SST mLen = gene.seq.length() - 2*FLANKING;
	ostringstream strmBind;
	if(target)	// If target, select the segment of sequence yielding the maximum likelihood.
	{
		double pSeq0;
		bool tag = true;	// true if there is no binding site at original sequence.
		// Consider original sequence first.
		for(unsigned i = 0; i < gene.locOrg.size(); i++)
		{
			unsigned loc = gene.locOrg[i].loc;
			string emSeq = emotif(loc, mLen, FLANKING, seq);
			double pTar, pNon;
			double pSeq = matscore(pTar, pNon, emSeq, tPSFM, nPSFM, 0, "", false);
			strmBind << " (" << pSeq << "," << loc << ",F)";
			// Find the sequence segment with maximum likelihood (target) OR minimum likelihood (non-target).
			if(i==0 || pSeq > pSeq0)
			{
				tag = false;
				pSeq0 = pSeq;
				gene.seq = emSeq;
				gene.loc = gene.locOrg.at(i).loc;
				gene.orien = 'F';
			}
		}
		// Consider reverse sequence then.
		string revseq = copyrev(seq);
		for(unsigned i = 0; i < gene.locRev.size(); i++)
		{
			unsigned loc = gene.locRev[i].loc;
			string emSeq = emotif(loc, mLen, FLANKING, revseq);
			double pTar, pNon;
			double pSeq = matscore(pTar, pNon, emSeq, tPSFM, nPSFM, 0, "", false);
			strmBind << " (" << pSeq << "," << loc << ",R)";
			if(tag || pSeq > pSeq0)
			{
				pSeq0 = pSeq;
				gene.seq = emSeq;
				gene.loc = gene.locRev.at(i).loc;
				gene.orien = 'R';
			}
		}
	}
	else	// If non-target, randomly pick one segment of sequence.
	{
		unsigned nBind = (unsigned)gene.locOrg.size() + (unsigned)gene.locRev.size();
		unsigned idx = (unsigned)((double)(rand()-1)/RAND_MAX*nBind);
		if(idx < gene.locOrg.size()) // Select binding site from forward strand.
		{
			gene.seq = emotif(gene.locOrg[idx].loc, mLen, FLANKING, seq);
			gene.loc = gene.locOrg.at(idx).loc;
			gene.orien = 'F';
		}
		else // Select binding site from reverse strand.
		{
			string revseq = copyrev(seq);
			idx -= (unsigned)gene.locOrg.size();
			gene.seq = emotif(gene.locRev[idx].loc, mLen, FLANKING, revseq);
			gene.loc = gene.locRev.at(idx).loc;
			gene.orien = 'R';
		}
	}

	return strmBind.str();
}

// Calculate the maximum likelihood for a gene being non-target given sequence and PSFM.
double seqlk(OrfMotif& gene, double psfm[][40], const unsigned FLANKING)
{
	double pSeq = 0.0;
	for(unsigned i = 0; i < gene.seq.length(); i++)
	{
		if(gene.seq[i] == 'N')
			pSeq += ran_freq(psfm, i);
		else
			pSeq += log(psfm[dnaidx(gene.seq[i])][i]);
	}
	
	return pSeq;
}

unsigned dnaidx(char base) // Give an index to DNA base.
{
	switch(base)
	{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		default:
			return 99;
	}
}

// Normalize DNA frequency matrix and store into the original matrix.
void normat(double matrix[][40], const unsigned len, const double sum)
{
	for(unsigned i = 0; i < 4; i++)
	{
		for(unsigned j = 0; j < len; j++)
			matrix[i][j] /= sum;
	}
}

// Normalize DNA frequency matrix and store into another matrix.
void normat(double tarmat[][40], const double orgmat[][40], const unsigned len, const double sum)
{
	for(unsigned i = 0; i < 4; i++)
	{
		for(unsigned j = 0; j < len; j++)
			tarmat[i][j] = orgmat[i][j]/sum;
	}
}

// Number of target genes.
unsigned ntar(const VGene& v, bool final)
{
	unsigned num = 0;
	for(unsigned j = 0; j < v.size(); j++)
	{
		if(final)
		{
			if(v[j].target)
				num++;
		}
		else
		{
			if(v[j].label == 1)
				num++;
		}
	}
	
	return num;
}

// Number of non-target genes.
unsigned nnon(const VGene& v, bool final)
{
	unsigned num = 0;
	for(unsigned j = 0; j < v.size(); j++)
	{
		if(final)
		{
			if(!v[j].target)
				num++;
		}
		else
		{
			if(v[j].label == 0)
				num++;
		}
	}
	
	return num;
}

/*double calc_lambda(const VGene& v) // Calculate lambda value.
{
	double sum = 0.0;
	for(unsigned i = 0; i < v.size(); i++)
		sum += v[i].weight;

	return sum/(double)v.size();
}*/

double sum_val(bool target, const VGene& v) // sum of gene expression values.
{
	double sum = 0.0;
	for(unsigned i = 0; i < v.size(); i++)
	{
		if(target)
			sum += v[i].val*v[i].label;
		else
			sum += v[i].val*(1-v[i].label);
	}

	return sum;
}

double sumsq_val(bool target, const VGene& v) // squared sum of gene expression values.
{
	double sum = 0.0;
	for(unsigned i = 0; i < v.size(); i++)
	{
		if(target)
			sum += (v[i].val*v[i].val)*v[i].label;
		else
			sum += (v[i].val*v[i].val)*(1-v[i].label);
	}

	return sum;
}

// Calculate predictive posterior parameter values for expression data.
inline void post_expr(const unsigned v, const double mu, const double sigma_2, 
					  const unsigned n, const double avg, const double s_2, 
					  double& shift, double& scale, unsigned& dof)
{
	double vn = v + n;
	shift = (v*mu + n*avg)/vn;
	double varn;
	varn = v*sigma_2 + (n-1)*s_2 + v*n*(avg-mu)*(avg-mu)/vn;
	varn /= vn;
	scale = sqrt((1+1/vn)*varn);
	dof = (unsigned)(vn - 1);
}

// Calculate posterior probability based on expression data.
inline double prob_expr(const double x, const double alpha, const double beta, const double mu)
{
	return log((alpha-1)/beta)/2 - (x - mu)*(x - mu)/2*(alpha-1)/beta;
}

// Set the values of hyper-parameters for background sequence. 
void prior_seq(const unsigned mismatch, const double dnaPor[4], const string motif, 
			   const unsigned eLen, const unsigned FLANKING, const unsigned V, double U[][40])
{
	unsigned nWrong = 0; // number of combinations if one position is mismatched.
	unsigned nTotal = 0; // number of all combinations with mismatches.
	if(mismatch > 0)
	{
		for(unsigned i = 0; i <= mismatch - 1; i++)
			nWrong += chose(i, (unsigned)motif.length()-numofN(motif)-1)*(unsigned)pow((double)3, (int)i);

		for(unsigned i = 0; i <= mismatch; i++)
			nTotal += chose(i, (unsigned)motif.length()-numofN(motif))*(unsigned)pow((double)3, (int)i);
	}

	for(unsigned i = 0; i < eLen; i++)
	{
		if(i >= FLANKING && i < FLANKING+motif.length()) // core motif region.
		{
			if(motif[i-FLANKING] != 'A' && motif[i-FLANKING] != 'C' && motif[i-FLANKING] != 'G' && motif[i-FLANKING] != 'T')
			{
				for(unsigned j = 0; j < 4; j++)
					U[j][i] = dnaPor[j]*V;
			}
			else
			{
				if(mismatch == 0)
				{
					U[0][i] = U[1][i] = U[2][i] = U[3][i] = 0;
					U[dnaidx(motif[i-FLANKING])][i] = V;
				}
				else // Calculate all possible combinations.
				{
					U[0][i] = U[1][i] = U[2][i] = U[3][i] = ((double)nWrong/nTotal)*V;
					U[dnaidx(motif[i-FLANKING])][i] = ((double)(nTotal - 3*nWrong)/nTotal)*V;
				}
			}
		}
		else // use pseudo-count for flanking regions and base 'N'.
		{
			for(unsigned j = 0; j < 4; j++)
				U[j][i] = dnaPor[j]*V;
		}
	}
}

// Set the values of hyper-parameters for target sequence. 
void prior_seq(const double dnaPor[4], const string motif, const unsigned eLen, const unsigned FLANKING, 
			   const unsigned V, double U[][40])
{
	for(unsigned i = 0; i < eLen; i++)
	{
		if(i >= FLANKING && i < FLANKING+motif.length()) // core motif region.
		{
			if(motif[i-FLANKING] != 'A' && motif[i-FLANKING] != 'C' && 
				motif[i-FLANKING] != 'G' && motif[i-FLANKING] != 'T')
			{
				for(unsigned j = 0; j < 4; j++)
					U[j][i] = dnaPor[j]*V;
			}
			else
			{
				U[0][i] = U[1][i] = U[2][i] = U[3][i] = 0;
				U[dnaidx(motif[i-FLANKING])][i] = V;
			}
		}
		else // use pseudo-count for flanking regions and base 'N'.
		{
			for(unsigned j = 0; j < 4; j++)
				U[j][i] = dnaPor[j]*V;
		}
	}
}

// Chose m items from n items.
unsigned chose(const unsigned m, const unsigned n)
{
	if(m == 0)
		return 1;
	
	return factorial(n)/(factorial(m)*factorial(n-m));
}

// Calculate the factorial for integer.
unsigned factorial(const unsigned n)
{
	if(n == 0 || n == 1)
		return 1;

	unsigned f = 1;
	for(unsigned i = 2; i <= n; i++)
		f *= i;

	return f;
}

// Sum all the elements in two matrices and store them into the target matrix.
void add2mat(double tar_mat[][40], const double mat1[][40], const double mat2[][40], const unsigned len)
{
	for(unsigned i = 0; i < 4; i++)
	{
		for(unsigned j = 0; j < len; j++)
			tar_mat[i][j] = mat1[i][j] + mat2[i][j];
	}
}

// Take average of weights after burn-in and store them into genes' information.
void avg_prob(VGene& vGene, const vector<double>& vProb)
{
	for(unsigned i = 0; i < vGene.size(); i++)
		vGene[i].prob = vProb[i]/WINDOW;
}

// Choose the minimum frequency of a base for a given column.
inline double ran_freq(const double psfm[][40], unsigned col)
{
	return log(psfm[(unsigned)floor(4*((double)rand()/RAND_MAX))][col]);
}

unsigned numofN(const string motif) // Count number of 'N' in motif.
{
	unsigned num = 0;
	for(unsigned i = 0; i < motif.length(); i++)
		if(motif[i] == 'N')
			num++;

	return num;
}

// sum up the results of all iterations.
void sum_iter(const unsigned i, VGene& vGene)
{
	if(i == 0)
	{
		for(unsigned j = 0; j < vGene.size(); j++)
			vGene[j].avg_prob = vGene[j].prob;
	}
	else
	{
		for(unsigned j = 0; j < vGene.size(); j++)
			vGene[j].avg_prob += vGene[j].prob;
	}
}

// take the average of all iterations and determine gene label.
void avg_iter(VGene& vGene)
{
	for(unsigned i = 0; i < vGene.size(); i++)
	{
		vGene[i].avg_prob /= STA_PNT;
		if(vGene[i].avg_prob > THRLD)
			vGene[i].target = true;
		else
			vGene[i].target = false;
	}
}

//******************************************************************************

double d_huge ( void )

//******************************************************************************
//
//  Purpose:
//
//    D_HUGE returns a "huge" real value, usually the largest legal real.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal real number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double D_HUGE, a "huge" real value.
//
{
  return HUGE_VAL;
}

//******************************************************************************

double d_epsilon ( void )

//******************************************************************************
//
//  Purpose:
//
//    D_EPSILON returns the round off unit for floating arithmetic.
//
//  Discussion:
//
//    D_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double D_EPSILON, the floating point round-off unit.
//
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}

//******************************************************************************

double gamma ( double x )

//******************************************************************************
//
//  Purpose:
//
//    GAMMA calculates the Gamma function for a real argument X.
//
//  Definition:
//
//    GAMMA(X) = Integral ( 0 <= T <= Infinity ) T**(X-1) EXP(-T) DT
//
//  Recursion:
//
//    GAMMA(X+1) = X * GAMMA(X)
//
//  Special values:
//
//    GAMMA(0.5) = SQRT(PI)
//    If N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
//
//  Discussion:
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the GAMMA
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for X .GE. 12 are from reference 2.
//    The accuracy achieved depends on the arithmetic system, the
//    compiler, the intrinsic functions, and proper selection of the
//    machine-dependent constants.
//
//  Machine-dependent constants:
//
//    BETA: radix for the floating-point representation.
//    MAXEXP: the smallest positive power of BETA that overflows.
//    XBIG: the largest argument for which GAMMA(X) is representable
//      in the machine, i.e., the solution to the equation
//      GAMMA(XBIG) = BETA**MAXEXP.
//    XMININ: the smallest positive floating-point number such that
//      1/XMININ is machine representable.
//
//    Approximate values for some important machines are:
//
//                               BETA       MAXEXP        XBIG
//
//    CRAY-1         (S.P.)        2         8191        966.961
//    Cyber 180/855
//      under NOS    (S.P.)        2         1070        177.803
//    IEEE (IBM/XT,
//      SUN, etc.)   (S.P.)        2          128        35.040
//    IEEE (IBM/XT,
//      SUN, etc.)   (D.P.)        2         1024        171.624
//    IBM 3033       (D.P.)       16           63        57.574
//    VAX D-Format   (D.P.)        2          127        34.844
//    VAX G-Format   (D.P.)        2         1023        171.489
//
//                              XMININ
//
//    CRAY-1         (S.P.)   1.84E-2466
//    Cyber 180/855
//      under NOS    (S.P.)   3.14E-294
//    IEEE (IBM/XT,
//      SUN, etc.)   (S.P.)   1.18E-38
//    IEEE (IBM/XT,
//      SUN, etc.)   (D.P.)   2.23D-308
//    IBM 3033       (D.P.)   1.39D-76
//    VAX D-Format   (D.P.)   5.88D-39
//    VAX G-Format   (D.P.)   1.12D-308
//
//  Author: 
//
//    William Cody and L. Stoltz,
//    Applied Mathematics Division,
//    Argonne National Laboratory,
//    Argonne, Illinois, 60439.
//
//    C++ translation by John Burkardt.
//
//  Reference: 
//
//    Willilam Cody,
//    "An Overview of Software Development for Special Functions", 
//    Lecture Notes in Mathematics, 506, 
//    Numerical Analysis Dundee, 1975, 
//    G. A. Watson (ed.),
//    Springer Verlag, Berlin, 1976.
//
//    Hart, Ward Cheney, Charles Lawson, Maehly, 
//    Charles Mesztenyi, John Rice, Thacher, Witzgall,
//    Computer Approximations, 
//    Wiley and sons, New York, 1968.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double GAMMA, the value of the function. 
//    The computation is believed to be free of underflow and overflow.
//
{
  double c[7] = {
    -1.910444077728E-03, 
     8.4171387781295E-04, 
    -5.952379913043012E-04, 
     7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
     8.333333333333333331554247E-02, 
     5.7083835261E-03 };
  double fact;
  int i;
  int n;
  double p[8] = {
    -1.71618513886549492533811, 
     2.47656508055759199108314E+01, 
    -3.79804256470945635097577E+02, 
     6.29331155312818442661052E+02, 
     8.66966202790413211295064E+02, 
    -3.14512729688483675254357E+04, 
    -3.61444134186911729807069E+04, 
     6.64561438202405440627855E+04 };
  bool parity;
  double pi = 
    3.141592653589793;
  double q[8] = {
    -3.08402300119738975254353E+01, 
     3.15350626979604161529144E+02,
    -1.01515636749021914166146E+03,
     -3.10777167157231109440444E+03, 
     2.25381184209801510330112E+04, 
     4.75584627752788110767815E+03, 
    -1.34659959864969306392456E+05, 
    -1.15132259675553483497211E+05 };
  double sqrtpi = 0.9189385332046727417803297;
  double sum2;
  double value;
  double xbig = 35.040;
  double xden;
  double xminin = 1.18E-38;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = false;
  fact = 1.0;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= 0.0 )
  {
    y = -x;
    y1 = ( double ) ( ( int ) ( y ) );
    value = y - y1;

    if ( value != 0.0 )
    {
      if ( y1 != ( double ) ( ( int ) ( y1 * 0.5 ) ) * 2.0 )
      {
        parity = true;
      }

      fact = -pi / sin ( pi * value );
      y = y + 1.0;
    }
    else
    {
      value = d_huge ( );
      return value;
    }

  }
//
//  Argument < EPS
//
  if ( y < d_epsilon ( ) )
  {
    if ( xminin <= y )
    {
      value = 1.0 / y;
    }
    else
    {
      value = d_huge ( );
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
//
//  0.0 < argument < 1.0
//
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
//
//  1.0 < argument < 12.0, reduce argument if necessary.
//
    else
    {
      n = int ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }

    value = xnum / xden + 1.0;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      value = value / y1;
    }
//
//  Adjust result for case  2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        value = value * y;
        y = y + 1.0;
      }

    }
  }
//
//  Evaluate for 12 <= argument.
//
  else
  {
    if ( y <= xbig )
    {
      ysq = y * y;
      sum2 = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum2 = sum2 / ysq + c[i];
      }
      sum2 = sum2 / y - y + sqrtpi;
      sum2 = sum2 + ( y - 0.5 ) * log ( y );
      value = exp ( sum2 );
    }
    else
    {
      value = d_huge ( );
      return value;

    }

  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    value = -value;
  }

  if ( fact != 1.0 )
  {
    value = fact / value;
  }

  return value;
}

//******************************************************************************

double student_pdf ( double x, double a, double b, double c )

//******************************************************************************
//
//  Purpose:
//
//    STUDENT_PDF evaluates the central Student T PDF.
//
//  Formula:
//
//    PDF(A,B,C;X) = Gamma ( (C+1)/2 ) /
//      ( Gamma ( C / 2 ) * Sqrt ( PI * C ) 
//      * ( 1 + ((X-A)/B)**2/C )**((C + 1)/2) )
//
//  Modified:
//
//    02 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, shape parameters of the PDF,
//    used to transform the argument X to a shifted and scaled 
//    value Y = ( X - A ) / B.  It is required that B be nonzero.
//    For the standard distribution, A = 0 and B = 1.
//
//    Input, double C, is usually called the number of 
//    degrees of freedom of the distribution.  C is typically an 
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
//    Output, double STUDENT_PDF, the value of the PDF.
//
{
# define PI 3.141592653589793

  double pdf;
  double y;

  y = ( x - a ) / b;

  pdf = gamma ( 0.5 * ( c + 1.0 ) ) / ( sqrt ( PI * c ) 
    * gamma ( 0.5 * c ) 
    * sqrt ( pow ( ( 1.0 + y * y / c ), ( c + 1.0 ) ) ) );
  
  return pdf;
# undef PI
}

//******************************************************************************

double normal_pdf ( double x, double a, double b )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL_PDF evaluates the Normal PDF.
//
//  Formula:
//
//    PDF(A,B;X) 
//      = exp ( - 0.5 * ( ( X - A ) / B )**2 )
//      / ( B * SQRT ( 2 * PI ) )
//
//  Discussion:
//
//    The normal PDF is also known as the Gaussian PDF.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
# define PI 3.141592653589793

  double pdf;
  double y;

  y = ( x - a ) / b;

  pdf = exp ( -0.5 * y * y )  / ( b * sqrt ( 2.0 * PI ) );

  return pdf;
# undef PI
}
