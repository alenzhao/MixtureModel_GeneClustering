#include "mysub.h"

// All subroutines used in Gibbs Sampling project.

// Transform a string to upper case.
string& str2upper(string& str)
{
	for(unsigned int i = 0; i < str.length(); i++)
		str[i] = toupper(str[i]);

	return str;
}

// Transform a string to lower case.
string& str2lower(string& str)
{
	for(unsigned int i = 0; i < str.length(); i++)
		str[i] = tolower(str[i]);

	return str;
}

// Find motif location with mismatch.
VLoc locmotif(const string motif, const unsigned FLANKING, const string seq, const unsigned int misallow)
{
	VLoc vLoc;
	for(unsigned i = FLANKING; i+motif.length()+FLANKING < seq.length(); i++)
	{
		unsigned j, mismatch = 0;
		for(j = 0; j < motif.length(); j++)
		{
			if(!iupac(seq[i+j], motif[j]))
				mismatch++;
			if(mismatch > misallow)
				break;
		}
		if(j == motif.length())
		{
			MatchSite match;
			match.loc = i;
			match.mis = mismatch;
			vLoc.push_back(match);
		}
	}

	return vLoc;	
}

// Extract the extended motif from DNA sequence.
string emotif(const int startLoc, const SST mLen, const unsigned FLANKING, const string& seq)
{
	SST eLen = mLen + 2*FLANKING; // extended motif length.
	int eStart = startLoc - FLANKING;

	if(eStart < 0 || eStart + eLen > seq.length())
	{
		cerr << "Flanking region exceeds boundary!" << endl;
		return "";
	}
	
	return seq.substr(eStart, eLen);
}

// Copy reversed sequence into another string.
string copyrev(const string seq)
{
	string rev = seq;
	for(unsigned int i = 0; i < seq.length(); i++)
		rev[i] = basepair(seq[seq.length()-1-i]);
	
	return rev;
}


char basepair(char base) // base-pairing function.
{
	switch(base)
	{
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		case 'a':
			return 't';
		case 't':
			return 'a';
		case 'g':
			return 'c';
		case 'c':
			return 'g';
		default:
			return 'N';
	}
}

bool cmp(OrfMotif a, OrfMotif b) // Comparison function for sorting target genes.
{
	if(a.avg_prob == b.avg_prob)
		return a.val > b.val;
	else
		return a.avg_prob > b.avg_prob;
}


void write_psfm(ofstream& h, const unsigned len, const double psfm[][40], const double info[40]) // Write PSFM into file.
{
	for(unsigned i = 0; i < len; i++)
	{
		h << i+1 << "\t";
		h.width(12); h.setf(ios::left); h.setf(ios::showpoint); // format for each element of matrix.
		h << psfm[dnaidx('A')][i];
		h.width(12); h.setf(ios::left); h.setf(ios::showpoint);
		h << psfm[dnaidx('C')][i];
		h.width(12); h.setf(ios::left); h.setf(ios::showpoint);
		h << psfm[dnaidx('G')][i];
		h.width(12); h.setf(ios::left); h.setf(ios::showpoint);
		h << psfm[dnaidx('T')][i];
		if(info != 0)
			h << info[i] << endl;
		else
			h << endl;
	}
}

void write_gene(ofstream& h, const VGene& vGene, const bool target, GeneAnno& geneAnno, map<string, string>& mbindloc) // Write gene information into file.
{
	if(target)
	{
		for(unsigned i = 0; i < vGene.size(); i++)
			if(vGene[i].target)
			{
				h.width(10); h.setf(ios::left); // format for ORF name.
				h << vGene[i].name << "\t";
				GeneInfo geneInfo;
				if(!geneAnno.empty())
				{
					geneInfo = geneAnno[vGene[i].name];
					string genNam = geneInfo.geneName;
					h.width(10); h.setf(ios::left); // format for gene name.
					h << genNam << "\t";
				}
				h.width(10); h.setf(ios::left); h.setf(ios::showpoint); // format for probability.
				h << vGene[i].avg_prob << "\t";
				h << vGene[i].seq << "\t";
				h << vGene[i].loc << '\t';
				h << vGene[i].orien << '\t';
				h << vGene[i].val << "\t\t";
				if(!geneAnno.empty())
				{
					string genDesc = geneInfo.geneDesc;
					string refSeq = geneInfo.refSeq;
					string locusLink = geneInfo.locusLink;
					h << genDesc << "\t\t\t";
					h << refSeq << "\t";
					h << locusLink << "\t";
				}
				h << mbindloc[vGene[i].name];
				h << endl;
			}
	}
	else
	{
		for(unsigned i = 0; i < vGene.size(); i++)
			if(!vGene[i].target)
			{
				h.width(10); h.setf(ios::left); // format for ORF name.
				h << vGene[i].name << "\t";
				GeneInfo geneInfo;
				if(!geneAnno.empty())
				{
					geneInfo = geneAnno[vGene[i].name];
					string genNam = geneInfo.geneName;
					h.width(10); h.setf(ios::left); // format for gene name.
					h << genNam << "\t";
				}
				h.width(10); h.setf(ios::left); h.setf(ios::showpoint); // format for probability.
				h << vGene[i].avg_prob << "\t";
				h << vGene[i].seq << "\t";
				h << vGene[i].loc << '\t';
				h << vGene[i].orien << '\t';
				h << vGene[i].val << "\t\t";
				if(!geneAnno.empty())
				{
					string genDesc = geneInfo.geneDesc;
					string refSeq = geneInfo.refSeq;
					string locusLink = geneInfo.locusLink;
					h << genDesc << "\t\t\t";
					h << refSeq << "\t";
					h << locusLink << "\t";
				}
				h << mbindloc[vGene[i].name];
				h << endl;
			}
	}
}
	
// Calculate information content for PSFM.
void entropy(const double matrix[][40], const unsigned len, double info[40])
{
	for(unsigned i = 0; i < len; i++)
	{
		info[i] = 2.0;
		for(unsigned j = 0; j < 4; j++)
		{
			if(matrix[j][i] == 0.0)
				continue;
			info[i] += matrix[j][i]*log(matrix[j][i])/log((double)2);
		}
	}
}

// Calculate profile frequency ratio.
void prof_ratio(const unsigned eLen, const double tMat[][40], const double nMat[][40], double rMat[][40])
{
	for(unsigned i = 0; i < 4; i++)
		for(unsigned j = 0; j < eLen; j++)
		{
			if(nMat[i][j] != 0.0)
				rMat[i][j] = tMat[i][j]/nMat[i][j];
			else if(tMat[i][j] != 0.0)
				rMat[i][j] = -1.0;
			else
				rMat[i][j] = 1.0;
		}
}

int load_data(const string& fMotif, const string& fExpr, const string& fSeq, const string& fGeneAnno, string& motifSeq, Expr& expr, Seq& seq, GeneAnno& geneAnno) // Load data into memory by data name.
{
	// Create handle to motif file and open it.
	ifstream hMotif(fMotif.data());
	if(!hMotif)
	{
		cerr << "Open " << fMotif << " failed!" << endl;
		return 1;
	}
	// Read motif into memory.
	if(hMotif.good())
	{
		hMotif >> motifSeq;
		str2upper(motifSeq);
	}
	else
	{
		cerr << "Read motif error!" << endl;
		return 1;
	}
	hMotif.close();

	// Create handle to expression file and open it.
	ifstream hExpr(fExpr.data());
	if(!hExpr)
	{
		cerr << "Open " << fExpr << " failed!" << endl;
		return 1;
	}
	// Create hash table for expression data.
	while(hExpr.good())
	{
		string orfNam;
		double exprVal;
		if(hExpr.peek() == -1)
			break;
		hExpr >> orfNam; str2upper(orfNam);
		if(hExpr.peek() == '\n')
			exprVal = 0;
		else
			hExpr >> exprVal;
		hExpr.ignore(1, '\n');
		if(orfNam != "")
		{
			Expr::iterator eIter = expr.find(orfNam); // Determine whether the gene has already existed.
			if(eIter == expr.end()) // If yes, insert probe name and value pair.
				expr[orfNam] = exprVal;
			else if(exprVal > expr[orfNam]) // Else, modify probe value if the new one is larger than the existing.
				expr[orfNam] = exprVal;
		}
	}
	hExpr.close();

	// Create handle to sequence file and open it.
	ifstream hSeq(fSeq.data());
	if(!hSeq)
	{
		cerr << "Open " << fSeq << " failed!" << endl;
		return 1;
	}
	// Create hash table for sequence data.
	while(hSeq.good())
	{
		string orfNam, orfSeq;
		if(hSeq.peek() == -1)
			break;
		hSeq >> orfNam; str2upper(orfNam);
		hSeq >> orfSeq; //str2upper(orfSeq);
		seq[orfNam] = orfSeq;

		hSeq.ignore(1, '\n');
	}
	hSeq.close();

	// Create handle to gene annotation file and open it.
	if(fGeneAnno == "")
		geneAnno.clear();
	else
	{
		ifstream hGeneAnno(fGeneAnno.data());
		if(!hGeneAnno)
		{
			cerr << "Open " << fGeneAnno << " failed!" << endl;
			return 1;
		}
		// Create hash table for gene annotation data.
		while(hGeneAnno.good())
		{
			string strLn;
			getline(hGeneAnno, strLn);
			// Data are separated by tabs. RefSeq\tGene name\tUniGeneID\tGene anno\tLocuslink\n.
			SST fstTab = strLn.find('\t');
			SST sndTab = strLn.find('\t', fstTab + 1);
			SST trdTab = strLn.find('\t', sndTab + 1);
			SST fohTab = strLn.find('\t', trdTab + 1);
			SST lineFd = strLn.find('\n', fohTab + 1);
			string refSeq = strLn.substr(0, fstTab);
			string genNam = strLn.substr(fstTab + 1, sndTab - fstTab - 1);
			string orfNam = strLn.substr(sndTab + 1, trdTab - sndTab - 1); str2upper(orfNam);
			string genDesc = strLn.substr(trdTab + 1, fohTab - trdTab - 1);
			string locusLink = strLn.substr(fohTab + 1, lineFd - fohTab - 1);
			if(orfNam != "")
			{
				GeneInfo geneInfo;
				geneInfo.refSeq = refSeq;
				geneInfo.geneName = genNam;
				geneInfo.geneDesc = genDesc;
				geneInfo.locusLink = locusLink;
				geneAnno[orfNam] = geneInfo;
			}
		}
		hGeneAnno.close();
	}

	return 0;
}

// Find candidate genes that contain motif with mismatch.
void motifgene(const string& motifSeq, const unsigned FLANKING, Expr& expr, const Seq& seq, const unsigned mis, VGene& vGene)
{
	OrfMotif aGene;
	unsigned nskip = 0;
	VLoc::size_type maxBinding = 0;
	string maxGene;
	for(Seq::const_iterator ci = seq.begin(); ci != seq.end(); ci++)
	{
		// ci->first: ORF name, ci->second: DNA sequence.
		aGene.locOrg = locmotif(motifSeq, FLANKING, ci->second, mis);
		string revseq = copyrev(ci->second);
		aGene.locRev = locmotif(motifSeq, FLANKING, revseq, mis);
		if(aGene.locOrg.empty() && aGene.locRev.empty()) // No motif contained, skip.
			continue;
		VLoc::size_type nBinding = aGene.locOrg.size() + aGene.locRev.size();
		if(nBinding > maxBinding)
		{
			maxBinding = nBinding;
			maxGene = aGene.name;
		}
		if(aGene.locOrg.empty()) // Only on reverse sequence.
			aGene.seq = emotif(aGene.locRev[0].loc, motifSeq.length(), FLANKING, revseq);	
		else // On original sequence.
			aGene.seq = emotif(aGene.locOrg[0].loc, motifSeq.length(), FLANKING, ci->second); // extract extended motif sequence.
		aGene.name = ci->first;
		if(expr.find(ci->first) == expr.end())
		{
			if(nskip == 0)
			{
				cerr << "Warning: expression ratios not found, genes skipped." << endl;
				nskip++;
			}
			else
				nskip++;
			continue;
		}
		aGene.val = expr[ci->first];
		vGene.push_back(aGene); // store an object of OrfMotif into vector.
	}
	if(nskip > 0)
		cerr << "Total genes skipped: " << nskip << endl;
	if(VERBOSE)
		cout << "Gene: " << maxGene << " has maximum number, " << maxBinding << ", of binding sites." << endl;
}

void motifgene(Expr& expr, const Seq& seq, VGene& vGene) // Load all genes.
{
	OrfMotif aGene;
	for(Seq::const_iterator ci = seq.begin(); ci != seq.end(); ci++)
	{
		aGene.name = ci->first;
		aGene.val = expr[ci->first];
		aGene.seq = ci->second;
		vGene.push_back(aGene);
	}
}

// generate initial labels for candidate genes.
void genlbl(VGene& vGene, const double bkm, const double bkstd)
{
	srand((unsigned)time(NULL));
	// Generate labels according to expression values.
	for(unsigned i = 0; i < vGene.size(); i++)
	{
		double p0, p1;
		p0 = normal_pdf(vGene[i].val, bkm, bkstd);
		p1 = normal_pdf(vGene[i].val, bkm+::INIV*bkstd, bkstd);
		vGene[i].label = Bernoulli(p1/(p0+p1));
	}
}

void sudocnt(const Seq& seq, double portion[4]) // Calculate pseudo-count of DNA bases.
{
	for(unsigned i = 0; i < 4; i++)
		portion[i] = 0.0;
	for(Seq::const_iterator ci = seq.begin(); ci != seq.end(); ci++)
	{
		for(unsigned i = 0; i < ci->second.length(); i++)
		{
			if(ci->second[i] == 'N' || islower(ci->second[i]))
				continue;
			portion[dnaidx(ci->second[i])]++;
		}
	}
	double sumBase = 0.0;
	for(unsigned i = 0; i < 4; i++)
		sumBase += portion[i];
	for(unsigned i = 0; i < 4; i++)
		portion[i] /= sumBase;
}

// Write both target and non-target PSFMs, and profile ratios into log file.
int write_log(const string logname, const unsigned len, double tPSFM[][40], double nPSFM[][40])
{
	ofstream h(logname.data()); // Open file for writing.
	if(!h)
	{
		cerr << "Open " << logname << " failed!" << endl;
		return 1;
	}

	// Write header.	
	h << "\tA\tC\tG\tT\tEntropy" << endl << endl;
	double tInfo[40], nInfo[40]; // Entropy.

	// Write target PSFM.
	h << "Target gene PSFM: " << endl << endl;
	entropy(tPSFM, len, tInfo); // Calculate entropy.
	write_psfm(h, len, tPSFM, tInfo); // Write target PSFM.
	h << endl << endl;

	// Write non-target PSFM.	
	h << "Non-target gene PSFM: " << endl << endl;
	entropy(nPSFM, len, nInfo); // Calculate entropy.
	write_psfm(h, len, nPSFM, nInfo); // Write non-target PSFM.
	h << endl << endl;
	
	// Write profile ratios.
	h << "Profile ratios: " << endl << endl;
	double profMat[4][40];
	prof_ratio(len, tPSFM, nPSFM, profMat);
	write_psfm(h, len, profMat); // Write non-target PSFM.
	h.close();
	
	return 0;
}
	
void write_log(ofstream& log_h, const unsigned len, double tPSFM[][40], double nPSFM[][40])
{
	// Write header.	
	log_h << "\n\tA\tC\tG\tT\tEntropy" << endl << endl;
	double tInfo[40], nInfo[40]; // Entropy.

	// Write target PSFM.
	log_h << "Target gene PSFM: " << endl << endl;
	entropy(tPSFM, len, tInfo); // Calculate entropy.
	write_psfm(log_h, len, tPSFM, tInfo); // Write target PSFM.
	log_h << endl << endl;

	// Write non-target PSFM.	
	log_h << "Non-target gene PSFM: " << endl << endl;
	entropy(nPSFM, len, nInfo); // Calculate entropy.
	write_psfm(log_h, len, nPSFM, nInfo); // Write non-target PSFM.
	log_h << endl << endl;
	
	// Write profile ratios.
	log_h << "Profile ratios: " << endl << endl;
	double profMat[4][40];
	prof_ratio(len, tPSFM, nPSFM, profMat);
	write_psfm(log_h, len, profMat); // Write profile ratios.
}

// Write both target and non-target genes information.
int write_info(const string infoname, const VGene& vGene, GeneAnno& geneAnno, const double tm, const double tstd, 
			   const double bkm, const double bks, const unsigned nTotal, const unsigned nTar, const unsigned nNon, 
			   const HypePrio& hypePrio, const string motif, map<string, string>& mbindloc)
{
	ofstream h(infoname.data()); // Open file for writing.
	if(!h)
	{
		cerr << "Open " << infoname << " failed!" << endl;
		return 1;
	}

	// Write statistics for expression data.
	h << "Total number of genes: " << nTotal << endl;
	h << "Total candidate genes: " << vGene.size() << endl;
	h << "Number of targets: " << nTar << endl;
	h << "Number of non-targets: " << nNon << endl;
	h << "All genes mean: " << bkm << ", std: " << bks << endl;
	h << "Target mean: " << tm << ", std: " << tstd << endl;
	h << "Core motif sequence: " << motif << endl << endl;
	h << "---- Information of hype-parameters ----" << endl;
	h << "Hypothetical samples: " << hypePrio.V << endl;
	h << "Hypothetical mean: " << hypePrio.MU << endl;
	h << "Hypothetical variance: " << hypePrio.SIGMA_2 << endl;
	h << "Flanking region: " << hypePrio.FLANKING << endl;
	h << "---- Gibbs Sampler Information ----" << endl;
	h << "Burn-in period: " << BURN_IN << ", Window size: " << WINDOW << endl;
	h << "Starting points: " << STA_PNT << endl << endl;
	h << "---- Information of some other settings ----" << endl;
	h << "Standard deviations of enforcement: " << ENFORN << endl;
	h << "Standard deviations of initial labels: " << INIV << endl;
	h << "Skip core motif region when calculate sequence likelihood? " << SKIPCORE << endl;
	h << "Randomly pick a binding site for low ratio genes? " << RANDBG << endl << endl;

	// Write header.
	if(geneAnno.empty())
		h << "ID\tProbability\tSeq\tLocation\tOrientation\tLog-ratio\tTentative binding sites" << endl << endl;
	else
		h << "ID\tName\tProbability\tSeq\tLocation\tOrientation\tLog-ratio\tAnnotation\tRefSeq\tLocusLink\tTentative binding sites" << endl << endl;

	// Write target PSFM.
	h << "Target gene information: " << endl << endl;
	write_gene(h, vGene, true, geneAnno, mbindloc);
	h << endl << endl;

	// Write non-target PSFM.
	h << "Non-target gene information: " << endl << endl;
	write_gene(h, vGene, false, geneAnno, mbindloc);
	h.close();
	
	return 0;
}

void cal_chars(VGene& vGene, Chars& charsInfo) // Calculate TP, FP etc. characteristics of Gibbs sampler.
{
	charsInfo.FN = 0;
	charsInfo.FP = 0;
	charsInfo.TN = 0;
	charsInfo.TP = 0;

	for(unsigned i = 0; i < vGene.size(); i++)
	{
		char* dumb;
		unsigned long label = strtoul(vGene[i].name.data(), &dumb, 10);
		if(label <= 1000 && vGene[i].target)
			charsInfo.TP++;
		else if(label <= 1000 && !vGene[i].target)
			charsInfo.FN++;
		else if(label > 1000 && vGene[i].target)
			charsInfo.FP++;
		else if(label > 1000 && !vGene[i].target)
			charsInfo.TN++;
	}
}

bool iupac(const char nuc, const char pac) // Test IUPAC symbol match.
{
	if(nuc == 'N')
		return false;

	switch(pac)
	{
		case 'A':
		case 'C':
		case 'G':
		case 'T':
			return (nuc == pac);
		case 'R':
			if(nuc == 'A' || nuc == 'G')
				return true;
			else
				break;
		case 'Y':
			if(nuc == 'C' || nuc == 'T')
				return true;
			else
				break;
		case 'M':
			if(nuc == 'A' || nuc == 'C')
				return true;
			else
				break;
		case 'K':
			if(nuc == 'G' || nuc == 'T')
				return true;
			else
				break;
		case 'S':
			if(nuc == 'C' || nuc == 'G')
				return true;
			else
				break;
		case 'W':
			if(nuc == 'A' || nuc == 'T')
				return true;
			else
				break;
		case 'B':
			if(nuc != 'A')
				return true;
			else
				break;
		case 'D':
			if(nuc != 'C')
				return true;
			else
				break;
		case 'H':
			if(nuc != 'G')
				return true;
			else
				break;
		case 'V':
			if(nuc != 'T')
				return true;
			else
				break;
		case 'N':
			return true;
		default:
			return false;
	}

	return false;
}

