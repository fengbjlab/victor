// body file: libfbj_genepi.cpp

#ifndef LIBFBJ_GENEPI
#define LIBFBJ_GENEPI

#include <vector>
#include <map>
#include <set>
#include <deque>
#include <string>

namespace genepi
{
	
	using namespace std;

	// ------------------ user defines ------------------
	// On 2017-07-17, changed chr_num from PLINK coding (1-22,X,Y,XY,M) to alphanumeric coding (1-22,M,X,Y,XY). XY is mostly not used so it's the last.
	// this section require set_path() or read_arguments()

	extern string MainPath; // {"./"}
	void set_path(string& MainPath, string Genome=""); // Genome can be hg19/hg38/GRCh37/GRCh38, case insensitive. Genome can be empty if MainPath contains one but only one of them.
	string get_genome();
	void read_arguments(vector<string>& srce_opt, bool to_set_path);
	void add_help_text_var();
	string help_text();
	
	int chrlen_bp(int chrNum);	// chrlen_bp[chrNum]
	int chrlen_bp_cumulative(int chrNum);
	int total_genome_len();
	int MaxChrNum();
	std::string cytoband(const int chrNum, const int bp1, const int bp2); // assumes bp1<=bp2
	int expected_ploidy(int chr_num, int bp, int sex);
	bool is_chrX		(int chr_num);
	bool is_chrY		(int chr_num);
	bool is_chrXY		(int chr_num);
	bool is_chrM		(int chr_num);
	bool is_autosomal	(int chr_num);
	bool is_autoOrPAR	(int chr_num, int bp);
	bool is_MHC			(int chr_num, int bp);
	string DNA_seq(int chrNumPlink, int bp, int len);
	bool read_chr_str(const string& input, string& cs);		// if success, cs = standard name (01-22,X,Y,XY,MT). return successfulness.
	int read_chr_num(const string& input);					// if success, return [1,kMaxChrNum]; otherwise, 0.
	bool read_chr_num(const string& input, int& chrNum);	// if success, chrNum=[1,kMaxChrNum]; otherwise, 0.
	string convert_chr_num(int chr_num);					// if success, return chr_str, otherwise exit_error
	
	// ------------------ ChrRegions ------------------

	// After reading, all positions are 1-based. Three ways to input chr regions:
	// 1) constructor with filename
	// 2) setup()
	// 3) add() and then finishing()
	// 4) set_whole_genome()
	// zero_based is only for input from a BED file (3 columns: Chr Start End). It's better to be true so that it's consistant with BED.
	// omit_unk conventions: false for internal functions in this library (exit w/ error msg), true for external useage by the user (no error msg, keep running).
	class ChrRegions {
	private:
		std::map< int, std::map< int,int> >	regions;
		std::map< int, std::vector< std::pair<int,int> > > unfinished;
		bool whole;
		bool did_setup;
	public:
		ChrRegions():whole(false),did_setup(false) {}
		ChrRegions(const             string&  filename,bool zero_based,int padding,bool omit_unk) { setup(filename,zero_based,padding,omit_unk); }
		ChrRegions(const std::vector<string>& filename,bool zero_based,int padding,bool omit_unk) { setup(filename,zero_based,padding,omit_unk); }
		ChrRegions(const std::deque <string>& filename,bool zero_based,int padding,bool omit_unk) { setup(filename,zero_based,padding,omit_unk); }
		ChrRegions(const std::set   <string>& filename,bool zero_based,int padding,bool omit_unk) { setup(filename,zero_based,padding,omit_unk); }
		template <typename T>  void setup(const T& filename,bool zero_based,int padding, bool omit_unk); // same T as above; whole=false
		void set_whole_genome(bool w) { whole=w; }
		void add(const int chrNum, const int start, const int end); // start,end is 1-based
		bool add(const std::vector<string>& input,bool zero_based,int padding, bool omit_unk);	// add to unfinished. Must call finishing() after adding all sites.
		void finishing();											// convert unfinished to regions, then clear unfinished.
		
		std::map< int, std::map< int,int> > data() { return regions; }
		bool not_read();
		bool empty();
		bool contain(const int chrNumPlink, const int bp);
		bool overlap(const int chrNumPlink, const int bp1, const int bp2);
		void write(std::ostream& file, bool zero_based); // add zero_based on 2018-04-24, previously always true
		void write(int desert_length_in_bp, const string& output_prefix);
		int  total_bp();
	};

	// ------------------ genes ------------------
	extern string GDB_name;		// gene database name (file is GDB_name.txt) {refGeneLite}
	extern bool	canonic;		// read canonical tx only {false}
	extern bool filt_tx;		// filter transcripts {false}
	extern set<string>	incl_g;	// gene symbol filter {}
	extern set<string>	incl_t;	// transcript ID filter {}

	// Same as UCSC, all start position is 0-based while end position are 1-based. chr is one character.
	// cdsStat is priority code about completeness of the gene. Bigger is better.
	// lenStat is priority by cdsLen & txLen.
	// name=tx_id name2=gene_id name3=tx_CHR#_txInit name4=gene_CHR#_txInit. for refGene, gene_id is also symbol.
	// exonStarts, exonEnds, xSrelative, xErelative start from the smallest genomic location, independent of the orientation of the gene.
	// for ncRNA cdsStart=cdsEnd=TxEnd. So I use (cdsStart==cdsEnd) to check whether it's a non-coding RNA. Also, for ncRNA u5Len=txLen.
	typedef std::pair<int,int> lenStat_t;
	struct gene_info {
		char strand; // + -
		int txStart, txEnd, txLen, txRgn, cdsStart, cdsEnd, cdsLen, cdsRgn, exonCount, u5Len, u3Len;
		int cdsStat;
		int chrNumPlink;
		int cdsOffset; // 0 does nothing, good for non-coding genes or coding genes with a start codon. For coding genes w/o a start codon, it's 0/1/2.
		string chr, name, name2, name3, name4, symbol, biotype, cdsStartStat, cdsEndStat, proteinID, alignID, col16;
		lenStat_t lenStat;
		vector<int> exonStarts, exonEnds;
		vector<int> xSrelative, xErelative; // exonStart/end relative to c.1 (A of the ATG) for coding genes / to n.1 for ncRNA
		gene_info(){}
	};

	typedef multimap< lenStat_t, gene_info, std::greater<lenStat_t> > tx_set_t;	// tx_set[lenStat]=tx
	extern map<string,gene_info> gene_db;				// by transcript name, bad because gene dup are randomly removed, but good if I need unique info for each TxID.
	extern map< string, tx_set_t > gene_bySymbol;		// [gene symbol]=set, tx_set_t may not be unique even with --canonical due to gene duplications (1 TxID 2+ Loc)
	extern map< string, tx_set_t > gene_byTranscript;	// [transcript]=set, tx_set_t may not be unique even with --canonical due to gene duplications (1 TxID 2+ Loc)
	
	void read_genes();
	void read_refGene   (const string& filename, bool read_tx=false, bool filter=true); // UCSC:refGene.txt or UCSC:ensGene.txt
	void read_knownGene (      string  pathname, bool read_tx=false); // read knownGene.txt, requires kgXref.txt knownCanonical.txt
	void read_UCSCgenes (	   string  pathname, const string& type="refGene", bool read_tx=false, bool filter=true);
	void write_knownGene(const string& filename, bool sort_tx_len=false);
	void write_universal(const string& filename, bool sort_tx_len=false);

	// annotation
	int genomic_to_cds_location(const gene_info& g, int chr, int bp); // 0:intronic, -std::numeric_limits<int>::max(): not in Tx or not protein-coding.
	int genomic_to_rna_location(const gene_info& g, int chr, int bp); // 0:intronic, -std::numeric_limits<int>::max(): not in Tx. No <0 for upstream for consistancy w/ cds.
	int cds_to_genomic_location(const gene_info& g, int loc, string& chr, int& bp, int& xS, int& xE); // input g loc, output chr bp xS xE
	int rna_to_genomic_location(const gene_info& g, int loc, string& chr, int& bp, int& xS, int& xE); // xS,xE is distance from exonStart,exonEnd (indep. of strand)
	inline int cds_to_genomic_location(const gene_info& g, int loc, string& chr, int& bp) { int xS, xE; return cds_to_genomic_location(g,loc,chr,bp,xS,xE); }
	inline int rna_to_genomic_location(const gene_info& g, int loc, string& chr, int& bp) { int xS, xE; return rna_to_genomic_location(g,loc,chr,bp,xS,xE); }
	bool overlap_StoE (const gene_info& g, int chr, int bp1, int bp2, int up=0, int dn=0);			// cover the whole gene from Start to End. Both bp1 bp2 are 1-based (same below).
	bool overlap_gene (const gene_info& g, int chr, int bp1, int bp2, int up=0, int dn=0);			// cover at least part of the gene.
	bool overlap_exon (const gene_info& g, int chr, int bp1, int bp2, int up=0, int dn=0, int in=0);// cover at least part of the exon.
	bool overlap_cds_ (const gene_info& g, int chr, int bp1, int bp2);								// cover at least part of the cds.
	bool in_first_exon(const gene_info& g, int chr, int bp);
	bool in_last_exon (const gene_info& g, int chr, int bp);
	bool nearest_copy(int chr, int bp, const string& TxID, gene_info* g);
	int read_hgvs(const string& hgvs, string& gene, char& code, string& chr, int& stt, int& end, string& fr, string& to, string& info);
	// hgvs is annotation for a gene. Format HGVS_Tx1,HGVS_Tx2[,]. HGVS_Tx# is Ann1:Ann2:etc. Ann# can be a Gene, TxID, [gmncrp].xxx.
	// If multiple Tx, read RefSeq. If multiple isoforms, read the canonical. If multiple Ann, read ncfrt3rgbtrvdzztqar4q431ee11q 6wyeysbg11v111r1ftyyyyyyyyyyyyiiiiuijrbgztydhgcgc5bgy3ghf6vgbvvb6ytgvgtxv rfrvtnr,ilbo. == c. >> p. == r. >> g. == m.
	// No multiple genes. To identify the gene, it reads RefSeq (NM_xxxx) > NCBI_Gene_Symbol > Ann1. Therefore, gene will become "g" if hgvs is g.A2311235T.
	// For p.xxx, reads SNV only.  If the gene is on the minor strand, fr/to not rc()!
	// For c.xxx, reads SNV/InDel. If the gene is on the minor strand, fr/to is rc(). It trys to check the RefGenome whenever possible.
	// It reads p.(=) p.0 p.0? p.Trp26Cys p.W26C p.W26= p.Trp26* p.W26* p.Trp26Ter p.*110Gln[ext*17] p.*110Q[ext*17] p.*321Arg[ext*?] ([] means optional).
	// It reads c.8158G>C c.-5G>C c.4909+1G>A c.8851-1G>T c.*5C>T c.7176_7177insT c.*8TA[4] c.*8_*9[4] c.1947dup c.185_186del c.203_506inv c.112_117delinsTG
	// It also reads the outdated HGVS c.IVS2+1G>C.
	
	string RNA_seq(const gene_info& g, bool convert_to_U); // previously mRNA_seq
	string UT5_seq(const gene_info& g, bool convert_to_U);
	string UT3_seq(const gene_info& g, bool convert_to_U);
	string CDS_seq(const gene_info& g);
	string AA_seq(const gene_info& g);

	// ------------------ DNA sequence ------------------
	
	extern std::map< string, char > AA_3to1;
	extern std::map< char, string > AA_1to3;
	extern std::map< string, char > genet_code;
	extern std::map< char, char > IUPAC_complement;
	string tr_AA_1to3(const string& in);
	string tr_AA_1to3(const char& in);
	string translate_cds(const string& dna_seq, bool stop_at_stop);

	// http://en.wikipedia.org/wiki/Nucleic_acid_notation exclude U because it's DNA code
	bool	dna_complement		(      string& seq, bool ExitIfErr);
	string	dna_complement_copy	(const string& seq, bool ExitIfErr);
	bool	dna_rc				(      string& seq, bool ExitIfErr);
	string	dna_rc_copy			(const string& seq, bool ExitIfErr);
	
	// ------------------ genotypes ------------------

	struct gtp_par {
		size_t	DP_col; // DP = read depth at this position for this sample (Integer)
		size_t	GQ_col; // GQ = conditional genotype quality = −10log10 p(genotype call is wrong, conditioned on the site’s being variant) (Integer)
		size_t	GP_col; // GP = genotype probability, three numbers in 0-1 per the latest VCF, not phred scale in the early VCF. Produced by IMPUTE2.
		gtp_par():DP_col(0),GQ_col(0),GP_col(0) {}
		void clear() { DP_col=0; GQ_col=0; GP_col=0; }
		void read(const string& s); // read location of DP,GQ,GP in FORMAT. Examples: GT:AD:DP:GQ:PL from RGC, GT:ADS:DS:GP from IMPUTE2
	};
	
	class genotype {
	private:
		static int		ReadOpt(vector<string>& opt, size_t& argi);// stop at the 1st unknown option; new argi point to its left
		static bool		is_missing(const char g);		// Don't use this, use !usable().
		static bool		has_missing(const string& in);	// Don't use this, use !usable().
		static char		to_char(const string& s, const bool to_filter, int& GQ, int& DP, double& GP, const gtp_par& par); // Don't use, use read(). 1=alt, a,0,2,3,..=ref. GQ,DP,GP=-1 if not read.
		static bool     vcf_gtp(const string& s);
		static int      GQ_tot;
		static int      DP_tot;
		static double	GP_tot;
		static int      GQ_cnt;
		static int      DP_cnt;
		static int      GP_cnt;

	public:
		// Below, characters should not overlap. Not even between diploid and haploid. Do not use a character >= 0x40.
		// A diploidy sequence looks like this: ===.==='===:===?===%===
		// A haploidy sequence looks like this: ---*---!-----
		static const char HoR = '=';	// diploidy 0|0
		static const char P01 = '.';	// diploidy 1|0
		static const char P10 = '\'';	// diploidy 0|1
		static const char HoA = ':';	// diploidy 1|1
		static const char Het = '%';	// diploidy 0/1
		static const char mss = '?';	// diploidy .|.
		static const char HpR = '-';	// haploidy 0
		static const char HpA = '*';	// haploidy 1
		static const char HpM = '!';	// haploidy .
		static const char CNA = ' ';	// chromosome not applicable

		static bool		NoMiss;
		static int		DP_cut; // DP default = 10
		static int		GQ_cut; // GQ default = 40
		static double	GP_cut; // GP default = 0.8
		static bool		GP_dom; // QC in dominant model,  can only be used for analysis under a dominant model
		static bool		GP_rec; // QC in recessive model, can only be used for analysis under a recessive model

		static int		mean_GQ();
		static int		mean_DP();
		static double	mean_GP();
		static void     rewind();
		
		static string	help_text();
		static void		read_arguments(vector<string>& arg, size_t start=1, bool expect_at_start=false, bool stop_at_unknown=false); // start is 0-based
		static double	prob_1(const char g);		// hap 1 (left  to |) is alt allele
		static double	prob1_(const char g);		// hap 2 (right to |) is alt allele
		static int		ploidy(const char g);		// return 0 (NA) 1 (haploid) 2 (diploid)
		static int		num_alt(const char g);				// number of alternative allele, HpA is 1
		static int		num_alt_recessive(const char g);	// number of alternative allele, HpA is 2
		static string	to_linkage(const char g);
		static string	to_vcf(const char g);
		static void		swap_allele(char& g);
		static void		set_all_ref(char& g);
		static void		set_missing(char& g);
		static void		dephase(char& g);
		
		static bool		valid(const char g);
		static bool		usable(const char g) { return !is_missing(g) &&  valid(g); }
		static bool		usable(const string& s) { return !s.empty() && !has_missing(s) && valid(s[0]); }
		static bool		phased(const string& s) { for (auto &c:s) if (c==Het) return false; return true; }
		static bool		denovo(const char g_ch, const char g_pa, const char g_ma, bool allow_hom=false);
		static char		fix_ploidy(const char g, const int to_ploidy, int& ploidy_err); // ++ploidy_err heterozygous haploidy
		static char		read(const string& s, const int chr_num, const int bp, const int sex, const gtp_par& par);	// make missing by DP,GQ,GP; fix ploidy; return g
		static char		read(const string& s, const gtp_par& par);													// make missing by DP,GQ,GP;             return g
		static char		read_vcf(const string& s, const int chr_num, const int bp, const int sex, const bool to_filter, int& GQ, int& DP, double& GP, const gtp_par& par);// same as read
		static char		read_vcf(const string& s, const bool to_filter, int& GQ, int& DP, double& GP, const gtp_par& par);												  // same as read
		static void		count(const char g, int& i, int& m, double& n, double& x, const double wt=1); // read g and update #ind, #mis, #chr, #alt
		static void     to_dn(string& s);	// make it alt/ref (diploidy) or alt (haploidy) or alt/ref/ref/.. (N-ploidy)
		static void     to_ref(string& s);	// make it ref|ref
		static void     to_void(string& s); // make it missing
		static void     to_rewr(char g, string& s); // make it missing if is_miss(g), or make it g if s is missing, g is not missing and NoMiss==true
		static void     to_swap(string& s); // VCF only. swap REF/ALT by swapping 0/1 before the first ":", so it is not robust to multi-ALT variants after "vSPLIT -a=0".
		//static void     clean(const char g, string& s) { if (valid(g)) to_void(s); else to_ref(s); } // (valid) should be (!valid)
		static bool		MenErr_auto(const char g_ch, const char g_pa, const char g_ma); // including denovo!
		static bool		MenErr_chrX(const char g_ch, const char g_pa, const char g_ma); // including denovo!
		static bool		MenErr_chrY(const char g_ch, const char g_pa, const char g_ma); // including denovo!
		static bool		MenErr_chrM(const char g_ch, const char g_pa, const char g_ma); // including denovo!
	};

	// The following two functions are over-simplified because it returns one MAf instead of probability of 3 genotypes.
	// However, they have the benefit of ease of use and doesn't inflate type-I-error for association tests require HWE.
	double up_af(const double& af, const vector<double>& penetr); // calculate MAF in cases    given population MAF (af) and a penetrance model (penetr)
	double dn_af(const double& af, const vector<double>& penetr); // calculate MAF in controls given population MAF (af) and a penetrance model (penetr)
	void up_dn_af(const double& af, const vector<double>& penetr, double& piA, double& piU);
	// The following two are more accurate in that it takes 3 probabilities and output another 3, not depending on HWE.
	std::tuple<double,double,double> up_af(const std::tuple<double,double,double>& af, const vector<double>& penetr);
	std::tuple<double,double,double> dn_af(const std::tuple<double,double,double>& af, const vector<double>& penetr);
	// The following works for multiallelic variants or haplotypes, in contrast to biallelic variants above. af & rr are for alleles, not genotypes!
	void up_dn_af(const vector<double>& af, const vector<double>& rr, const double& prevalence, vector<double>& cs, vector<double>& ct, char model);
	// The following works for genotypes, in contrast to haplotypes above. af & rr & pn are for genotypes, not alleles or haplotypes!
	void up_dn_af(const vector<double>& af, const vector<double>& rr, const double& prevalence, vector<double>& cs, vector<double>& ct);
	void up_dn_af(const vector<double>& af, const vector<double>& pn,                           vector<double>& cs, vector<double>& ct);

	// ------------------ linkage analysis ------------------
	
	struct HLOD_result {
		double LOD,alpha,HLOD;
		HLOD_result():LOD(0),alpha(0),HLOD(0){}
		HLOD_result(double l, double a, double h):LOD(l),alpha(a),HLOD(h){}
	};
	
	HLOD_result cal_HLOD(const std::multiset<double>& lodmset);
	
	HLOD_result cal_HLOD(const std::vector<double>& lodvec);
	
	// calculate HLOD conditional on another locus, prbvec is the probability of linkage for the other locus/loci
	// note: L+/L- = LR = 10^LOD ; a1=alpha_when_link_to_other_loci ; a2=alpha_when_NOT_link_to_other_loci
	HLOD_result cal_HLOD_conditional(const std::vector<double>& lodvec, const std::vector<double>& prbvec);
	
	// calculate HLOD with individual alpha for each fam, prbvec is the probability of linkage for the target locus
	HLOD_result cal_HLOD_iAlpha(const std::vector<double>& lodvec, const std::vector<double>& prbvec);

	// calculate summation HLOD (overall HLOD for several loci as a whole, not useful for testing a specific locus)
	HLOD_result cal_HLOD_summation(const std::vector<double>& lodvec, const std::vector< std::vector<double> >& lodoth, const std::vector<double>& aj);
	
	HLOD_result cal_HLOD_summation(const std::vector<double>& lodvec, const std::vector< std::vector<double> >& lodoth, std::vector<double>& result_aj, double incr_aj);

	double LOD_to_P(const double LOD);
	// Province MA. (2001) The Signiﬁcance of Not Finding a Gen. AJHG 69:660-3
	// careful: this assume df=1, i.e. only theta varied to maximize likelihood. If >1 parameters varied (eg, alpha), use the other function.
	
	double P_to_LOD(const double P); // reverse of the above
	
	double LOD_to_P(const double LOD, const unsigned df);
	//http://www.sph.umich.edu/csg/abecasis/LAMP/tour/association.html
	// Convert LOD score to p-value: LOD*2ln(10) = chi_square(df=2)
	
	struct NPL_result {
		double delta, chisq, lod_chisq, pval_chisq, Zmean, pval_Zmean;
		NPL_result():delta(0),chisq(0),lod_chisq(0),pval_chisq(1),Zmean(0),pval_Zmean(1){}
		NPL_result(double d, double c, double z);
	};

	NPL_result cal_NPL(const std::vector<double>& Z_scores, double minDelta, double maxDelta); // return <delta,chisq>;
	NPL_result cal_NPL(const std::multiset<double>& Z_scores, double minDelta, double maxDelta); // return <delta,chisq>
	double kosambi2rf(double k);
	double haldane2rf(double h);
	
	// ------------------ linkage disequilibrium ------------------
	
	// http://en.wikipedia.org/wiki/Linkage_disequilibrium (below allele 2 is 1, 1 is 0)
	// This class calculate allele frequency of SNP-A conditional on SNP-B, ie, P(A=1|B=0/1).
	class LDmodel {
	private:
		bool is_set;
		double p1,p2,q1,q2,Dprime;					// input data
		double D,Dmax,X11,X12,X21,X22;				// internal data
	public:
				LDmodel();
		void	set(double p, double q, double Dp); // p=P(A=1), q=P(B=1); D'=Dp
		void	set(double Dp);						// keep original p, q; D'=Dp
		double	cond_p(int B);						// P(A=1|B=0/1)
	};
	
	template <typename T1,typename T2>
	void cal_LD(T1 x11,T1 x12, T1 x21, T1 x22,T2& Dprime,T2& Rsquare);
	// http://en.wikipedia.org/wiki/Linkage_disequilibrium
	
	template <typename T1,typename T2>
	void cal_LD_from_hap_counts(T1 x11,T1 x12, T1 x21, T1 x22,T2& Dprime,T2& Rsquare);
	
	template <typename T>
	void cubex (const T inputarray[3][3], map<string,double>& result);
	// http://www.oege.org/software/cubex/ modified so that result["Dprime"] & result["Rsquared"] is the most probable solution.
	// When no solution exists, result["Dprime"] & result["Rsquared"] doesn't exist. Access to them will return 0, which are good default values.
	// Therefore, no need to check their existence. This happens when at least one of the SNPs is non-polymorphic.
	// The most probable solution = alpha/beta/gamma that has the lowest xxxx-chisq (min_chisq). But sig chisq means genotype out of HWE!
	// Not all alpha/beta/gamma solution exist: it exist only when xxxxposs=1
	// I also changed result["xxxx-e####"] to a double variable for faster computation.
	// result["Error"]=1 seems never will happen, from my simulations.
	
	template <typename T>
	void cubex (const T inputarray[9], map<string,double>& result);

	// ------------------ end ------------------

}

#endif
