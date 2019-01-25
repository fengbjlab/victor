// body file: victor_par.cpp

#ifndef PERCH_PARAMETERS
#define PERCH_PARAMETERS

#include <set>
#include <string>
#include <vector>
#include <map>
#include <tft/libfbj_file.hpp>

namespace perch
{
	// data structure
	struct trio {
		std::string id[3];		// IID DAD MOM
		std::vector<int> sex;	// IID DAD MOM
		std::vector<int> col;	// 0-based col in VCF
		int aff;				// aff status of child
		trio():sex(3,0),col(3,0),aff(0) {}
	};
	
	// basic parameters
	extern std::string				TMPDIR;	// environment variable TMPDIR
	extern std::set<std::string>	h_col1;	// header of the first column
	extern std::string				h_MxAF;	// header of the "MaxAF" column
	extern std::string				h_symb;	// header of the "Gene Symbol" column
	extern std::string				h_func;	// header of the "Function Type" column
	extern std::string				h_fdet;	// header of the "Functional Details" column
	extern std::string				h_SEGb;	// header of vSEG analysis results as log10 Bayes factors (default)
	extern std::string				h_HLRb;	// header of vAAA analysis results as log10 Bayes factors (default)
	extern std::string				h_GLRb;	// header of vAAA analysis results as log10 Bayes factors (--glr)
	extern std::string				h_mLgP;	// header of vAAA analysis results as -log10p's (--out-mlp)
	extern std::string				h_pVal;	// header of vAAA analysis results as p-values (--out-pv)
	extern std::string				h_Pcrt;	// header of vAAA analysis results as corrected p (--permute)
	extern std::string				h_o_r_;	// header of vAAA analysis results as odds ratios (--out-or)
	extern std::string				h_Adet; // header of vAAA analysis results as details (--detail)
	extern std::string				h_Acol; // header of vAAA analysis results as collapse (--collapse)
	extern std::string				h_FINb;	// header of vFIN analysis results as log10 Bayes factors (default)
	extern std::string				h_Axtr; // header of vAAA vSEG results in addition to the main output
	extern std::string				i_func; // INFO sub-field for functional consequence annotation. vQC looks for Func_Type in it.
	extern std::string				gnuplot;// path to gnuplot
	extern std::string				CADDdb;	// CADD database file
	extern std::string				LoFtol;	// panel of LoF-toerated genes
	extern double					VQSsnv;	// skip variants if VQSLOD < VQSsnv, for SNV
	extern double					VQSidl;	// skip variants if VQSLOD < VQSidl, for InDel
	extern double					MisCut;	// skip variants if missing% > MisCut in either cases or controls, 1 = noQC
	extern double					FltDel;	// exclude variants if BayesDel<FltDel (LJB26 1%=-0.562088, 5%=-0.259269 best_div=-0.1 none=-INFINITY)
	extern double					FltDel_hsAF; // with allele frq
	extern double					FltDel_noAF; // w/o  allele frq
	extern double					filXAF;	// exclude variants if MaxAF>filXAF. 0 means no filtering. Default=0.001 assuming MaxAF populations are unselected (UK10K/G1K).
	extern double					filSAF;	// exclude variants if SplAF>filSAF. 0 means no filtering. Default=0 as it may remove causal var if sample size is small.
	extern double					filPAF;	// exclude variants if PopAF>filPAF. 0 means no filtering. Default=0 as it may remove causal var if sample size is small.
	extern double					filFAF;	// exclude variants if FdrAF>filFAF. 0 means no filtering. Default=0 as it may remove causal var if sample size is small.
	extern int						filFAFminAN;
	extern double					FiltQD;	// exclude variants if    QD<FiltQD. 0 means no filtering. GATK default = 2      QD
	extern double					FiltMQ;	// exclude variants if    MQ<FiltMQ. 0 means no filtering. GATK default = 40     MQ
	extern double					FiltFS;	// exclude variants if    FS>FiltFS. 0 means no filtering. GATK default = 60     FS
	extern double					FiltHS;	// exclude variants if    HS>FiltHS. 0 means no filtering. GATK default = 13     HaplotypeScore
	extern double					FiltMR;	// exclude variants if    MR<FiltMR. 0 means no filtering. GATK default = -12.5  MQRankSum
	extern double					FiltRP; // exclude variants if    RR<FiltRR. 0 means no filtering. GATK default = -8     ReadPosRankSum
	extern std::set<std::string>	filflt; // exclude variants if FILTER is not one of the filflt
	extern std::set<std::string>	rm_ind; // remove individuals
	extern double					BDELge;	// keep var if BDel >= BDELge. (nan trun it off) recommend the same as FltDel
	extern double					AgrRAF;	// aggregated risk allele frequency in a gene. previously 0.01, which is arbitrary
	extern double					preval; // prevalence
	extern std::vector<double>		penetr; // penetrance
	extern bool						Mis_ea; // filter missing rate (MisCut) in each sample set (cs/ct). But if cohort is defined, it always filter in each cohort.
	extern bool						rf_del;	// reverse filter BayseDel
	extern bool						rf_XAF;	// reverse filter MaxAF
	extern bool						rf_SAF;	// reverse filter SplAF
	extern bool						rf_PAF;	// reverse filter PopAF
	extern bool						rf_FAF;	// reverse filter FdrAF
	extern bool						VarCla; // variant classification
	extern bool						CDonly;	// analysis restrict to CDS variants
	extern bool						LFonly;	// analysis restrict to LoF variants
	extern bool						DomNeg;	// analysis restrict to dominant negative variants
	extern bool						no_MHC;	// analysis restrict to non-MHC regions
	extern bool						MHCsol;	// analysis restrict to MHC region only
	extern bool						VQSnan; // remove variants if VQSLOD annotation is missing
	extern bool						hardft;	// do hard filtering (QD MQ FS HaplotypeScore MQRankSum ReadPosRankSum) no matter what             . These two cannot be buth true.
	extern bool						HFnoVQ;	// do hard filtering (QD MQ FS HaplotypeScore MQRankSum ReadPosRankSum) only when there's no VQSLOD. These two cannot be buth true.
	extern bool						_Debug; // more output
	extern bool						a_iPop;	// adjust for inferred population instead of principle components
	extern std::vector<std::string>	h_afID;	// other header names for allele frequencies
	// constances
	extern std::set<std::string>	h_pid;	// header of PedID, all lower case
	extern std::set<std::string>	h_iid;	// header of IndID, all lower case
	extern std::set<std::string>	h_sid;	// header of SeqID, all lower case
	
	// basic functions
	void prelude(int argc, char * const argv[]);
	void read_configure(const std::string& filename);	// will not check_parameters()
	void read_arguments();								// will     check_parameters()
	void check_arguments();
	void add_help_text_var();
	std::string help_text();
	bool		within_covered_region(int chr_num, int bp);
	std::string DBpath();
	std::string DBname();
	std::string	find_file(const std::string& fn);	// return new name
	int			read_sex(std::string& input);	// standardize input and return 0=unknown 1=male 2=female
	double		read_aff(std::string& input);	// standardize input and return >0 or nan, no 0! 2=aff 1=unaff nan=unknown. If there're other values, it's a QTL.
	void		read_SeqID(std::string& SeqID, bool& is_proband); // IndID in a pedigree file may contain [p] but SeqID in VCF doesn't. Remove it if necessary.
	bool		read_variable(tabular_file& in, const int field, const std::vector<std::string>& INFO, const std::string& header, double& result);
	bool		read_variable(const std::vector<std::string>& in, const int field, const std::vector<std::string>& INFO, const std::string& header, double& result);
	void		read_meta(const std::string& line);
	void 		read_spl(const std::string& 							spl_in,		// input Sample File
						 bool											rmUnAf,		// remove samples with missing affection status
						 bool											rmUnCv,		// remove samples with missing covariate
						 bool											no_cov,		// do not read covariate
						 double											unfAff,		// unified AFF: 0 means read from file, other number means all samples have AFF=unfAff
						 std::set< std::string >&						h_csID,		// SeqID of cases
						 std::set< std::string >&						h_ctID,		// SeqID of controls
						 std::set< std::string >&						h_ukID,		// SeqID of unknowns
						 std::map< std::string, int >&					SexMap,		// SeqID => gender (1 for male, 2 for female)
						 std::map< std::string, double >&				DepMap,		// SeqID => dependent variable
						 std::map< std::string, std::string>& 			PopMap,		// SeqID => population
						 std::map< std::string, std::string >&			StrMap,		// SeqID => strata
						 std::map< std::string, std::vector<double> >&	CovMap,		// SeqID => covariate values. [c] is the same as CovNID.
						 std::vector<std::string>&						CovNID);	//          covariate new ID. [c] is the same as CovMap.
	bool		is_qtl(); // whether aff is qtl. Call it after read_aff() for all lines.
	bool 		is_LoFtol(const std::string& GeneSymbol);
	bool 		is_coding(const std::string& FuncConseq); // is in coding regions
	bool 		is_DomNeg(const std::string& FuncConseq); // is dominant negative
	bool 		is_RedPrd(const std::string& FuncConseq); // is reduced production
	void		clear_fltdel(bool reverse_filtering);
	void		clear_flt_af(bool reverse_filter);
	bool 		filter_AnnAF(const double& BayesDel, const std::string& GeneSymbol, bool is_lof, bool is_vks);
	void		get_trios(const std::string& filename, tabular_file& vcf, std::map<std::string,trio>& trios);
};

#endif
