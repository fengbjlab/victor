// body file: vAAA_GLR.cpp

#ifndef PERCH_AAA_GLR
#define PERCH_AAA_GLR
#include <map>

namespace ExtCT
{
	struct ExAC_Var {
		int    chr;
		int    bp;
		std::string ref;
		std::string alt;
		double xiU;
		double niU;
		double xiA;
		double niA;
		double del;
		double vqs;
		bool added;
		bool DoNotUse;
		std::string idx;
		ExAC_Var():chr(0),bp(0),xiU(0),niU(0),xiA(0),niA(0),del(0),vqs(0),added(false),DoNotUse(false) {}
		void clear() { DoNotUse=true; }
	};
	extern std::map<std::string, std::map<std::string,ExAC_Var > > ExAC_dta;	// ExAC_dta[GeneSymb][VarIdx]=data
	extern std::map< int, std::map<int,int> >	ExAC_spl;		// ExAC_spl[chr_numb][bp]=#SamplesCovered
	extern std::map< int, std::map<int,int> >	StudySpl;		// StudySpl[chr_numb][bp]=#SamplesCovered
	extern std::set<std::string>				ExAC_QC;		// QC removed variants in ExAC
	extern std::set<std::string>				studyQC;		// QC removed variants in study
	extern double					ExAC_fem;		// proportion of females in ExAC
	extern double					Case_fem;		// proportion of females in case
	extern double					Ctrl_fem;		// proportion of females in controls
	extern double					ExAC_sub;		// proportion of subset in ExAC
	extern double					ExAC_pCv;		// proportion of samples covered by Nx sequencing
	extern double					ExAC_fqc;		// frequency QC
	extern std::string				ExAC_IAC;		// AC in ExAC. _adj means only include individuals with genotype quality (GQ) >= 20 and depth (DP) >= 10.
	extern std::string				ExAC_IAN;		// AN in ExAC. _adj means only include individuals with genotype quality (GQ) >= 20 and depth (DP) >= 10.
	extern std::vector<std::string>	ExAC_pfx;		// input ref,study prefix. correspond to 5 files: 0.cov 0.qc.log 0.ann.del 1.qc.log 2.cov
	extern std::string				ExAC_pop;		// ExAC population. default _adj means only include ind with genotype quality (GQ) >= 20 and depth (DP) >= 10.
	extern std::string				ExAC_VAG;
	extern std::string				ExAC_log;
	extern boost::iostreams::filtering_ostream	ExAC_out;				// file, set by main if !log_fn.empty()
	extern int						ExAC_pad;
	extern bool						ExAC_cso;
	
	void wr_log(const std::string& GeneSymbol, const std::string& variant, const std::string& message);
	void read_files();
	double covered_samples( std::map< int, std::map<int,int> >& CovDB, const int& chr_num, const int& pos, const std::string& ref, const std::string& alt);
	void clear_added(const std::string& GeneSymb);
};

#endif
