#include "victor_par.hpp"
#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_math.hpp>
#include <cmath>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

/*
 VQSLOD used by ExAC:
 SNV:        -2.632
 InDel:       1.2168
 
 VQSLOD by Kimber
 SNV:   ExAC -2.78, gnomADwes -3.80, PERSPECTIVE-German -5.368, PERSPECTIVE-Dutch -4.626, RGC -0.16, 600Gene+RGC+TCGA::600Gene -3.744, ::TCGA+RGC -3.92
 InDel: ExAC -2.15, gnomADwes -2.68, PERSPECTIVE-German -3.61,  PERSPECTIVE-Dutch -4.208, RGC -3.94, 600Gene+RGC+TCGA::600Gene -3.754, ::TCGA+RGC -3.97

 */
namespace perch
{
	using namespace std;
	typedef genepi::genotype GTP;

	string				TMPDIR = "/tmp/";	// environment variable TMPDIR
	set<string>			rm_ind ;
	set<string>			h_col1 = {"#CHROM","#Chrom","CHROM","Chrom","#CHR","#Chr","#chr","Chr","CHR","chr","Chromosome","#CHROM:POS:REF:ALT","#CHROM:POS:REF:ALT:Func_Gene"};
	string				h_MxAF = "MaxAF";
	string				h_symb = "Func_Gene";
	string				h_func = "Func_Type";
	string				h_fdet = "Func_Detail";
	string				h_SEGb = "SEGlbf";
	string				h_HLRb = "AAAlbf";
	string				h_GLRb = "GLRlbf";
	string				h_mLgP = "AAAmlp";
	string				h_pVal = "AAApvl";
	string				h_Pcrt = "Pcorr";
	string				h_o_r_ = "O.R.";
	string				h_Adet = "AAA_details";
	string				h_Acol = "AAA_collapse";
	string				h_FINb = "FINlbf";
	string				h_Axtr = "AAA_info";
	string				i_func = "vAnnGene";
	string				gnuplot= "gnuplot";
	string				LoFtol = "panel_LoFtol";	// panel of LoF-toerated genes
	string				gsbdel = "BayesDel_GST"; 	// gene-specific BayesDel threshold
	double				VQSsnv = -5.368;	// the minimum among serveral studies
	double				VQSidl = -4.208;	// the minimum among serveral studies
	double				MisCut = 0.01;		// 0.02 is optimized by Mia's H&N project.
	double				FltDel 		= -0.0570105;	// filter by BayesDel. Default to BayesDel_nsfp33a_noAF, will be changed to BayesDel_nsfp33a if meta data has Included_MaxAF_in_BayesDel
	double				FltDel_hsAF = 0.0692655;	// filter by BayesDel_nsfp33a      (none=-INFINITY. nan also work the same but I won't rely on that.) previously -0.0592577, then -0.0570105
	double				FltDel_noAF = -0.0570105;	// filter by BayesDel_nsfp33a_noAF (none=-INFINITY. nan also work the same but I won't rely on that.)
	double				filXAF = 0.01;		// filter by MaxAF
	double				filSAF = 0;			// filter by SplAF, not as good as filFAF (see below)
	double				filPAF = 0;			// filter by PopAF, not ideal becaue it is almost equal to filter by AF in controls
	double				filFAF = 0.05;		// filter by FdrAF, better than filSAF becaue 1) it considers founders and 2) it includes unknown status samples
	int					filFAFminAN=0;		// apply filFAF only when FdrAN>=filFAFminAN. 0 means this criteria not applicable.
	double				SplSzP = 0.000001;	// p-value for the calculation of sample size required for filFAF
	double				FiltQD = 0;
	double				FiltMQ = 0;
	double				FiltFS = 0;
	double				FiltHS = 0;
	double				FiltMR = 0;
	double				FiltRP = 0;
	set<string>			filflt = {".","PASS"};	// must be empty for ExAC
	double				BDELge = std::numeric_limits<double>::signaling_NaN();	// keep var if BDel >= BDELge. (NaN truns it off) recommend the same as FltDel
	double				AgrRAF = 0;				// aggregated risk allele frq. Use it to infer preval. Previously 0.0025 (carrier prob 0.005), which is arbitrary.
	double				preval = 0.025;
	vector<double>		penetr = {0.02,0.1,0.5};// My simulation showed that it is more powerful to use a higher phenocopy and penetrance than the true parameters in simulation.
	bool				Mis_ea = false;
	bool				rf_del = false;			// reverse filter BayseDel
	bool				rf_XAF = false;			// reverse filter MaxAF
	bool				rf_SAF = false;			// reverse filter SplAF
	bool				rf_PAF = false;			// reverse filter PopAF
	bool				rf_FAF = false;			// reverse filter FdrAF
	bool				VarCla = false;			// run mode is variant classification
	bool				CDonly = false;			// analysis restrict to CDS variants
	bool				LFonly = false;			// analysis restrict to LoF variants
	bool				DomNeg = false;			// analysis restrict to dominant negative variants
	bool				no_MHC = false;			// analysis restrict to non-MHC regions
	bool				MHCsol = false;			// analysis restrict to MHC region solely
	bool				VQSnan = false;			// remove variants if VQSLOD annotation is missing
	bool				hardft = false;			// do hard filtering (QD MQ FS)
	bool				HFnoVQ = true;			// do hard filtering (QD MQ FS) only when there's no VQSLOD.
	bool				a_iPop = false;			// adjust for inferred population instead of principle components
	bool				flipaf = false;			// flip affection status (for one-sided test to find protecitve genes)
	bool				_Debug = false;			// more output
	vector<string>		h_afID;					// header of allele frequency IDs
	// below shoud be all lowercase
	set<string>			h_pid = { "pid","fid","ped","fam","ped_id","fam_id","pedigree","family","pedid","famid","pedigreeid","familyid","pedigree_id","family_id" };
	set<string>			h_iid = { "iid","individ","ind","ind_id","indid","individual","individual_id","individualid","person","person_id","personid","subject","subject_id","subjectid" };
	set<string>			h_sid = { "seqid", "seq_id", "sequencingid", "sequencing_id", "sequenceid", "sequence_id", "exomeid", "exome_id", "genomeid", "genome_id", "experimentid", "experiment_id", "wes_id", "wgs_id", "ngs_id" };
	set<string>			h_gid = { "gene","symbol","genesymbol",boost::to_lower_copy(h_symb) };
	
	namespace {
		vector<string>				conf_r;		// remaining options in configure file
		set<double>					afNo12;		// affection status that's not 1 or 2. If not empty, it's a QTL.
		set<string>					rmi_in;
		set<string>					exclLF;		// genes to be excluded in --lof-only analysis
		set<string>					inclFN;
		set<string>					exclFN;
		vector<genepi::ChrRegions>	InclCR;
		vector<genepi::ChrRegions>	ExclCR;
		string						DelInp;		// --filt-del input
		map<string, double>			bgrDel; 	// bgrDel[GeneSymbol]=max_BayesDel
		set<string>					bgrLoF;		// bgrLoF=GeneSymbol that has LoF variants with MaxAF>
		bool						_noWeb = false;	// do not check version
		bool						_Chk_V = false; // check version
		map<string, double>			gs_Del; 		// gs_Del[GeneSymbol]=BayesDel_cutoff
		map<string, double>			gs_Del_wiAF; 	// gs_Del_wiAF[GeneSymbol]=BayesDel_cutoff with allele frequency
		map<string, double>			gs_Del_woAF; 	// gs_Del_woAF[GeneSymbol]=BayesDel_cutoff without allele frequency

		void _add_bed(const set<string>& files, const string& s2, vector<genepi::ChrRegions>& regions)
		{
			int bp; ReadStr(s2,bp);
			for (auto &fn:files)
			{
				if (fn.empty()) exit_error("Region File name cannot be empty");
				regions.push_back(genepi::ChrRegions());
				regions.back().setup(fn,true,bp,true);
			}
		}
		
		string version_number;
		string version_string;
		void check_prog_version()
		{
			lns<<showl<<"Programs version: "<<version_string<<flush_logger;
			if (!_noWeb && !program.no_web())
			{
				string url1="http://www.fengbj-laboratory.org/download/victor_version";
				string url2="https://www.dropbox.com/s/dp2hyxuikacx1zd/victor_version";
				string read_version;
				int chk_net_vers_result = chk_net_vers(url1,url2,version_string,read_version);
				if		(chk_net_vers_result==-1) lns<<showl<<"Programs are outdated. The new version is "+read_version<<flush_logger;
				else if (chk_net_vers_result== 1) lns<<showl<<"Programs are up-to-date."<<flush_logger;
			}
		}
		void check_data_version()
		{
			string this_version;
			for (Rows_in_File(in,DBpath()+"version",1))
			{
				this_version=in[0];
				break;
			}
			lns<<showl<<DBpath()<<" version: "<<this_version<<flush_logger;
			if (!_noWeb && !program.no_web())
			{
				string url1="http://www.fengbj-laboratory.org/download/"+DBname()+"_version";
				string url2="https://www.dropbox.com/s/x3raywmt5zk5wp0/"+DBname()+"_version";
				string read_version;
				int chk_net_vers_result = chk_net_vers(url1,url2,this_version,read_version);
				if		(chk_net_vers_result==-1) lns<<showl<<DBpath()<<" is outdated. The new version is "+read_version<<flush_logger;
				else if (chk_net_vers_result== 1) lns<<showl<<DBpath()<<" is up-to-date."<<flush_logger;
			}
		}
		void check_update_version()
		{
			if (FileExists(DBpath()+"update_version"))
			{
				string this_version;
				for (Rows_in_File(in,DBpath()+"update_version",1))
				{
					this_version=in[0];
					break;
				}
				lns<<showl<<DBpath()<<" LatestUpdates version: "<<this_version<<flush_logger;
				if (!_noWeb && !program.no_web())
				{
					string url1="http://www.fengbj-laboratory.org/download/"+DBname()+"_update";
					string url2="https://www.dropbox.com/s/x3raywmt5zk5wp0/"+DBname()+"_update";
					string read_version;
					int chk_net_vers_result = chk_net_vers(url1,url2,this_version,read_version);
					if		(chk_net_vers_result==-1) lns<<showl<<DBpath()<<" LatestUpdates is outdated. The new version is "+read_version<<flush_logger;
					else if (chk_net_vers_result== 1) lns<<showl<<DBpath()<<" LatestUpdates is up-to-date."<<flush_logger;
				}
			}
			else
			{
				lns<<showl<<DBpath()<<" LatestUpdates has not been installed."<<flush_logger;
			}
		}

		void _check_parameters()
		{
			// database path
			filename_change_home_path(genepi::MainPath);
			if (str_endsw(genepi::MainPath,"/")) genepi::MainPath.pop_back();
			if 		(DirExists(genepi::MainPath)) genepi::MainPath.push_back('/');
			else if (DirExists(program.exe_path()+"data/"+genepi::MainPath)) genepi::MainPath=program.exe_path()+"data/"+genepi::MainPath+"/";
			else exit_error("cannot find data folder "+genepi::MainPath);
			genepi::set_path(genepi::MainPath);
			lns<<showl<<"Set data folder to "<<genepi::MainPath<<flush_logger;
			if (_Chk_V) { check_prog_version(); check_data_version(); check_update_version(); exit(0); }

			// disease model
			if (preval<0 || preval>=1)				exit_error("prevalence must be within (0,1).");
			if (!penetr.empty())
			{
				if (penetr.size()!=3)				exit_error("penetrance should be 3 numbers, eg: 0.05,0.5,0.5");
				if (penetr[0]<=0 || penetr[0]>=1)	exit_error("penetrance must be within (0,1).");
				if (penetr[1]<=0 || penetr[1]>=1)	exit_error("penetrance must be within (0,1).");
				if (penetr[2]<=0 || penetr[2]>=1)	exit_error("penetrance must be within (0,1).");
				if (preval==0 && AgrRAF) preval = (1-AgrRAF)*(1-AgrRAF)*penetr[0] + 2*(1-AgrRAF)*AgrRAF*penetr[1] + AgrRAF*AgrRAF*penetr[2];
				// if (preval>0 && preval<=penetr[0])	exit_error("prevalence should be > pen[0]."); // no more checking between prevalence and penetrance, becaue
				// if (preval>0 && preval>=penetr[2])	exit_error("prevalence should be < pen[2]."); // prv is for general population, pen is for case-control sample
			}
			// if (preval==0)		exit_error("prevalence is not set"); // removed for lazy checking
			// if (penetr.empty())	exit_error("penetrance is not set"); // removed for lazy checking
			
			// variant filters
			set<string> tmp1; for (auto &fn:inclFN) tmp1.insert(find_file(fn)); inclFN=tmp1;
			set<string> tmp2; for (auto &fn:exclFN) tmp2.insert(find_file(fn)); exclFN=tmp2;
			if (!inclFN.empty()) _add_bed(inclFN,"0",InclCR);
			if (!exclFN.empty()) _add_bed(exclFN,"0",ExclCR);
			if (hardft || HFnoVQ)
			{
				if (FiltQD==0) FiltQD=2;
				if (FiltMQ==0) FiltMQ=40;
				if (FiltFS==0) FiltFS=60;
				if (FiltHS==0) FiltHS=13;
				if (FiltMR==0) FiltMR=-12.5;
				if (FiltRP==0) FiltRP=-8;
			}
			if (hardft && HFnoVQ) exit_error("--hard-filter and --HardFiltIfNoVQ cannot be both true.");
			if (HFnoVQ) VQSnan=false;
			if (filFAF)
			{
				if (filFAF>filXAF)
				{
					// MAF<filXAF is the varaint I want to include. I made filFAF>filXAF so that I will not exclude the desired RV by filFAF.
					// The sample size should be large enough so that the probability of seeing >=k variants from n chromosomes given p=filXAF is low, where k/n=filFAF.
					// I should get filFAFminAN=300 if filFAF=0.05 filXAF=0.01 SplSzP=0.000001
					// I should get filFAFminAN=90  if filFAF=0.10 filXAF=0.01 SplSzP=0.000001
					filFAFminAN=0;
					double factor = 1/filFAF;
					for (double k=1;k<10000;++k)
					{
						double n=k*factor;
						if (cdf_binomial_ge(k,filXAF,n)<SplSzP) { filFAFminAN=n; break; }
					}
					if (!filFAFminAN)
						lns<<showw<<"Couldn't find a solution for the sample size requirement for --filt-FdrAF."<<flush_logger; // Please increase the difference between the --filt-FdrAF and --filt-MaxAF cutoff values.
				}
				//else
				//	lns<<showw<<"--filt-FdrAF cutoff is not higher than the --filt-MaxAF cutoff"<<flush_logger;
			}
			if (LFonly && !LoFtol.empty())
			{
				for (Rows_in_File(in,find_file(LoFtol),1)) exclLF.insert(in[0]);
			}
			if (!DelInp.empty())
			{
				if (str_has(DelInp,".ann.del"))
				{
					if (!filXAF)	exit_error("Adaptive BayesDel filter should be used together with a MaxAF filter (--filt-MaxAF)");
					if (rf_XAF)		exit_error("Adaptive BayesDel filter should not be used together with a reversed MaxAF filter (--rev-filt-MaxAF)");
					string filAnn=DelInp;
					// read file
					{
						string			h_dels;				// header of BayesDel defined by user
						field_numbers	FldSym(false,true);	// field numb for GeneSymbol
						field_numbers	FldFun(false,true);	// field numb for FuncConseq
						field_numbers	FldXAF(false,true);	// field numb for MaxAF
						field_numbers	FldInf(false,true);	// field numb for INFO
						int	ColDel = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
						tfile_format det_format;
						det_format.set_delimiters("\t");
						det_format.set_option(SKIP_NOTES,false);
						for (Rows_in_File(in,find_file(filAnn),&det_format))
						{
							if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
							{
								in.clear_nf();
								if ( h_dels.empty() && str_startsw(in[0],"##INFO=<ID=BayesDel")) ColDel=-1;
								if (!h_dels.empty() && str_startsw(in[0],"##INFO=<ID="+h_dels+",")) ColDel=-1;
								continue;
							}
							if (exist_any(h_col1, in.contents()))
							{
								FldSym.clear();
								FldFun.clear();
								FldXAF.clear();
								FldInf.clear();
								for (int i=0;i<in.NumFields();++i)
								{
									if (in[i]==h_symb)	FldSym.push_back(i+1);
									if (in[i]==h_func)	FldFun.push_back(i+1);
									if (in[i]==h_MxAF)	FldXAF.push_back(i+1);
									if (in[i]=="INFO")	FldInf.push_back(i+1);
									if ( h_dels.empty() && str_startsw(in[i],"BayesDel"))	{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for BayesDel"); }
									if (!h_dels.empty() && in[i]==h_dels)					{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for "+h_dels); }
								}
								if (FldXAF.no_input()) ;														else det_format.set_field_nums(FldXAF,"",tfile_format::Exit);
								if (FldInf.no_input()) exit_error("The INFO  column is missing.");				else det_format.set_field_nums(FldInf,"",tfile_format::Exit);
								if (FldSym.no_input()) exit_error("The "+h_symb+" column is missing.");	else det_format.set_field_nums(FldSym,"",tfile_format::Exit);
								if (FldFun.no_input()) exit_error("The "+h_func+" column is missing.");	else det_format.set_field_nums(FldFun,"",tfile_format::Exit);
								if (ColDel==-2) exit_error("BayesDel annotation is missing."); else if (ColDel>=0) det_format.set_field_nums(ColDel+1,"",tfile_format::Exit);
								continue;
							}
							if (FldSym.no_input()) exit_error("header row not read from "+filAnn);
							
							vector<string> INFO;
							if (!FldInf.no_input())
							{
								if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
							}
							
							double MaxAF = std::numeric_limits<double>::signaling_NaN();
							if (!FldXAF.no_input())	read_val(in[FldXAF[0]],MaxAF);
							else					MaxAF=get_value(INFO,h_MxAF);
							if (std::isnan(MaxAF)) MaxAF=0;
							double f = (MaxAF>0.5 ? 1-MaxAF : MaxAF);
							if (f<=filXAF) { continue; }
							
							double BayesDel = std::numeric_limits<double>::signaling_NaN();
							if		(ColDel==-1)	BayesDel = get_value_sw(INFO,"BayesDel");
							else if (ColDel>=0)		read_val(in[ColDel],BayesDel);
							
							string GeneSymbol;
							if (!FldSym.no_input()) GeneSymbol=in[FldSym[0]];
							if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
							if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
							
							if (!std::isnan(BayesDel))
							{
								if (exist_element(bgrDel,GeneSymbol)) { if (BayesDel>bgrDel[GeneSymbol]) bgrDel[GeneSymbol]=BayesDel; }
								else bgrDel[GeneSymbol]=BayesDel;
							}
							
							bool is_lof = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"LoF")||str_has(in[FldFun[0]],"NMD") );
							if (is_lof) bgrLoF.insert(GeneSymbol);
						}
						lns<<showl<<"Read annotated variants from "<<filAnn<<flush_logger;
						
						// check healthy
						Values<double> cutoffs;
						for (auto &x:bgrDel) cutoffs.push_back(x.second);
						lns<<showl<<"BayesDel cutoffs have median="<<cutoffs.get(STAT::MEDIAN)<<" mean="<<cutoffs.get(STAT::MEAN)<<flush_logger;
						
						// for (auto &x:bgrDel) lns<<showl<<" "<<x.first<<' '<<x.second<<flush_logger;
					}
				}
				else
				{
					if (!read_val_noNaN(DelInp,FltDel)) exit_error("--filt-del argument cannot be NaN");
					FltDel_hsAF=FltDel;
					FltDel_noAF=FltDel;
					if (!rf_del && FltDel== INFINITY)	exit_error("--filt-del cutoff cannot be inf");
					if ( rf_del && FltDel==-INFINITY)	exit_error("--rev-filt-del cutoff cannot be -inf");
				}
			}
			if (!gsbdel.empty())
			{
				for (Rows_in_File(in,find_file(gsbdel),3))
				{
					if (exist_any_lower(h_gid,in.contents())) continue;
					
					string& GeneSymbol=in[0];
					if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
					if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
					if (GeneSymbol.empty()) continue;
					
					double BayesDel_woAF = std::numeric_limits<double>::signaling_NaN();
					read_val(in[1],BayesDel_woAF);
					if (std::isnan(BayesDel_woAF)) continue;
					
					double BayesDel_wiAF = std::numeric_limits<double>::signaling_NaN();
					read_val(in[2],BayesDel_wiAF);
					if (std::isnan(BayesDel_wiAF)) continue;
					
					gs_Del_woAF[GeneSymbol]=BayesDel_woAF;
					gs_Del_wiAF[GeneSymbol]=BayesDel_wiAF;
				}
				gs_Del=gs_Del_woAF;
			}

			// sample filters
			if (!rmi_in.empty())
			{
				tfile_format format;
				format.forbid_nf_rpt();
				format.set_field_nums(1,"missing column 1",tfile_format::Continue);
				for (Rows_in_File(in,rmi_in,&format)) rm_ind.insert(in[0]);
			}

		}
	}
	
	void _NeedArg(const vector<string>& in, int how_many)
	{
		if ((int)in.size()<(how_many+1)) exit_error("Missing argument for "+in[0]+" in parameter file.");
	}
	
	template <typename T>
	void _ReadArg(const vector<string>& in, T& value, int ErrCode=-1)
	{
		std::size_t found = in[0].find('=');
		std::string input;
		if (found!=std::string::npos)	{ _NeedArg(in,0); input = in[0].substr(found+1); }
		else							{ _NeedArg(in,1); input = in[1]; }
		ReadStr(input,value,ErrCode);
	}
	
	template <typename Container>
	void _ReadSet(const vector<string>& in, Container& values, std::string d=",:;|\\ \t\n", int ErrCode=-1)
	{
		std::string input;
		std::size_t found = in[0].find('=');
		char b4='\0';
		if (found!=std::string::npos)	{ b4 = in[0][found-1];	input = in[0].substr(found+1); }
		else							{ _NeedArg(in,1);		input = in[1]; }
		if (b4=='+')		{                 ReadSet_add(input,values,true, d,ErrCode); } // for +=
		else				{ values.clear(); ReadSet_add(input,values,false,d,ErrCode); } // for  =
	}

	void read_configure(const string& filename)
	{
		if (!FileExists(filename)) return;
		lns<<showl<<"Read parameters from "<<filename<<flush_logger;
		string opt,arg;
		tfile_format format;
		format.set_option(SUCCESSIVE_DELIMITERS_AS_ONE,true);
		format.set_option(SKIP_BLANKS,true);
		format.set_field_nums(1,"missing column 1",tfile_format::Continue);
		format.forbid_nf_rpt();
		for (Rows_in_File(in,filename,&format))
		{
			if (in.RowNumber()==0) { if (!str_startsw(in[0],"VICTOR_parameters_1.0")) break; else continue; }
			if		(str_startsw(in[0],"--col-1"))			_ReadSet(in.contents(),h_col1);
			else if	(str_startsw(in[0],"--col-one"))		_ReadSet(in.contents(),h_col1);
			else if	(str_startsw(in[0],"--col-func"))		_ReadArg(in.contents(),h_func);
			else if	(str_startsw(in[0],"--col-gene"))		_ReadArg(in.contents(),h_symb);
			else if	(str_startsw(in[0],"--col-fdet"))		_ReadArg(in.contents(),h_fdet);
			else if	(str_startsw(in[0],"--info-func"))		_ReadArg(in.contents(),i_func);
			else if	(str_startsw(in[0],"--filt-vqs-snv"))	_ReadArg(in.contents(),VQSsnv);
			else if	(str_startsw(in[0],"--filt-vqs-indel"))	_ReadArg(in.contents(),VQSidl);
			else if	(str_startsw(in[0],"--filt-miss-rate"))	_ReadArg(in.contents(),MisCut);
			else if (str_startsw(in[0],"--filt-miss-ea"))	_ReadArg(in.contents(),Mis_ea);
			else if	(str_startsw(in[0],"--filt-filter"))	_ReadSet(in.contents(),filflt);
			else if	(str_startsw(in[0],"--gene-wise-del"))	_ReadArg(in.contents(),gsbdel);
			else if (str_startsw(in[0],"--filt-del"))	 {	_ReadArg(in.contents(),DelInp); rf_del=false; }
			else if (str_startsw(in[0],"--rev-filt-del")){	_ReadArg(in.contents(),DelInp); rf_del=true; }
			else if (str_startsw(in[0],"--filt-MaxAF"))  {	_ReadArg(in.contents(),filXAF); rf_XAF=false; }
			else if (str_startsw(in[0],"--filt-SplAF"))  {	_ReadArg(in.contents(),filSAF); rf_SAF=false; }
			else if (str_startsw(in[0],"--filt-PopAF"))  {	_ReadArg(in.contents(),filPAF); rf_PAF=false; }
			else if (str_startsw(in[0],"--filt-FdrAF"))  {	_ReadArg(in.contents(),filFAF); rf_FAF=false; }
			else if (str_startsw(in[0],"--rev-filt-MaxAF")){_ReadArg(in.contents(),filXAF); rf_XAF=true; }
			else if (str_startsw(in[0],"--rev-filt-SplAF")){_ReadArg(in.contents(),filSAF); rf_SAF=true; }
			else if (str_startsw(in[0],"--rev-filt-PopAF")){_ReadArg(in.contents(),filPAF); rf_PAF=true; }
			else if (str_startsw(in[0],"--rev-filt-FdrAF")){_ReadArg(in.contents(),filFAF); rf_FAF=true; }
			else if (str_startsw(in[0],"--SplSzP"))			_ReadArg(in.contents(),SplSzP);
			else if	(str_startsw(in[0],"--filt-QD"))		_ReadArg(in.contents(),FiltQD);
			else if	(str_startsw(in[0],"--filt-MQ"))		_ReadArg(in.contents(),FiltMQ);
			else if	(str_startsw(in[0],"--filt-FS"))		_ReadArg(in.contents(),FiltFS);
			else if	(str_startsw(in[0],"--filt-HS"))		_ReadArg(in.contents(),FiltHS);
			else if	(str_startsw(in[0],"--filt-MR"))		_ReadArg(in.contents(),FiltMR);
			else if	(str_startsw(in[0],"--filt-RP"))		_ReadArg(in.contents(),FiltRP);
			else if	(str_startsw(in[0],"--hard-filter"))	_ReadArg(in.contents(),hardft);
			else if	(str_startsw(in[0],"--HardFiltIfNoVQ"))	_ReadArg(in.contents(),HFnoVQ);
			else if (str_startsw(in[0],"--filt-DP"))		_ReadArg(in.contents(),GTP::DP_cut);
			else if (str_startsw(in[0],"--filt-GQ"))		_ReadArg(in.contents(),GTP::GQ_cut);
			else if (str_startsw(in[0],"--filt-GP"))		_ReadArg(in.contents(),GTP::GP_cut);
			else if (str_startsw(in[0],"--filt-domGP"))		_ReadArg(in.contents(),GTP::GP_dom);
			else if (str_startsw(in[0],"--filt-recGP"))		_ReadArg(in.contents(),GTP::GP_rec);
			else if (str_startsw(in[0],"--no-miss"))		_ReadArg(in.contents(),GTP::NoMiss);
			else if	(str_startsw(in[0],"--araf"))			_ReadArg(in.contents(),AgrRAF);
			else if	(str_startsw(in[0],"--prevalence"))		_ReadArg(in.contents(),preval);
			else if	(str_startsw(in[0],"--penetrance"))		_ReadSet(in.contents(),penetr);
			else if (str_startsw(in[0],"--include"))		_ReadSet(in.contents(),inclFN);
			else if (str_startsw(in[0],"--exclude"))		_ReadSet(in.contents(),exclFN);
			else if (str_startsw(in[0],"--vc"))				_ReadArg(in.contents(),VarCla);
			else if (str_startsw(in[0],"--gnuplot"))		_ReadArg(in.contents(),gnuplot);
			else if (str_startsw(in[0],"--rm-ind"))			_ReadSet(in.contents(),rmi_in);
			else if (str_startsw(in[0],"--flip-aff"))		_ReadArg(in.contents(),flipaf);
			else if (str_startsw(in[0],"--cds-only"))		_ReadArg(in.contents(),CDonly);
			else if (str_startsw(in[0],"--lof-only"))		_ReadArg(in.contents(),LFonly);
			else if (str_startsw(in[0],"--lof-tol"))		_ReadArg(in.contents(),LoFtol);
			else if (str_startsw(in[0],"--lof-no-del"))	{	_ReadArg(in.contents(),LFonly); if (LFonly) { FltDel=-INFINITY; rf_del=false; } }
			else if (str_startsw(in[0],"--dom-neg"))		_ReadArg(in.contents(),DomNeg);
			else if (str_startsw(in[0],"--no-mhc"))			_ReadArg(in.contents(),no_MHC);
			else if (str_startsw(in[0],"--mhc-only"))		_ReadArg(in.contents(),MHCsol);
			else if (str_startsw(in[0],"--filt-vqs-nan"))	_ReadArg(in.contents(),VQSnan);
			else if (str_startsw(in[0],"--no-web"))		{	_ReadArg(in.contents(),_noWeb); program.set_no_web(_noWeb); }
			else if (str_startsw(in[0],"--bdel-cutoff"))	_ReadArg(in.contents(),BDELge);
			else if (str_startsw(in[0],"--quiet"))			_ReadArg(in.contents(),program.quiet);
			else if (str_startsw(in[0],"-q"))				_ReadArg(in.contents(),program.quiet);
			else if (str_startsw(in[0],"--debug"))			_ReadArg(in.contents(),_Debug);
			else if (str_startsw(in[0],"--af"))				_ReadSet(in.contents(),h_afID);
			else if (str_startsw(in[0],"--use-iPop"))		_ReadArg(in.contents(),a_iPop);
			else { conf_r.push_back(in[0]); if (!str_has(in[0],"=")) conf_r.push_back(in[1]); }
			// else exit_error("Unknown option "+in[0]+" in "+filename);
		}
		if (!conf_r.empty())
			lns<<showl<<"Remaining options from "<<filename<<": "<<str_of_container(conf_r,' ',false)<<flush_logger;
	}
	
	string help_text() {
		string result="BASIC PARAMETERES\n";
		result += "--af Ss            Header of allele frequencies {"+str_of_container(h_afID,',',false)+"}\n";
		result += "--col-one Ss       Header of the first column {"+str_of_container(h_col1,',',false)+"}\n";
		result += "--col-func S       Header of the column for functional consequence {"+h_func+"}\n";
		result += "--col-gene S       Header of the column for gene symbol {"+h_symb+"}\n";
		result += "--col-fdet S       Header of the column for functional details {"+h_fdet+"}\n";
		result += "--genome S         Genome assembly (GRCh37, GRCh38) {"+DBname()+"}\n";
		result += "--gdb F            Gene database (refGeneLite, refgeneFull, ensGeneLite, ensGeneFull) {"+genepi::GDB_name+"}\n";
		result += "--info-func S      The INFO sub-field for functional consequence annotation {"+i_func+"}\n";
		result += "--filt-vqs-nan B   Exclude variants if VQSLOD is missing {"+str_YesOrNo(VQSnan)+"}\n";
		result += "--filt-vqs-snv D   Exclude variants if VQSLOD<D (-inf = no filtering) for SNVs {"+ftos(VQSsnv)+"}\n";
		result += "--filt-vqs-indel D Exclude variants if VQSLOD<D (-inf = no filtering) for InDels {"+ftos(VQSidl)+"}\n";
		result += "--filt-miss-rate D Exclude variants if missing rate in cases and controls is > D (1 = no filtering) {"+ftos(MisCut)+"}\n";
		result += "--filt-miss-ea B   Exclude variants by missing rate in cases and controls separately {"+str_YesOrNo(Mis_ea)+"}\n";
		result += "--filt-filter Ss   Exclude variants if FILTER is not one of the Ss {"+str_of_container(filflt,',',false)+"}\n";
		result += "--filt-del D       Exclude variants if BayesDel < D (-inf = no filtering) {"+ftos(FltDel)+"}\n";
		result += "--filt-MaxAF D     Exclude variants if MaxAF > D (0 = no filtering) {"+ftos(filXAF)+"}\n";
		result += "--filt-SplAF D     Exclude variants if SplAF > D (0 = no filtering) {"+ftos(filSAF)+"}\n";
		result += "--filt-PopAF D     Exclude variants if PopAF > D (0 = no filtering) {"+ftos(filPAF)+"}\n";
		result += "--filt-FdrAF D     Exclude variants if FdrAF > D (0 = no filtering) {"+ftos(filFAF)+"}\n";
		result += "--SplSzP D         Exclude variants by FdrAF only when sample size is so large that the probability of excluding a desired variant (MAF < MaxAF cutoff) is <D {"+ftos(SplSzP)+"}\n";
		result += "--filt-QD D        Exclude variants if QD (Qual. By Depth.) < D (0 = no filtering, GATK default = 2) {"+ftos(FiltQD)+"}\n";
		result += "--filt-MQ D        Exclude variants if MQ (Mapping Quality) < D (0 = no filtering, GATK default = 40) {"+ftos(FiltMQ)+"}\n";
		result += "--filt-FS D        Exclude variants if FS (Fisher Strand)   > D (0 = no filtering, GATK default = 60) {"+ftos(FiltFS)+"}\n";
		result += "--filt-HS D        Exclude variants if HaplotypeScore       > D (0 = no filtering, GATK default = 13) {"+ftos(FiltHS)+"}\n";
		result += "--filt-MR D        Exclude variants if MQRankSum            < D (0 = no filtering, GATK default = -12.5) {"+ftos(FiltMR)+"}\n";
		result += "--filt-RP D        Exclude variants if ReadPosRankSum       < D (0 = no filtering, GATK default = -8) {"+ftos(FiltRP)+"}\n";
		result += "--hard-filter B    Hard filter by QD,MQ,FS,HS,MR,RP. Will populate thresholds with GATK defaults. {"+str_YesOrNo(hardft)+"}\n";
		result += "--HardFiltIfNoVQ B Hard filter by QD,MQ,FS,HS,MR,RP if there's no VQSLOD. Will populate thresholds with GATK defaults. {"+str_YesOrNo(HFnoVQ)+"}\n";
		result += "--bdel-cutoff D    keep variant if BayesDel is >= D (NaN turns it off). Will be ued by vQC --rm-ii-MaxAF-lt and vAnnGene --filter. {"+ftos(BDELge)+"}\n";
		result += "--prevalence D     prevalence {"+ftos(preval)+"}\n";
		result += "--penetrance D,D,D penetrance {"+str_of_container(penetr,',',false)+"}\n";
		result += "--include Fs       Restrict analysis to regions     in FILEs. Use 2+ files to define intersection. {"+str_of_container(inclFN,',',false)+"}\n";
		result += "--exclude Fs       Restrict analysis to regions not in FILEs. Use 2+ files to define union. {"+str_of_container(exclFN,',',false)+"}\n";
		result += "--filt-DP I        Exclude genotypes if DP<I (0 = no filter) {"+itos(GTP::DP_cut)+"}\n";
		result += "--filt-GQ I        Exclude genotypes if GQ<I (0 = no filter) {"+itos(GTP::GQ_cut)+"}\n";
		result += "--filt-GP D        Exclude genotypes if GP<D (0 = no filter) {"+ftos(GTP::GP_cut)+"}\n";
		result += "--filt-domGP B     Exclude genotypes by GP with a dominant model {"+str_YesOrNo(GTP::GP_dom)+"}\n";
		result += "--filt-recGP B     Exclude genotypes by GP with a recessive model {"+str_YesOrNo(GTP::GP_rec)+"}\n";
		result += "--no-miss B        Treat missing genotype as homozygous REF {"+str_YesOrNo(GTP::NoMiss)+"}\n";
		result += "--cds-only B       Analysis includes only variants that change the protein sequence {"+str_YesOrNo(CDonly)+"}\n";
		result += "--lof-only B       Analysis includes loss-of-function variants only {"+str_YesOrNo(LFonly)+"}\n";
		result += "--lof-tol F        Analysis excludes LoF-tolerated genes in file F when --lof-only=yes {"+LoFtol+"}\n";
		result += "--dom-neg B        Analysis includes dominant-negative variants only {"+str_YesOrNo(DomNeg)+"}\n";
		result += "--no-mhc B         Analysis excludes MHC regions {"+str_YesOrNo(no_MHC)+"}\n";
		result += "--mhc-only B       Analysis includes the MHC region only {"+str_YesOrNo(MHCsol)+"}\n";
		result += "--gene-wise-del F  Obtain gene-specific BayesDel cutoff from file F {"+gsbdel+"}\n";
		result += "--gnuplot S        Path to gnuplot 4.6.0 or later {"+gnuplot+"}\n";
		result += "--rm-ind Fs        Remove individuals listed in Fs column 1 {"+str_of_container(rmi_in,',',false)+"}\n";
		result += "--flip-af B        Flip affection status {"+str_YesOrNo(flipaf)+"}\n";
		result += "--vc B             Run mode is Variant Classification {"+str_YesOrNo(VarCla)+"}\n";
		return result;
	}
	
	void add_help_text_var()
	{
		genepi::add_help_text_var();
		program.help_text_var("_Default_af",str_of_container(h_afID,',',false));
		program.help_text_var("_Default_col_one",str_of_container(h_col1,',',false));
		program.help_text_var("_Default_col_func",h_func);
		program.help_text_var("_Default_col_gene",h_symb);
		program.help_text_var("_Default_col_fdet",h_fdet);
		program.help_text_var("_Default_info_func",i_func);
		program.help_text_var("_Default_filt_VQS_SNV",ftos(VQSsnv));
		program.help_text_var("_Default_filt_VQS_InDel",ftos(VQSidl));
		program.help_text_var("_Default_filt_miss_t",ftos(MisCut));
		program.help_text_var("_Default_filt_miss_ea",str_YesOrNo(Mis_ea));
		program.help_text_var("_Default_filt_FILTER",str_of_container(filflt,',',false));
		program.help_text_var("_Default_hard_filter",str_YesOrNo(hardft));
		program.help_text_var("_Default_hardft_novq",str_YesOrNo(HFnoVQ));
		program.help_text_var("_Default_araf",ftos(AgrRAF));
		program.help_text_var("_Default_prevalence",ftos(preval));
		program.help_text_var("_Default_penetrance",str_of_container(penetr,',',false));
		program.help_text_var("_Default_include",str_of_container(inclFN,',',false));
		program.help_text_var("_Default_exclude",str_of_container(exclFN,',',false));
		program.help_text_var("_Default_gnuplot",gnuplot);
		program.help_text_var("_Default_rm_ind",str_of_container(rmi_in,',',false));
		program.help_text_var("_Default_flip_aff",str_YesOrNo(flipaf));
		program.help_text_var("_Default_var_cla",str_YesOrNo(VarCla));
		program.help_text_var("_Default_filt_QD",ftos(FiltQD));
		program.help_text_var("_Default_filt_MQ",ftos(FiltMQ));
		program.help_text_var("_Default_filt_FS",ftos(FiltFS));
		program.help_text_var("_Default_filt_HS",ftos(FiltHS));
		program.help_text_var("_Default_filt_MR",ftos(FiltMR));
		program.help_text_var("_Default_filt_RP",ftos(FiltRP));
		program.help_text_var("_Default_cds_only",str_YesOrNo(CDonly));
		program.help_text_var("_Default_lof_only",str_YesOrNo(LFonly));
		program.help_text_var("_Default_lof_tol",LoFtol);
		program.help_text_var("_Default_dom_neg",str_YesOrNo(DomNeg));
		program.help_text_var("_Default_no_mhc",str_YesOrNo(no_MHC));
		program.help_text_var("_Default_mhc_only",str_YesOrNo(MHCsol));
		program.help_text_var("_Default_gene_wise_del",gsbdel);
		program.help_text_var("_Default_vqs_nan",str_YesOrNo(VQSnan));
		program.help_text_var("_Default_FltDel_hsAF",ftos(FltDel_hsAF));
		program.help_text_var("_Default_FltDel_noAF",ftos(FltDel_noAF));
		program.help_text_var("_Default_filt_MaxAF",ftos(filXAF));
		program.help_text_var("_Default_filt_SplAF",ftos(filSAF));
		program.help_text_var("_Default_filt_PopAF",ftos(filPAF));
		program.help_text_var("_Default_filt_FdrAF",ftos(filFAF));
		program.help_text_var("_Default_minAN",itos(filFAFminAN)); // prv _Default_minAN_for_filt_FdrAF
		program.help_text_var("_Default_SplSzP",ftos(SplSzP));
		program.help_text_var("_Default_bdel_cutoff",ftos(BDELge));
		program.help_text_var("_Default_use_iPop",str_YesOrNo(a_iPop));
		program.help_text_var("_Default_filt_DP",itos(GTP::DP_cut));
		program.help_text_var("_Default_filt_GQ",itos(GTP::GQ_cut));
		program.help_text_var("_Default_filt_GP",ftos(GTP::GP_cut));
		program.help_text_var("_Default_filt_domGP",str_YesOrNo(GTP::GP_dom));
		program.help_text_var("_Default_filt_recGP",str_YesOrNo(GTP::GP_rec));
		program.help_text_var("_Default_no_miss",str_YesOrNo(GTP::NoMiss));
	}
	
	void read_arguments() {
		vector<string> srce_opt;
		srce_opt.insert( srce_opt.end(), program.arg().begin(), program.arg().end() );
		srce_opt.insert( srce_opt.end(), conf_r.begin(), conf_r.end() );
		genepi::read_arguments(srce_opt,false);

		vector<int> to_rm(srce_opt.size(),0);
		for (size_t argi=1; argi<srce_opt.size(); ++argi)
		{
			if		(str_startsw(srce_opt[argi],"--col-1"))			{	to_rm[argi]=1; ReadSet(srce_opt,argi,h_col1); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--col-one"))		{	to_rm[argi]=1; ReadSet(srce_opt,argi,h_col1); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--col-func"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,h_func); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--col-gene"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,h_symb); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--col-fdet"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,h_fdet); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--info-func"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,i_func); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-vqs-snv"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,VQSsnv); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-vqs-indel")){	to_rm[argi]=1; ReadArg(srce_opt,argi,VQSidl); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-miss-rate")){	to_rm[argi]=1; ReadArg(srce_opt,argi,MisCut); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-miss-ea"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,Mis_ea); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-filter"))	{	to_rm[argi]=1; ReadSet(srce_opt,argi,filflt); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--hard-filter"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,hardft); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--HardFiltIfNoVQ")){	to_rm[argi]=1; ReadArg(srce_opt,argi,HFnoVQ); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--gene-wise-del"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,gsbdel); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-del"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,DelInp); to_rm[argi]=1; rf_del=false; }
			else if (str_startsw(srce_opt[argi],"--rev-filt-del"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,DelInp); to_rm[argi]=1; rf_del=true; }
			else if (str_startsw(srce_opt[argi],"--filt-MaxAF"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,filXAF); to_rm[argi]=1; rf_XAF=false; }
			else if (str_startsw(srce_opt[argi],"--filt-SplAF"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,filSAF); to_rm[argi]=1; rf_SAF=false; }
			else if (str_startsw(srce_opt[argi],"--filt-PopAF"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,filPAF); to_rm[argi]=1; rf_PAF=false; }
			else if (str_startsw(srce_opt[argi],"--filt-FdrAF"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,filFAF); to_rm[argi]=1; rf_FAF=false; }
			else if (str_startsw(srce_opt[argi],"--rev-filt-MaxAF")){	to_rm[argi]=1; ReadArg(srce_opt,argi,filXAF); to_rm[argi]=1; rf_XAF=true; }
			else if (str_startsw(srce_opt[argi],"--rev-filt-SplAF")){	to_rm[argi]=1; ReadArg(srce_opt,argi,filSAF); to_rm[argi]=1; rf_SAF=true; }
			else if (str_startsw(srce_opt[argi],"--rev-filt-PopAF")){	to_rm[argi]=1; ReadArg(srce_opt,argi,filPAF); to_rm[argi]=1; rf_PAF=true; }
			else if (str_startsw(srce_opt[argi],"--rev-filt-FdrAF")){	to_rm[argi]=1; ReadArg(srce_opt,argi,filFAF); to_rm[argi]=1; rf_FAF=true; }
			else if (str_startsw(srce_opt[argi],"--SplSzP"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,SplSzP); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--araf"))			{	to_rm[argi]=1; ReadArg(srce_opt,argi,AgrRAF); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--prevalence"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,preval); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--penetrance"))	{	to_rm[argi]=1; ReadSet(srce_opt,argi,penetr); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--include"))		{	to_rm[argi]=1; ReadSet(srce_opt,argi,inclFN); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--exclude"))		{	to_rm[argi]=1; ReadSet(srce_opt,argi,exclFN); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--gnuplot"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,gnuplot);to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--rm-ind"))		{	to_rm[argi]=1; ReadSet(srce_opt,argi,rmi_in); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--flip-aff"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,flipaf); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--vc"))			{	to_rm[argi]=1; ReadArg(srce_opt,argi,VarCla); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-DP"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,GTP::DP_cut); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-GQ"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,GTP::GQ_cut); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-GP"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,GTP::GP_cut); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-domGP"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,GTP::GP_dom); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-recGP"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,GTP::GP_rec); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--no-miss"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,GTP::NoMiss); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-QD"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,FiltQD); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-MQ"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,FiltMQ); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-FS"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,FiltFS); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-HS"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,FiltHS); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-MR"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,FiltMR); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-RP"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,FiltRP); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--cds-only"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,CDonly); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--lof-only"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,LFonly); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--lof-tol"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,LoFtol); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--lof-no-del"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,LFonly); to_rm[argi]=1; if (LFonly) { FltDel=-INFINITY; rf_del=false; } }
			else if (str_startsw(srce_opt[argi],"--dom-neg"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,DomNeg); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--no-mhc"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,no_MHC); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--mhc-only"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,MHCsol); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-vqs-nan"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,VQSnan); to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--bdel-cutoff"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,BDELge); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--af"))			{	to_rm[argi]=1; ReadSet(srce_opt,argi,h_afID); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--quiet"))			{	to_rm[argi]=1; ReadArg(srce_opt,argi,program.quiet); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"-q"))				{	to_rm[argi]=1; ReadArg(srce_opt,argi,program.quiet); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--debug"))			{	to_rm[argi]=1; ReadArg(srce_opt,argi,_Debug); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--version"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,_Chk_V); to_rm[argi]=1; }
			else if	(str_startsw(srce_opt[argi],"--use-iPop"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,a_iPop); to_rm[argi]=1; }
		}
		vector<string> tmp;
		for (size_t argi=0; argi<srce_opt.size(); ++argi) { if (!to_rm[argi]) tmp.push_back(srce_opt[argi]); }
		program.arg()=tmp;
	}
	
	void check_arguments()
	{
		_check_parameters();
		add_help_text_var();
		genepi::add_help_text_var();
		program.push_back_help(genepi::help_text());
		program.check_help_request();
	}
	
	bool within_covered_region(int chr_num, int bp)
	{
		if (no_MHC &&  genepi::is_MHC(chr_num,bp)) return false;
		if (MHCsol && !genepi::is_MHC(chr_num,bp)) return false;
		bool incl=true;	 for (auto &i:InclCR) if (!i.contain(chr_num,bp)) incl=false;
		bool excl=false; for (auto &i:ExclCR) if ( i.contain(chr_num,bp)) excl=true;
		return (incl && !excl);
	}
	string DBpath()
	{
		return genepi::MainPath;
	}
	string DBname()
	{
		string res=genepi::MainPath;
		if (str_endsw(res,"/")) res.pop_back();
		res=substr_after_rfind(res,"/");
		return res;
	}
	string find_file(const string& fn) // return new name, does not allow "" or label_stdin()
	{
		string cvt_name=fn;
		if		(str_startsw(genepi::GDB_name,"refGene")) boost::algorithm::replace_all(cvt_name,"@GDB","refGene");
		else if (str_startsw(genepi::GDB_name,"ensGene")) boost::algorithm::replace_all(cvt_name,"@GDB","ensGene");
		else exit_error("unknown gene database "+genepi::GDB_name);
		string new_name;
		if (!find_file_zippedOrNot(cvt_name,new_name))
			if (!find_file_zippedOrNot(DBpath()+cvt_name,new_name))
				exit_error("cannot find the file "+fn+" in the working directory or in "+DBpath());
		return new_name;
	}
	int		read_sex(string& s) // support UnknSex
	{
		boost::to_lower(s);
		if		(s=="2" || s=="f" || s=="female" || s=="woman" || s=="girl") { s="2"; return 2; }
		else if (s=="1" || s=="m" || s=="male"   || s=="man"   || s=="boy")  { s="1"; return 1; }
		else																 { s="0"; return 0; }
	}
	double	read_aff(string& s, const string& default_UnknAff) // support UnknAff
	{
		boost::to_lower(s);
		double res = std::numeric_limits<double>::signaling_NaN();
		if		(s=="case")			res = 2;
		else if	(s=="control")		res = 1;
		else if	(s=="affected")		res = 2;
		else if	(s=="unaffected")	res = 1;
		else if	(s=="aff")			res = 2;
		else if	(s=="unaff")		res = 1;
		else read_val(s,res); // prv: read_val_gt(s,res,0.0); PRO: 0 is unknown, Pedigree file needs "0", vAAA "if (DepVar)" works for both AFF and QTL. CON: 0 or negative not allowed for QTL.
		if (!std::isnan(res) && res!=1 && res!=2) { if (flipaf) exit_error("--flip-aff does not work with QTL"); afNo12.insert(res); }
		if (!std::isnan(res) && flipaf) { if (res==1) res=2; else if (res==2) res=1; }
		if (!std::isnan(res)) s=ftos(res); else s=default_UnknAff; // "0" is compatible with read_val_gt above but not read_val.
		return res;
	}
	bool is_qtl()
	{
		return !afNo12.empty();
	}
	bool is_LoFtol(const string& GeneSymbol)
	{
		return exist_element(exclLF,GeneSymbol);
	}
	void read_SeqID(string& SeqID, bool& is_proband) // IndID in a pedigree file may contain [p] but SeqID in VCF doesn't. Remove it if necessary.
	{
		is_proband=false;
		string NewID = boost::to_lower_copy(SeqID);
		if (NewID=="(p)" || NewID=="[p]" || \
			NewID=="(proband)" || NewID=="[proband]" || \
			NewID=="_(proband)" || NewID=="_[proband]" || \
			NewID=="proband")
		{
			is_proband=true;
			return;
		}
		if (str_endsw(NewID,"(p)"))			{ SeqID.resize(SeqID.size()-3); is_proband=true; }
		if (str_endsw(NewID,"[p]"))			{ SeqID.resize(SeqID.size()-3); is_proband=true; }
		if (str_endsw(NewID,"(proband)"))	{ SeqID.resize(SeqID.size()-9); is_proband=true; }
		if (str_endsw(NewID,"[proband]"))	{ SeqID.resize(SeqID.size()-9); is_proband=true; }
		if (str_endsw(NewID,"_(proband)"))	{ SeqID.resize(SeqID.size()-10);is_proband=true; }
		if (str_endsw(NewID,"_[proband]"))	{ SeqID.resize(SeqID.size()-10);is_proband=true; }
	}
	void get_trios(const string& ped_in, tabular_file& in, map<string,trio>& trios)
	{
		set<string> ped_aff;// all affected individuals in pedigrees, even those removed by QC
		for (Rows_in_File(pi,ped_in,6))
		{
			// skip header
			if (exist_element(h_pid,boost::to_lower_copy(pi[0]))) continue;
			
			// read affection status
			int aff = read_aff(pi[5],"UnknAff");
			if (std::isnan(aff)) continue; // skip samples with missing outcome
			if (aff==2) ped_aff.insert(pi[1]); // make sure no "continue" before this except isnan(aff)
			
			// skip unwanted
			if (pi[1]=="0" || pi[2]=="0" || pi[3]=="0") continue;
			if (exist_element(rm_ind,pi[1])) continue;
			if (exist_element(rm_ind,pi[2])) continue;
			if (exist_element(rm_ind,pi[3])) continue;
			
			// read sex
			int sex = read_sex(pi[4]);
			if (sex==0) continue; // skip sampels with missing sex
			
			trio t;
			t.id[0]=pi[1];
			t.id[1]=pi[2];
			t.id[2]=pi[3];
			t.sex[0]=sex;
			t.sex[1]=1;
			t.sex[2]=2;
			t.aff=aff;
			
			// see whether sequenced
			int sequenced=0;
			for (int i=0;i<in.NumFields();++i)
			{
				if (in[i]==pi[1]) { ++sequenced; t.col[0]=i; }
				if (in[i]==pi[2]) { ++sequenced; t.col[1]=i; }
				if (in[i]==pi[3]) { ++sequenced; t.col[2]=i; }
			}
			if (sequenced!=3) continue; // skip unsequenced
			
			trios[t.id[0]]=t;
		}
		// removed trios that at least one parent is affected
		map<string,trio> bkup(trios);
		trios.clear();
		for (auto &i:bkup)
		{
			trio& t=i.second;
			if (exist_element(ped_aff,t.id[1])||exist_element(ped_aff,t.id[2])) continue;
			trios[t.id[0]]=t;
		}
		lns<<showl<<"got "<<trios.size()<<" sequenced trios from "<<ped_in<<flush_logger;
	}
	bool read_variable(tabular_file& in, const int field, const vector<string>& INFO, const string& header, double& result)
	{
		result=std::numeric_limits<double>::signaling_NaN();
		if		(field==-1)	result = get_value(INFO,header);
		else if (field>=0)	read_val_noNaN(in[field],result);
		if (std::isnan(result)) return false;
		return true;
	}
	bool read_variable(const vector<string>& in, const int field, const vector<string>& INFO, const string& header, double& result)
	{
		result=std::numeric_limits<double>::signaling_NaN();
		if		(field==-1)	result = get_value(INFO,header);
		else if (field>=0)	read_val_noNaN(in[field],result);
		if (std::isnan(result)) return false;
		return true;
	}
	void read_meta(const string& line)
	{
		if (str_startsw(line,"##INFO=<ID=Included_MaxAF_in_BayesDel,")) { FltDel=FltDel_hsAF; gs_Del=gs_Del_wiAF; }
	}
	void read_spl(const string& 					spl_in,		// input Sample File
				  bool								rmUnAf,		// remove samples with missing affection status
				  bool								rmUnCv,		// remove samples with missing covariate
				  bool								no_cov,		// do not read covariate
				  double							unfAff,		// unified affection status: 0 means read from file, other number means all samples have AFF=unfAff
				  set< string >&					h_csID,		// SeqID of cases
				  set< string >&					h_ctID,		// SeqID of controls
				  set< string >&					h_ukID,		// SeqID of unknowns
				  map< string, int >&				SexMap,		// SeqID => gender (1 for male, 2 for female)
				  map< string, double >&			DepMap,		// SeqID => dependent variable
				  map< string, string>& 			PopMap,		// SeqID => population
				  map< string, string >&			StrMap,		// SeqID => strata
				  map< string, vector<double> >&	CovMap,		// SeqID => covariate values. [c] is the same as CovNID.
				  vector<string>&					CovNID)		//          covariate new ID. [c] is the same as CovMap.
	{
		map< string, vector<string> >	CovRaw;		// SeqID => covariate column => value string
		vector<string>					CovNam;		//			covariate column => name
		vector< int >					CovIaN;		// 			covariate column => number of numbers
		vector< map<string,int> >		CovNaN;		// 			covariate column => string => number of observations
		vector< int >					CovIgn;		// 			covariate column => to ignore if 1
		int spl_wrn_1 = elog.get_token("samples omitted due to a missing value for the outcome variable.");
		int spl_wrn_2 = elog.get_token("samples omitted due to a missing value for a covariate.");
		
		int			cols=0;		// number of columns
		const int	covBgn=3;	// covariates start from this column (0-based)
		const int	minCol=3;	// minimum number of columns in a sample file
		const int	Spl_ID=0;	// 0-based column number for SeqID
		const int	SplSex=1;	// 0-based column number for Sex
		const int	SplAff=2;	// 0-based column number for Aff
		int			SplPop=0;	// 0-based column number for Pop (0=not provided, since SeqID is 0)
		bool		HasPop=false; // Sample File contains the population origin information provided by the user
		
		for (Rows_in_File(in,spl_in,0)) // required columns: SeqID Sex Outcome
		{
			if (in.RowNumber()==0)
			{
				cols=in.NumFields();
				if (cols<minCol)															exit_error("Insufficient number of columns in the Sample File "+spl_in);
				boost::to_lower(in[Spl_ID]); if (in[Spl_ID]!="seqid"&&in[Spl_ID]!="sample")	exit_error("The first  column of a Sample File should be SeqID/sample.");
				boost::to_lower(in[SplSex]); if (in[SplSex]!="sex"&&in[SplSex]!="gender")	exit_error("The second column of a Sample File should be sex/gender.");
				for (int c=covBgn;c<cols;++c) { string t=boost::to_lower_copy(in[c]); if (t=="_ipop")				{ 				SplPop=c; break; } }
				for (int c=covBgn;c<cols;++c) { string t=boost::to_lower_copy(in[c]); if (t=="pop"||t=="population"){ HasPop=true;	SplPop=c; break; } }
				for (int c=covBgn;c<cols;++c)	CovNam.push_back(in[c]);
				if (cols>covBgn) { CovIaN.assign(cols-covBgn,0); CovNaN.assign(cols-covBgn,map<string,int>()); CovIgn.assign(cols-covBgn,0); }
				if (HasPop)				{ for (int c=covBgn;c<cols;++c) if (str_startsw(in[c],"_PC_")||in[c]=="_iPop")	CovIgn[c-covBgn]=1; }
				else if (perch::a_iPop)	{ for (int c=covBgn;c<cols;++c) if (str_startsw(in[c],"_PC_"))					CovIgn[c-covBgn]=1; }
				else					{ for (int c=covBgn;c<cols;++c) if (in[c]=="_iPop")								CovIgn[c-covBgn]=1; }
				continue;
			}
			if (cols!=in.NumFields()) exit_error("inconsistent number of columns in "+spl_in);
			if (exist_element(SexMap,in[Spl_ID])) exit_error("duplicated sample "+in[Spl_ID]+" in "+spl_in);
			if (exist_element(perch::rm_ind,in[Spl_ID])) continue; // skip samples to be removed
			
			// read sex
			int		sex = perch::read_sex(in[SplSex]);
			SexMap[in[Spl_ID]] = sex;

			// read aff
			double	dep = std::numeric_limits<double>::signaling_NaN();
			if (unfAff)
			{
				dep = unfAff;
			}
			else
			{
				dep = perch::read_aff(in[SplAff],"UnknAff");
				if (rmUnAf && std::isnan(dep)) { elog.add(spl_wrn_1); continue; } // skip samples with missing outcome
			}

			// read covariates, allow strings, DO ignore samples with missing values
			if (!no_cov && cols>covBgn)
			{
				vector<string> cov;
				bool has_missing=false;	// this individual has at least 1 missing value
				for (int c=covBgn;c<cols;++c)
				{
					if (CovIgn[c-covBgn]) in[c]="1";
					boost::to_upper(in[c]);
					if (in[c]=="" || in[c]=="." || in[c]=="UNKNOWN" || in[c]=="UNKNCOV") has_missing=true;
					else
					{
						cov.push_back(in[c]);
						double val=std::numeric_limits<double>::signaling_NaN();
						if (ReadStr(in[c],val,0) && !std::isnan(val)) CovIaN[c-covBgn]++;
						else CovNaN[c-covBgn][in[c]]++;
					}
				}
				if (has_missing) 	{ elog.add(spl_wrn_2); if (rmUnCv) continue; } // skip samples with missing covariate
				else				CovRaw[in[Spl_ID]]=cov;
			}

			// read pop
			if (!no_cov && SplPop)
			{
				boost::to_upper(in[SplPop]);
				if (in[SplPop]!="" && in[SplPop]!="." && in[SplPop]!="UNKNOWN" && in[SplPop]!="UNKNPOP") PopMap[in[Spl_ID]]=in[SplPop];
			}
			
			// update
			//if (std::isnan(dep)) h_ukID.insert(in[Spl_ID]);
			DepMap[in[Spl_ID]]=dep;
		}
		
		// convert string covariates to numeric covariates
		if (!CovRaw.empty())
		{
			vector<string> effective_covariates;
			for (size_t i=0;i<CovIaN.size();++i)
			{
				if (CovIaN[i]==0 &&  CovNaN[i].empty()) exit_error("covariates in a Sample File cannot be missing for all samples.");
				if (CovIaN[i]!=0 && !CovNaN[i].empty()) exit_error("covariates in a Sample File must be either strings or numbers, but cannot be both.");
				set<string> uniq_strings;
				for (auto &j:CovRaw) uniq_strings.insert(j.second[i]);
				if (uniq_strings.size()==1) continue;
				effective_covariates.push_back(CovNam[i]);
				if (CovIaN[i]!=0)
				{
					if (str_startsw(CovNam[i],"_")) CovNID.push_back("X"+CovNam[i]); else CovNID.push_back(CovNam[i]);
					for (auto &j:CovRaw)
					{
						double val=std::numeric_limits<double>::signaling_NaN();
						ReadStr(j.second[i],val,0);
						CovMap[j.first].push_back(val);
					}
				}
				else
				{
					string 	max_str;
					int 	max_val=-1;
					for (auto &v:CovNaN[i])
						if (v.second>max_val) { max_str=v.first; max_val=v.second; }
					for (auto &v:CovNaN[i])
					{
						if (v.first==max_str) continue;
						if (str_startsw(CovNam[i],"_")) CovNID.push_back("X"+CovNam[i]+"_"+v.first); else CovNID.push_back(CovNam[i]+"_"+v.first);
						for (auto &j:CovRaw)
						{
							if (j.second[i]==v.first) 	CovMap[j.first].push_back(1);
							else 						CovMap[j.first].push_back(0);
						}
					}
				}
			}
			if (!effective_covariates.empty()) lns<<showl<<"Effective covariates include "<<str_of_container(effective_covariates,',')<<flush_logger;
		}
		
		// convert aff=1/2 to aff=0/1 for case-control data
		if (!perch::is_qtl())
		{
			for (auto &i:DepMap)
			{
				if		(i.second==2)		{	h_csID.insert(i.first); i.second-=1; }
				else if (i.second==1)		{	h_ctID.insert(i.first); i.second-=1; }
				else if (std::isnan(i.second))	h_ukID.insert(i.first);
				else 	exit_error("wrong affection status "+itos(i.second));
			}
			lns<<showl<<h_csID.size()<<" cases, "<<h_ctID.size()<<" controls, and "<<h_ukID.size()<<" UnknAff subjects read from "<<spl_in<<flush_logger;
		}
		else
		{
			double min_qtl=std::numeric_limits<double>::max();
			for (auto &i:DepMap) if (i.second<min_qtl) min_qtl=i.second;
			if (min_qtl<=0) { for (auto &i:DepMap) if (!std::isnan(i.second)) i.second=i.second-min_qtl+1; }
			for (auto &i:DepMap) if (!std::isnan(i.second)) h_csID.insert(i.first);
			lns<<showl<<h_csID.size()<<" subjects with a non-missing quantitative trait read from "<<spl_in<<flush_logger;
			if (min_qtl<=0) lns<<showl<<"converted the quantitative trait value (v) by v=v+1-"<<min_qtl<<flush_logger;
		}
		
		// make StrMap
		if (!CovMap.empty())
		{
			for (auto &x:CovMap)
				StrMap[x.first]=str_of_container(x.second,',');
		}
	}
	void clear_fltdel(bool reverse_filter)
	{
		if (reverse_filter) { FltDel=FltDel_hsAF=FltDel_noAF=INFINITY;  rf_del=true;  gs_Del.clear(); bgrDel.clear(); }
		else				{ FltDel=FltDel_hsAF=FltDel_noAF=-INFINITY; rf_del=false; gs_Del.clear(); bgrDel.clear(); }
	}
	void clear_flt_af(bool reverse_filter)
	{
		perch::filXAF=perch::filSAF=perch::filPAF=perch::filFAF=0;
		perch::rf_XAF=perch::rf_SAF=perch::rf_PAF=perch::rf_FAF=reverse_filter;
	}
	bool filter_AnnAF(const double& BayesDel, const string& GeneSymbol, bool is_lof, bool is_vks) // filter out if return true
	{
		if (is_vks) return false;
		if (is_lof)
		{
			// if (exist_element(bgrLoF,GeneSymbol)) return true;
		}
		else
		{
			if (!std::isnan(BayesDel))
			{
				double new_cutoff=FltDel;
				if (exist_element(bgrDel,GeneSymbol))
				{
					if 		(std::isinf(FltDel))		new_cutoff=bgrDel[GeneSymbol];
					else if (bgrDel[GeneSymbol]>FltDel)	new_cutoff=bgrDel[GeneSymbol];
				}
				if (exist_element(gs_Del,GeneSymbol))
				{
					new_cutoff=gs_Del[GeneSymbol];
				}
				if (!std::isinf(new_cutoff))
				{
					if (!rf_del)	{ if (BayesDel< new_cutoff) return true; else return false; }
					else			{ if (BayesDel>=new_cutoff) return true; else return false; }
				}
			}
		}
		return false;
	}
	bool is_coding(const string& FuncConseq)
	{
		if (str_startsw(FuncConseq,"Missense"))		return true;
		if (str_startsw(FuncConseq,"InFrame"))		return true;
		if (str_startsw(FuncConseq,"StopGain"))		return true;
		if (str_startsw(FuncConseq,"StopLoss"))		return true;
		if (str_startsw(FuncConseq,"Frameshift"))	return true;
		if (str_startsw(FuncConseq,"SpliceSite"))	return true;
		if (str_startsw(FuncConseq,"Translation"))	return true;
		// if (str_startsw(FuncConseq,"SpliceAltering")) return true; // Same as LoF, I don't use it because it's predicted.
		return false;
	}
	bool is_DomNeg(const string& FuncConseq)
	{
		if (str_has(FuncConseq,"NMD"))				return false;
		if (str_startsw(FuncConseq,"Missense"))		return true;
		if (str_startsw(FuncConseq,"InFrame"))		return true;
		if (str_startsw(FuncConseq,"StopGain"))		return true;
		if (str_startsw(FuncConseq,"StopLoss"))		return true;
		if (str_startsw(FuncConseq,"Frameshift"))	return true;
		if (str_startsw(FuncConseq,"SpliceSite"))	return true;
		// if (str_startsw(FuncConseq,"SpliceAltering")) return true; // Same as LoF, I don't use it because it's predicted.
		return false;
	}
	bool is_RedPrd(const string& FuncConseq)
	{
		if (str_has(FuncConseq,"NMD"))				return true;
		if (str_startsw(FuncConseq,"Translation"))	return true;
		if (str_startsw(FuncConseq,"Transcription"))return true;
		if (str_startsw(FuncConseq,"StructuralDel"))return true;
		return false;
	}
	void prelude(int argc, char * const argv[])
	{
		version_number = "1.2beta";
		version_string = "Version "+version_number+" Build " + macro_date_to_boost(string(__DATE__));
		program.trademark = "("+version_string+")";
		program.enable_option("--no-web"); // not used, but keep it.
		program.forbid_option("--version");
		string env_TMPDIR = str_env("TMPDIR").value; if (!env_TMPDIR.empty()) TMPDIR=env_TMPDIR; if (!str_endsw(TMPDIR,"/")) TMPDIR.push_back('/');
		string cerr_mutex_file_name = str_env("STDERR_MUTEX").value;
		if (!cerr_mutex_file_name.empty())
		{
			if (!FileExists(cerr_mutex_file_name))
			{
				std::ofstream mutex_file(cerr_mutex_file_name);
				mutex_file.close();
			}
			lns.set_lock(cerr_mutex_file_name);
			elog.set_lock(cerr_mutex_file_name);
		}
		//lns<<showl<<"Program started:";
		//for (int argi=0;argi<argc;++argi) { lns<<' '<<argv[argi]; }
		//lns<<flush_logger;
		boost::filesystem::path full_path(boost::filesystem::current_path());
		string wd=boost::to_lower_copy(full_path.string());
		int found_genome=0;
		string GenomeName;
		if (str_has(wd,"hg19")                      )	{ ++found_genome; GenomeName="hg19"; }
		if (                    str_has(wd,"grch37"))	{ ++found_genome; GenomeName="GRCh37"; }
		if (str_has(wd,"hg38")||str_has(wd,"grch38"))	{ ++found_genome; GenomeName="GRCh38"; }
		if (found_genome==1) genepi::MainPath=GenomeName;
		//else genepi::MainPath="hg19"; // removed this line to force users to name folder by genome or to use the --genome option in command or par.txt
		genepi::GDB_name="refGeneLite";
		genepi::canonic=false;
		genepi::filt_tx=false;
		read_configure(program.exe_path()+"etc/par.txt");
		read_configure("par.txt");
	}
};
