#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_genepi.hpp>
#include <mutex> // requires c++11. g++47 needs it! clang++4.0/g++44 doesn't.
#include "vAAA_par.hpp"
#include "vAAA_GLR.hpp"
#include "victor_par.hpp"

namespace ExtCT
{
	using namespace std;
	
	map<string, map<string,ExAC_Var > > ExAC_dta;	// ExAC_dta[GeneSymb][VarIdx]=data
	map< int, map<int,int> >			ExAC_spl;	// ExAC_spl[chr_numb][bp]=#SamplesCovered
	map< int, map<int,int> >			StudySpl;	// StudySpl[chr_numb][bp]=#SamplesCovered
	set<string>							ExAC_QC;	// QC removed variants in ExAC
	set<string> 						studyQC;
	double								ExAC_fem = 0.4319;		// proportion of females in ExAC
	double								Case_fem = 0;			// proportion of females in cases
	double								Ctrl_fem = 0;			// proportion of females in controls
	double								ExAC_sub = 53105;		// sample size of subset in ExAC
	double								ExAC_pCv = 0.9;			// proportion of samples covered by Nx sequencing
	double								ExAC_fqc = 0.000001;	// frequency QC
	string								ExAC_set = "nonTCGA";	// which subset of ExAC (nonTCGA, all, nonPsych)
	string								ExAC_pop = "Adj";		// population. _adj means only include individuals with genotype quality (GQ) >= 20 and depth (DP) >= 10.
	string								ExAC_IAC = "AC_Adj";	// AC in ExAC. _adj means only include individuals with genotype quality (GQ) >= 20 and depth (DP) >= 10.
	string								ExAC_IAN = "AN_Adj";	// AN in ExAC. _adj means only include individuals with genotype quality (GQ) >= 20 and depth (DP) >= 10.
	vector<string>						ExAC_pfx;				// prefix. correspond to 5 files: 0.cov 0.qc.log 0.ann.del (ExAC files) 1.qc.log (Study QC) 2.cov (Study coverage)
	string								ExAC_VAG;
	string								ExAC_log;
	mutex								ExAC_lmt;
	boost::iostreams::filtering_ostream ExAC_out;				// file, set by main if !log_fn.empty()
	int									ExAC_pad = 10;			// padding
	bool								ExAC_cso = false;		// use case from study only

	#define ExAC_wr(msg) wr_log(GeneSymbol,v.idx,msg);
	
	void wr_log(const string& GeneSymbol, const string& variant, const string& message)
	{
		if (!ExAC_log.empty())
		{
			ExAC_lmt.lock();
			ExAC_out<<GeneSymbol<<DLMTR<<variant<<DLMTR<<message<<endl;
			ExAC_lmt.unlock();
		}
	}

	void read_files()
	{
		// check errors
		if (ExAC_pfx.empty()) return;
		if (ExAC_pfx.size()!=3) exit_error("--xct-pfx reqires 3 strings.");

		// prepare parameters
		// to get sample size:
		//   AN=`df.cut -f 8 ExAC.r1.sites.vep.vcf.gz | head -10000 | df.replace -f 1 --strg 'value[tf[],;,AN_SAS=]' | df.grep --max 1 --first`
		//   echo $AN/2 | bc
		if 		(str_has(ExAC_pfx[0],"nonTCGA")||str_has(ExAC_pfx[0],"noTCGA")||str_has(ExAC_pfx[0],"non-TCGA")||str_has(ExAC_pfx[0],"non_TCGA")) ExAC_set="nonTCGA";
		else if (str_has(ExAC_pfx[0],"nonPsych")||str_has(ExAC_pfx[0],"noPsych")||str_has(ExAC_pfx[0],"non-Psych")||str_has(ExAC_pfx[0],"non_Psych")) ExAC_set="nonPsych";
		else if (str_has(ExAC_pfx[0],"all")||str_has(ExAC_pfx[0],"whole")||str_has(ExAC_pfx[0],"entire")) ExAC_set="all";
		if (ExAC_set == "nonTCGA") // 53105 samples
		{
			if 		(ExAC_pop=="Adj")	 {	ExAC_sub = 53105;	ExAC_fem = 0.431917898502966; }
			else if (ExAC_pop=="AFR")	 {	ExAC_sub = 4533;	ExAC_fem = 0.431917898502966; }	// 4533 samples
			else if (ExAC_pop=="AMR")	 {	ExAC_sub = 5608;	ExAC_fem = 0.431917898502966; }	// 5608 samles
			else if (ExAC_pop=="EAS")	 {	ExAC_sub = 3933;	ExAC_fem = 0.431917898502966; }	// 3933 samles
			else if (ExAC_pop=="FIN")	 {	ExAC_sub = 3307;	ExAC_fem = 0.431917898502966; }	// 3307 samles
			else if (ExAC_pop=="NFE")	 {	ExAC_sub = 27173;	ExAC_fem = 0.431917898502966; } // 27173 samles
			else if (ExAC_pop=="OTH")	 {	ExAC_sub = 347;		ExAC_fem = 0.431917898502966; }	// 347 samles
			else if (ExAC_pop=="SAS")	 {	ExAC_sub = 8204;	ExAC_fem = 0.431917898502966; }	// 8204 samles
			else if (ExAC_pop=="MALE")	 {	ExAC_sub = 30168;	ExAC_fem = 0; }	// 30168 samples
			else if (ExAC_pop=="FEMALE") {	ExAC_sub = 22937;	ExAC_fem = 1; }	// 22937 samples
			else exit_error("wrong ExAC population "+ExAC_pop);
		}
		else if (ExAC_set == "all") // 60702 samples
		{
			if 		(ExAC_pop=="Adj")	 {	ExAC_sub = 60702;	ExAC_fem = 0.445767849494251; }
			else if (ExAC_pop=="AFR")	 {	ExAC_sub = 5201;	ExAC_fem = 0.445767849494251; }	// 5201 samples
			else if (ExAC_pop=="AMR")	 {	ExAC_sub = 5789;	ExAC_fem = 0.445767849494251; }	// 5789 samles
			else if (ExAC_pop=="EAS")	 {	ExAC_sub = 4327;	ExAC_fem = 0.445767849494251; }	// 4327 samles
			else if (ExAC_pop=="FIN")	 {	ExAC_sub = 3307;	ExAC_fem = 0.445767849494251; }	// 3307 samles
			else if (ExAC_pop=="NFE")	 {	ExAC_sub = 33369;	ExAC_fem = 0.445767849494251; } // 33369 samles
			else if (ExAC_pop=="OTH")	 {	ExAC_sub = 454;		ExAC_fem = 0.445767849494251; }	// 454 samles
			else if (ExAC_pop=="SAS")	 {	ExAC_sub = 8256;	ExAC_fem = 0.445767849494251; }	// 8256 samles
			else if (ExAC_pop=="MALE")	 {	ExAC_sub = 33644;	ExAC_fem = 0; }	// 33644 samples
			else if (ExAC_pop=="FEMALE") {	ExAC_sub = 27059;	ExAC_fem = 1; }	// 27059 samples
			else exit_error("wrong ExAC population "+ExAC_pop);
		}
		else exit_error("cannot infer ExAC subset from prefix.");
		ExAC_IAC="AC_"+ExAC_pop;
		ExAC_IAN="AN_"+ExAC_pop;

		if (!ExAC_log.empty()) { if (!openOutFile(ExAC_out, ExAC_log)) exit_error("cannot create file "+ExAC_log); }
		
		// read QC log files
		for (Rows_in_File(in,perch::find_file(ExAC_pfx[0]+".qc.log"),4))
		{
			ExAC_QC.insert(in[0]+'_'+in[1]+'_'+in[2]+'_'+in[3]);
		}
		
		for (Rows_in_File(in,ExAC_pfx[1]+".qc.log",4))
		{
			studyQC.insert(in[0]+'_'+in[1]+'_'+in[2]+'_'+in[3]);
		}
		
		// read coverage files
		{
			tfile_format cov_fmt;
			cov_fmt.set_option(SKIP_NOTES,false);
			double ExAC_DP_CUT=0;
			double ExAC_PC_CUT=0;
			double StudyDP_CUT=0;
			double StudyPC_CUT=0;

			for (Rows_in_File(in,perch::find_file(ExAC_pfx[0]+".cov"),&cov_fmt)) // ExAC coverage
			{
				if 		(str_startsw(in[0],"#DP_CUT=")) { if (!read_val_noNaN(in[0].substr(8),ExAC_DP_CUT)) exit_error("Cannot read meta data #DP_CUT= from "+ExAC_pfx[0]+".cov"); in.clear_nf(); continue; }
				else if	(str_startsw(in[0],"#PC_CUT=")) { if (!read_val_noNaN(in[0].substr(8),ExAC_PC_CUT)) exit_error("Cannot read meta data #PC_CUT= from "+ExAC_pfx[0]+".cov"); in.clear_nf(); continue; }
				if (in.NumFields()<3) exit_error("not enough number of fields in "+ExAC_pfx[0]+".cov");
				int chr_num = genepi::read_chr_num(in[0]); if (chr_num<1) continue;
				int bp  = 0; if (!read_val_gt(in[1],bp,0)) continue;
				double spl = 0; if (!read_val_ge_le(in[2],spl,0.0,1.0)) exit_error("the 3rd column of a coverage file should be the proportion of samples (>=0 and <=1) covered");
				if (spl<ExAC_pCv) continue;
				ExAC_spl[chr_num][bp]=spl;
			}
			
			for (Rows_in_File(in,ExAC_pfx[2]+".cov",&cov_fmt)) // Study coverage
			{
				if 		(str_startsw(in[0],"#DP_CUT=")) { if (!read_val_noNaN(in[0].substr(8),StudyDP_CUT)) exit_error("Cannot read meta data #DP_CUT= from "+ExAC_pfx[2]+".cov"); in.clear_nf(); continue; }
				else if	(str_startsw(in[0],"#PC_CUT=")) { if (!read_val_noNaN(in[0].substr(8),StudyPC_CUT)) exit_error("Cannot read meta data #PC_CUT= from "+ExAC_pfx[2]+".cov"); in.clear_nf(); continue; }
				if (in.NumFields()<3) exit_error("not enough number of fields in "+ExAC_pfx[2]+".cov");
				int chr_num = genepi::read_chr_num(in[0]); if (chr_num<1) continue;
				int bp  = 0; if (!read_val_gt(in[1],bp,0)) continue;
				double spl = 0; if (!read_val_ge_le(in[2],spl,0.0,1.0)) exit_error("the 3rd column of a coverage file should be the proportion of samples (>=0 and <=1) covered");
				if (spl<ExAC_pCv) continue;
				StudySpl[chr_num][bp]=spl;
			}
			
			if (!ExAC_DP_CUT) exit_error("No #DP_CUT meta data in "+ExAC_pfx[0]+".cov");
			if (!ExAC_PC_CUT) exit_error("No #PC_CUT meta data in "+ExAC_pfx[0]+".cov");
			if (!StudyDP_CUT) exit_error("No #DP_CUT meta data in "+ExAC_pfx[2]+".cov");
			if (!StudyPC_CUT) exit_error("No #PC_CUT meta data in "+ExAC_pfx[2]+".cov");
			if (ExAC_DP_CUT!=StudyDP_CUT) exit_error("The coverage files for the Study and ExAC were prepared with different DP_CUT parameters.");
			if (ExAC_PC_CUT>ExAC_pCv) exit_error("The coverage files for ExAC were prepared with a PC_CUT parameters higher than the --xct-cov-pc argument.");
			if (StudyPC_CUT>ExAC_pCv) exit_error("The coverage files for the Study were prepared with a PC_CUT parameters higher than the --xct-cov-pc argument.");
		}
		
		// read .ann.del files
		int count0=elog.get_token("variants from external controls added");
		int count1=elog.get_token("variants from external controls removed due to overlapping with QC filtered variants in the study");
		//int count2=elog.get_token("variants from external controls removed due to VQSLOD");
		//int count3=elog.get_token("variants from external controls removed due to hard filters");
		int count4=elog.get_token("variants from external controls removed due to FILTER");
		int count5=elog.get_token("variants from external controls removed due to wrong $CHROM");
		int count6=elog.get_token("variants from external controls removed due to wrong POS");
		int count7=elog.get_token("variants from external controls removed due to chr regions");
		int count8=elog.get_token("variants from external controls removed due to non-SNV");
		int count9=elog.get_token("variants from external controls removed due to non-LoF");
		int count10=elog.get_token("variants from external controls removed due to MaxAF");
		int count11=elog.get_token("variants from external controls removed due to BayesDel");
		int count12=elog.get_token("variants from external controls removed due to MissingRate");
		int count13=elog.get_token("variants from external controls removed due to coverage");
		int count14=elog.get_token("variants from external controls removed due to no AC/AN");
		field_numbers ExAC_info(false,true);
		field_numbers ExAC_symb(false,true);
		field_numbers ExAC_MxAF(false,true);
		field_numbers ExAC_BDel(false,true);
		field_numbers ExAC_func(false,true);
		field_numbers ExAC_fltr(false,true);
		tfile_format ExAC_fmt;
		ExAC_fmt.set_delimiters("\t");
		ExAC_fmt.set_option(SKIP_NOTES,false);
		for (Rows_in_File(in,perch::find_file(ExAC_pfx[0]+".ann.del"),&ExAC_fmt))
		{
			if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
			{
				if (str_has(in[0],"##VictorCommandLine=<ID=vAnnGene,")) { ExAC_VAG=in[0].substr(52); ExAC_VAG.pop_back(); }
				in.clear_nf();
				continue;
			}
			if (exist_any(perch::h_col1, in.contents()))
			{
				for (int i=0;i<in.NumFields();++i)
				{
					if (in[i]=="FILTER")				ExAC_fltr.push_back(i+1);
					if (in[i]=="INFO")					ExAC_info.push_back(i+1);
					if (in[i]==perch::h_MxAF)			ExAC_MxAF.push_back(i+1);
					if (in[i]==perch::h_symb)			ExAC_symb.push_back(i+1);
					if (in[i]==perch::h_func)			ExAC_func.push_back(i+1);
					if (str_startsw(in[i],"BayesDel"))	ExAC_BDel.push_back(i+1);
				}
				if ( ExAC_info.no_input()) exit_error("The INFO column is missing in"+ExAC_pfx[0]+".ann.del");
				if ( ExAC_symb.no_input()) exit_error("The Gene Symbol column is missing in "+ExAC_pfx[0]+".ann.del");
				continue;
			}
			if ( ExAC_info.no_input()) exit_error("The header row is missing in"+ExAC_pfx[0]+".ann.del");
			vector<string> INFO; boost::split(INFO,in[ExAC_info[0]],boost::is_any_of(";"));
			string& ref = in[3];
			string& alt = in[4];
			bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
			bool is_lof = ( ExAC_func.no_input() ? false : str_has(in[ExAC_func[0]],"LoF")||str_has(in[ExAC_func[0]],"NMD") );
			bool is_vks = ( ExAC_func.no_input() ? false : str_has(in[ExAC_func[0]],"knClinSig=1") );
			bool is_DNg = ( ExAC_func.no_input() ? false : perch::is_DomNeg(in[ExAC_func[0]]) );
			bool is_cds = ( ExAC_func.no_input() ? false : perch::is_coding(in[ExAC_func[0]]) );
			string GeneSymbol;
			if (!ExAC_symb.no_input()) GeneSymbol=in[ExAC_symb[0]];
			if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
			if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
			ExAC_Var v;
			
			// skip by exclusion
			v.idx = in[0]+'_'+in[1]+'_'+in[3]+'_'+in[4];
			if (exist_element(studyQC,v.idx)) { elog.add(count1); ExAC_wr("SKIPPED=ExAC_var_removed_by_study_QC"); continue; }
			// skip by VQSLOD or hard filter
			v.vqs = get_value(INFO,"VQSLOD");
			/*if (std::isnan(VQSLOD))
			{
				if (perch::VQSnan)	{ ExAC_QC.insert(v.idx); elog.add(count2); ExAC_wr("SKIPPED=ExAC_var_QC_VQSLOD"); continue; }
				if (perch::HFnoVQ || perch::hardft)
				{
					if (perch::FiltQD) { if (get_value(INFO,"QD")				<perch::FiltQD) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltMQ) { if (get_value(INFO,"MQ")				<perch::FiltMQ) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltFS) { if (get_value(INFO,"FS")				>perch::FiltFS) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltHS) { if (get_value(INFO,"HaplotypeScore")	>perch::FiltHS) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltMR) { if (get_value(INFO,"MQRankSum")		<perch::FiltMR) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltRP) { if (get_value(INFO,"ReadPosRankSum")	<perch::FiltRP) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
				}
			}
			else
			{
				if (is_snv) { if (v.vqs<perch::VQSsnv) { ExAC_QC.insert(v.idx); elog.add(count2); ExAC_wr("SKIPPED=ExAC_var_QC_VQSLOD"); continue; } }
				else		{ if (v.vqs<perch::VQSidl) { ExAC_QC.insert(v.idx); elog.add(count2); ExAC_wr("SKIPPED=ExAC_var_QC_VQSLOD"); continue; } }
				if (perch::hardft)
				{
					if (perch::FiltQD) { if (get_value(INFO,"QD")				<perch::FiltQD) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltMQ) { if (get_value(INFO,"MQ")				<perch::FiltMQ) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltFS) { if (get_value(INFO,"FS")				>perch::FiltFS) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltHS) { if (get_value(INFO,"HaplotypeScore")	>perch::FiltHS) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltMR) { if (get_value(INFO,"MQRankSum")		<perch::FiltMR) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
					if (perch::FiltRP) { if (get_value(INFO,"ReadPosRankSum")	<perch::FiltRP) { ExAC_QC.insert(v.idx); elog.add(count3); ExAC_wr("SKIPPED=ExAC_var_HardFilter"); continue; } }
				}
			} //*/
			if (std::isnan(v.vqs)) v.vqs=0; // it's important for GLR, although it should not happen
			// skip by FILTER
			if (!ExAC_fltr.no_input()) if (!perch::filflt.empty() && !exist_element(perch::filflt,in[ExAC_fltr[0]])) { ExAC_QC.insert(v.idx); elog.add(count4); ExAC_wr("SKIPPED=ExAC_var_QC_FILTER"); continue; }
			// skip by chr,pos
			int chr_num = genepi::read_chr_num(in[0]); if (chr_num<1) { ExAC_QC.insert(v.idx); elog.add(count5); ExAC_wr("SKIPPED=ExAC_var_Chr_wrong"); continue; }
			int bp=-1; if (!read_val_gt(in[1],bp,0)) { ExAC_QC.insert(v.idx); elog.add(count6); ExAC_wr("SKIPPED=ExAC_var_POS_wrong"); continue; }
			if (!perch::within_covered_region(chr_num,bp)) { ExAC_QC.insert(v.idx); elog.add(count7); ExAC_wr("SKIPPED=ExAC_var_out_of_region"); continue; }
			// skip by variant type
			if (AAA::SNonly && !is_snv)									{ ExAC_QC.insert(v.idx); elog.add(count8); ExAC_wr("SKIPPED=ExAC_var_Not_SNV"); continue; }
			if (perch::LFonly && is_lof && perch::is_LoFtol(GeneSymbol)){ ExAC_QC.insert(v.idx); elog.add(count9); ExAC_wr("SKIPPED=ExAC_var_LoF_Tol"); continue; }
			if (perch::LFonly && !is_lof)								{ ExAC_QC.insert(v.idx); elog.add(count9); ExAC_wr("SKIPPED=ExAC_var_Not_LoF"); continue; }
			if (perch::CDonly && !is_cds)								{ ExAC_QC.insert(v.idx); elog.add(count9); ExAC_wr("SKIPPED=ExAC_var_Not_CDS"); continue; }
			if (perch::DomNeg && !is_DNg)								{ ExAC_QC.insert(v.idx); elog.add(count9); ExAC_wr("SKIPPED=ExAC_var_Not_DomNeg"); continue; }
			// skip by MaxAF
			if (perch::filXAF)
			{
				double MaxAF = std::numeric_limits<double>::signaling_NaN();
				if (!ExAC_MxAF.no_input())	read_val(in[ExAC_MxAF[0]],MaxAF);
				else						MaxAF=get_value(INFO,perch::h_MxAF);
				if (std::isnan(MaxAF)) MaxAF=0;
				double f = (MaxAF>0.5 ? 1-MaxAF : MaxAF);
				if (!perch::rf_XAF && f> perch::filXAF) { ExAC_QC.insert(v.idx); elog.add(count10); ExAC_wr("SKIPPED=ExAC_var_filter_MaxAF"); continue; }
				if ( perch::rf_XAF && f<=perch::filXAF) { ExAC_QC.insert(v.idx); elog.add(count10); ExAC_wr("SKIPPED=ExAC_var_filter_MaxAF"); continue; }
			}
			// skip by BayesDel
			v.del = std::numeric_limits<double>::signaling_NaN();
			if (!ExAC_BDel.no_input())	read_val(in[ExAC_BDel[0]],v.del);
			else						v.del=get_value(INFO,"BayesDel");
			if (perch::filter_AnnAF(v.del,GeneSymbol,is_lof,is_vks)) { ExAC_QC.insert(v.idx); elog.add(count11); ExAC_wr("SKIPPED=ExAC_var_filter_BayesDel"); continue; }
			// skip by coverage
			bool covered=true;
			int beg=bp, end=bp;
			if (ref.size()>1) end=end+ref.size()-1;
			if (alt.size()>1) end=end+1;
			for (int j=beg;j<=end;++j)
			{
				if (!exist_element(StudySpl,chr_num) || !exist_element(StudySpl[chr_num],j)) { covered=false; break; }
				if (!exist_element(ExAC_spl,chr_num) || !exist_element(ExAC_spl[chr_num],j)) { covered=false; break; }
			}
			if (!covered) { ExAC_QC.insert(v.idx); elog.add(count13); ExAC_wr("SKIPPED=ExAC_var_not_well_covered"); continue; }
			// skip by availability of xiU niU
			v.xiU=get_value(INFO,ExAC_IAC); if (std::isnan(v.xiU)) { ExAC_QC.insert(v.idx); elog.add(count14); ExAC_wr("SKIPPED=ExAC_var_no_AC"); continue; }
			v.niU=get_value(INFO,ExAC_IAN); if (std::isnan(v.niU)) { ExAC_QC.insert(v.idx); elog.add(count14); ExAC_wr("SKIPPED=ExAC_var_no_AN"); continue; }
			// skip by missing rate
			if (ExAC_pCv)
			{
				double ExAC_ftr=0; // expected number of chromosomes
				if		(genepi::is_autosomal(chr_num)) { ExAC_ftr=2;			}
				else if (genepi::is_chrX	 (chr_num)) { ExAC_ftr=1+ExAC_fem;	}
				else if (genepi::is_chrY	 (chr_num)) { ExAC_ftr=1-ExAC_fem;	}
				else								 	{ ExAC_ftr=1;			}
				double call_rate = v.niU/ExAC_ftr/ExAC_sub;
				if (call_rate<ExAC_pCv) { ExAC_QC.insert(v.idx); elog.add(count12); ExAC_wr("SKIPPED=ExAC_var_QC_call_rate"); continue; }
			}
			// add data
			v.chr=chr_num;
			v.bp=bp;
			v.ref=ref;
			v.alt=alt;
			ExAC_dta[in[ExAC_symb[0]]][v.idx]=v;
			elog.add(count0);
		}
	}
	
	double covered_samples( map< int, map<int,int> >& CovDB, const int& chr_num, const int& pos, const string& ref, const string& alt)
	{
		int b1=pos;
		int b2=pos+ref.size()-1;
		if (ref.size()==1&& alt.size()>1) { ++b2; }
		b1-=ExAC_pad;
		b2+=ExAC_pad;
		double min_covered=std::numeric_limits<double>::max();
		for (int bp=b1; bp<=b2; ++bp)
		{
			if (exist_element(CovDB,chr_num) && exist_element(CovDB[chr_num],bp))
			{
				if (CovDB[chr_num][bp]<min_covered) min_covered=CovDB[chr_num][bp];
			}
			else
			{
				min_covered=0;
				break;
			}
		}
		if (min_covered<ExAC_pCv) min_covered=0;
		return min_covered;
	}
	
	void clear_added(const string& GeneSymb)
	{
		if (!exist_element(ExAC_dta,GeneSymb)) return;
		for (auto &var:ExAC_dta[GeneSymb]) var.second.added=false;
	}
};
