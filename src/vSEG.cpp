// link libs: -lptoc -lnlopt_cxx

/*
 Features:
 1) Only pedigrees with observed variants and 2+ non-missing individuals are analyzed. Previously all pedigrees were analyzed, which has 3 problems:
    A. If the pedigree file is big, it will go out of mlink's limitation
    B. Those pedigrees have a negative LOD score, so the overall LOD could be very negative, although it still has some discrimination power.
    C. The discrimination power is low due to the dilution by these pedigrees
 2) correct for ascertainment, otherwise the LOD also measure association
 3) gene-wise LOD is calcuated as weighted sum of LRs, where variant weight is the posterior probability of deleteriousness (DEL+VQS)
 4) Output either gene-wise or variant-wise LOD scores. For variant-wise, LOD=nan if it has a Mendelian error. In both cases, LOD=0 if not calculated.
 5) Only work for rare variants.

 Removed options:
 --raf STRs        in input, header of ref. allele frq. [EspEA_R]
 --wc FILE         well-covered regions
 --no-wt           Equal to "--wt-del=no --wt-vqs=no"
 --no-del          Do not weight or filter variants by deleteriousness
 */

#include "linkage51.hpp"
#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_genepi.hpp>
#include "victor_par.hpp"

using namespace std;
typedef genepi::genotype GTP;

enum WtMethodType { Weighting_LODs, Weighting_LRs, Weighting_Ls } WtMethod = Weighting_LRs;
vector<double> variant_lnkLH1;
vector<double> variant_lnkLH0;
vector<double> variant_lnkRes;
vector<double> variant_weight;
vector<double> variant_BayesD;
vector<string> variant_filter;
bool single = false; // otherwise gene-based

void write_results(bool wrINFO, field_numbers& FldRes, const string& wrWhat, const string& MisStr, bool is_the_last_group)
{
	string result;
	
	if (!single) // gene-based
	{
		double overall_weight=0;
		double overall_result=0;
		for (size_t	i=0;i<variant_lnkRes.size();++i)
			if (std::isnan(variant_lnkRes[i]) || variant_lnkRes[i]==0) variant_weight[i]=0; else overall_weight+=variant_weight[i];
		if (overall_weight)
		{
			if		(WtMethod == Weighting_LRs)
			{
				MakeFraction(variant_weight);
				for (size_t	i=0;i<variant_lnkRes.size();++i)
					if (!std::isnan(variant_lnkRes[i])) overall_result += variant_weight[i] * pow(10,variant_lnkRes[i]);
				overall_result = log10(overall_result);
			}
			else if	(WtMethod == Weighting_Ls)
			{
				// MakeFraction(variant_weight); // no need; made no difference
				double overall_H1=0;
				double overall_H0=0;
				for (size_t	i=0;i<variant_lnkRes.size();++i)
				{
					overall_H1 += variant_weight[i] * pow(10,variant_lnkLH1[i]);
					overall_H0 += variant_weight[i] * pow(10,variant_lnkLH0[i]);
				}
				overall_result = log10(overall_H1)-log10(overall_H0);
			}
			else if (WtMethod == Weighting_LODs)
			{
				MakeFraction(variant_weight);
				for (size_t	i=0;i<variant_lnkRes.size();++i)
					if (!std::isnan(variant_lnkRes[i])) overall_result += variant_weight[i] * variant_lnkRes[i];
			}
			else exit_error("Wrong WtMethod");
		}
		result = ftos(overall_result,MisStr);
	}

	int i=0;
	if (is_the_last_group)
		for (each_element(program.main_data(),it1),++i)
		{
			string output=result;
			if (std::isnan(variant_lnkRes[i])) output="nan";
			else if (single) output = variant_filter[i].empty() ? ftos(variant_lnkRes[i],MisStr) : "SKIPPED:"+variant_filter[i];
			if (wrINFO) (*it1)[FldRes[0]] += ";" + wrWhat + "=" + output;
			else		(*it1)[FldRes[0]]  = output;
			print_container(*it1,program.outf,DLMTR,true);
		}
	else
		for (each_interval(program.main_data(),it1,it2),++i)
		{
			string output=result;
			if (std::isnan(variant_lnkRes[i])) output="nan";
			else if (single) output = variant_filter[i].empty() ? ftos(variant_lnkRes[i],MisStr) : "SKIPPED:"+variant_filter[i];
			if (wrINFO) (*it1)[FldRes[0]] += ";" + wrWhat + "=" + output;
			else		(*it1)[FldRes[0]]  = output;
			print_container(*it1,program.outf,DLMTR,true);
		}
}

// remember to call ss.clear() before ss.seekg since the stream ss is in the eof state after reading.
void rewind(stringstream& ss)
{
	ss.clear();
	ss.seekg(ss.beg);
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	string			wrWhat = perch::h_SEGb;
	bool			AscertCorr = true;		// Do ascertainment correction (w/o it linkage is also a burden test, overlapping with vAAA)
	bool			SkipPV_mis = false;		// Do not analyze variants with missing values, so that bad-QC variants do not have a low score
	bool			force_defMAF = false;	// force the program to use the defMAF even MaxAF annotation exist
	double			defMAF = 0.00001;		// MaxAF if MaxAF=0. ExAC has 60706 samples, so 1/60706/2=0.00000823641815.
	double			mut_m = 0;				// de novo mutation male
	double			mut_f = 0;				// de novo mutation female
	string			vcf_file;				// input genotype file in VCF format
	string			ped_in;					// input pedigree file before makeped
	string			liFile;					// input liability class file
	bool			SNonly = false;			// analysis restricted to SNVs
	bool			AddInf = false;			// write result to column 8 (INFO)
	bool			Wt_DEL = true;			// weight variants by deleteriousness
	bool			Wt_VQS = true;			// weight variants by variant call quality
	bool			cap_pg = false;			// put a cap on the number of variants
	string			MisStr = "NA";			// missing value string
	string			h_dels;					// header of BayesDel defined by user
	set<string>		h_grID = {perch::h_symb};	// header of group IDs
	bool			do_nothing = false;
	// handle program arguments
	perch::MisCut=1;
	perch::VQSsnv=-INFINITY;
	perch::VQSidl=-INFINITY;
	perch::filflt.clear();
	perch::h_afID.push_back("PopAF");
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"--ped"))				ReadArg(program.arg(),argi,ped_in);
		else if (str_startsw(program.arg()[argi],"--default-maf"))		ReadArg(program.arg(),argi,defMAF);
		else if (str_startsw(program.arg()[argi],"--allele-freq"))	{	ReadArg(program.arg(),argi,defMAF); force_defMAF=true; }
		else if (str_startsw(program.arg()[argi],"--mut-male"))			ReadArg(program.arg(),argi,mut_m);
		else if (str_startsw(program.arg()[argi],"--mut-female"))		ReadArg(program.arg(),argi,mut_f);
		else if (str_startsw(program.arg()[argi],"--liab"))				ReadArg(program.arg(),argi,liFile);
		else if	(str_startsw(program.arg()[argi],"--group"))			ReadSet(program.arg(),argi,h_grID);
		else if (str_startsw(program.arg()[argi],"--wr-info"))			ReadArg(program.arg(),argi,AddInf);
		else if (str_startsw(program.arg()[argi],"--wt-del"))			ReadArg(program.arg(),argi,Wt_DEL);
		else if (str_startsw(program.arg()[argi],"--wt-vqs"))			ReadArg(program.arg(),argi,Wt_VQS);
		else if (str_startsw(program.arg()[argi],"--cap"))				ReadArg(program.arg(),argi,cap_pg);
		else if (str_startsw(program.arg()[argi],"--snv-only"))			ReadArg(program.arg(),argi,SNonly);
		else if (str_startsw(program.arg()[argi],"--single"))			ReadArg(program.arg(),argi,single);
		else if (str_startsw(program.arg()[argi],"--detail"))		{	ReadArg(program.arg(),argi,single); if (single) wrWhat="vSEG_details"; }
		else if (str_startsw(program.arg()[argi],"--do-nothing"))		ReadArg(program.arg(),argi,do_nothing);
		else if (program.arg()[argi]=="--no-asc")		{	AscertCorr=false; }
		else if (program.arg()[argi]=="--no-vqs")		{	Wt_VQS=false; }
		else if (program.arg()[argi]=="--no-del")		{	Wt_DEL=false; perch::clear_fltdel(false); }
		else if (program.arg()[argi]=="--no-wt")		{	Wt_DEL=false; Wt_VQS=false; }
		else if (program.arg()[argi]=="--gene-wise-wtL")	WtMethod = Weighting_Ls;
		else if (program.arg()[argi]=="--gene-wise-wtLR")	WtMethod = Weighting_LRs;
		else if (program.arg()[argi]=="--gene-wise-wtLOD")	WtMethod = Weighting_LODs;
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else if (vcf_file.empty())	vcf_file=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	program.help_text_var("_Default_groups",str_of_container(h_grID,',',false));
	program.help_text_var("_Default_default_maf",ftos(defMAF));
	program.help_text_var("_Default_mut_male",ftos(mut_m));
	program.help_text_var("_Default_mut_female",ftos(mut_f));
	program.help_text_var("_Default_wt_del",str_YesOrNo(Wt_DEL));
	program.help_text_var("_Default_wt_vqs",str_YesOrNo(Wt_VQS));
	program.help_text_var("_Default_liab",liFile);
	program.help_text_var("_Default_cap",str_YesOrNo(cap_pg));
	program.help_text_var("_Default_detail",str_YesOrNo(single));
	program.help_text_var("_Default_snv_only",str_YesOrNo(SNonly));
	program.help_text_var("_Default_single",str_YesOrNo(single));
	program.help_text_var("_Default_do_nothing",str_YesOrNo(do_nothing));
	program.push_back_help(GTP::help_text());
	perch::check_arguments();
	
	if (do_nothing)
	{
		tfile_format format;
		format.set_delimiters("\t");
		format.set_option(SKIP_NOTES,false);
		for (Lines_in_File(in, vcf_file, &format)) program.outf << in[0] << endl;
		return 0;
	}

	// check errors
	if (ped_in.empty())			exit_error("Pedigree File not set.");
	if (perch::penetr.empty())	exit_error("penetrance is not set.");
	if (perch::filXAF>0.1)		exit_error("vSEG only works for low-frequency variants.");
	if (perch::rf_XAF)			exit_error("vSEG only works for low-frequency variants.");
	if (perch::rf_SAF)			exit_error("vSEG only works for low-frequency variants.");
	if (perch::rf_PAF)			exit_error("vSEG only works for low-frequency variants.");
	if (perch::rf_FAF)			exit_error("vSEG only works for low-frequency variants.");
	if (force_defMAF && !perch::VarCla) exit_error("--allele-freq only works with --vc");
	
	linkage::defaultResultType = linkage::LOG10L;
	int	min_AnyCar = 0;			// if min_AnyCar=1 & min_GtpPer=2, must  correct ascertainment
	int	min_GtpPer = 0;			// if min_AnyCar=0 & min_GtpPer=1, don't correct ascertainment
	if (AscertCorr) { min_AnyCar=1; min_GtpPer=2; }
	else			{ min_AnyCar=0; min_GtpPer=1; }
	if (perch::VarCla)
	{
		single=true;
		Wt_DEL=false;
		Wt_VQS=false;
		perch::clear_fltdel(false);
		perch::clear_flt_af(false);
		perch::filflt.clear();
	}
	
	// read liability
	string liability;
	string liab_chrX;
	vector<double> relative_risk = {1};
	if (!liFile.empty())
	{
		{
			openInpFile_or_exit(file,liFile);
			while (next_row(file)) { safeGetline_add(file,liability); liability+='\n'; }
			while (file.peek()!=EOF && is_blank_row(file)) skip_blank_lines(file);
			while (next_row(file)) { safeGetline_add(file,liab_chrX); liab_chrX+='\n'; }
			closefile(file);
			if (liability.empty()) exit_error("The liability class file "+liFile+" does not contain penetrance for autosomal chromosome.");
			if (liab_chrX.empty()) liab_chrX=liability;
			//if (liab_chrX.empty()) exit_error("The liability class file "+liFile+" does not contain penetrance for chromosome X.");
			if (!str_endsw(liability,"\n")) liability += '\n';
			if (!str_endsw(liab_chrX,"\n")) liab_chrX += '\n';
		}
		{
			int num_lia=0;
			openInpFile_or_exit(file,liFile);
			for (;next_row(file);skip_one_line(file))
			{
				if (num_lia)
				{
					double p0,p1,p2;
					file>>p0>>p1>>p2;
					if (p0)			relative_risk.push_back(p2/p0);
					else if (p2)	relative_risk.push_back(INFINITY);
					else			relative_risk.push_back(1);
				}
				else
				{
					file>>num_lia;
					if (num_lia<=0) exit_error("the number of liability class is wrong");
				}
			}
			closefile(file);
			if ((int)(relative_risk.size())!=(num_lia+1))
				exit_error("the number of liability class "+itos(num_lia)+" does not match with the number of lines for penetrance "+itos(relative_risk.size()-1));
		}
	}

	// read pedigree
	vector<string>	L2FID;		// L2FID[line#] = FamilyID
	vector<string>	L2SID;		// L2SID[line#] = SeqID
	vector<int>		L2AFF;		// L2AFF[line#] = Affection status (0/1)
	vector<int>		L2SEX;		// L2SEX[line#] = gender (1/2)
	vector<int>		L2LIA;		// L2LIA[line#] = liability class
	vector<string>	part1;		// part1[line#] = pedigree file column 1-6/7
	vector<int>		part2;		// part2[line#] = 0-based column number in VCF (0 means not sequenced, so never put genotype on column 0)
	map<string,int> fld4ind;	// fld4ind[SID] = 0-based column number in VCF (0 means not sequenced, so never put genotype on column 0)
	map<string,int> fid2prb;	// fid2prb[FID] = line#
	map<string,int> pid_seq;	// pid_seq[FID] = number of sequenced samples
	bool			p_ped = false;	// pedigree file is in PERCH format
	int				p1_sz = 0;		// part1 size: 6 w/o liability, 7 w/ liability, 0 unknown.
	lns<<showl<<"read Pedigree File "<<ped_in<<flush_logger;
	for (Rows_in_File(in, ped_in, 6))
	{
		// header
		if (in.RowNumber()==0 && exist_any_lower(perch::h_sid,in.contents()))
		{
			lns<<showw<<"your Pedigree File "<<ped_in<<" is in the old PERCH format. I can still run for now but please update to the new format."<<flush_logger;
			p_ped=true;
		}
		if (!p1_sz)  p1_sz  = in.NumFields()-(p_ped?1:0);
		else if		(p1_sz != in.NumFields()-(p_ped?1:0)) exit_error("Inconsistant number of columns in "+ped_in);
		if (p1_sz!=6 && p1_sz!=7) exit_error("pedigree file "+ped_in+" format not right ");
		if (exist_any_lower(perch::h_pid,in.contents()) || exist_any_lower(perch::h_sid,in.contents()))
		{
			// check header
			if (p1_sz==7) { boost::to_lower(in[6]); if (!str_startsw(in[6],"liab")) exit_error("pedigree file "+ped_in+" column 7 is not liability class (column head starts with liab)"); }
			// skip header
			continue;
		}
		
		// SeqID
		bool is_proband = false;
		string SeqID = ( p_ped ? in[p1_sz] : in[1] );
		perch::read_SeqID(SeqID,is_proband);
		if (p_ped) in[p1_sz]=SeqID; else in[1]=SeqID;
		if (is_proband)
		{
			if (exist_element(fid2prb,in[0])) exit_error("Error in "+ped_in+": >1 probands in "+in[0]);
			// if (in[5]!="2") exit_error("Error in "+ped_in+": the proband in "+in[0]+" is unaffected");
			fid2prb[in[0]]=part1.size();
		}
		if (SeqID.empty()) SeqID="0";
		if (exist_element(perch::rm_ind,SeqID)) SeqID="0"; // skip samples to be removed

		// sex and aff
		int sex=perch::read_sex(in[4]); if (sex==0) exit_error("Sex in Pedigree File should be either 1 (male) or 2 (female).");
		double aff=perch::read_aff(in[5]);
		L2SEX.push_back(sex);
		if (aff==2) L2AFF.push_back(1);
		else		L2AFF.push_back(0);
		
		// lia
		int lia=0;
		if (p1_sz==7) { if (!read_val_ge(in[6],lia,1)) exit_error("Liability class wrong in Pedigree File"); }
		L2LIA.push_back(lia);
		L2FID.push_back(in[0]);
		L2SID.push_back(SeqID);
		part1.push_back(string()); container_head_to_str(in.contents(),p1_sz,part1.back(),DLMTR);
		if (SeqID!="0")
		{
			fld4ind[SeqID]=0; // 0=unknown, will change when reading VCF (not necessarily all of them)
			++pid_seq[in[0]];
		}
	}
	lns<<showl<<fld4ind.size()<<" individuals in Pedigree File "<<ped_in<<flush_logger;
	int valid_pid=0;
	for (auto &p:pid_seq) if (p.second>1) ++valid_pid;
	if (valid_pid==0) exit_error("No pedigree has more than 1 sequenced samples.");
	lns<<showl<<valid_pid<<" pedigrees have more than 1 sequenced samples."<<flush_logger;
	if (p1_sz==7 && liFile.empty()) exit_error("If the pedigree file contain liability classes, you must use the --liab option to input penetrances.");
	if (perch::VarCla)
	{
		//for (auto &f:L2FID)
		//	if (!exist_element(fid2prb,f)) exit_error("Error in "+ped_in+": no proband in "+f);
	}
	else
	{
		if (!fid2prb.empty())
		{
			lns<<showl<<"Probands are ignored if the analysis is not for variant classificaiton."<<flush_logger;
			fid2prb.clear();
		}
	}
	
	// read VCF and do linkage analysis line by line
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldFmt(false,true);	// field numb for FORMAT
	field_numbers	FldFlt(false,true);	// field numb for FILTER
	field_numbers	FldXAF(false,true);	// field numb for MaxAF
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldSym(false,true);	// field numb for "Gene Symbol"
	field_numbers	FldFun(false,true);	// field numb for "Function Type"
	field_numbers	FldGrp(false,true);	// field numb for group ID
	field_numbers	FldIdx(false,true);	// field numb for chr,start,end,ref,alt (ANNOVAR's 1st 5 columns)
	field_numbers	FldRes(false,true);	// field numb for result
	tfile_format	format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);
	format.set_storage_to(program.main_data());
	program.main_data().keep_all_input();
	string		prev_group;	// the group right before this group, see see whether this is a new group
	set<string> old_groups; // all previous groups, to check whether rows are sorted by groups
	set<string> prev_index;	// chr-pos-id-ref-alt
	bool		header_not_read = true;
	int warning1 = elog.get_token("duplicated variants within a group were skipped.");
	int warning2 = elog.get_token("duplicated groups.");
	size_t	max_pg=0;	// max number of variants for analysis, 0=no_limit
	int	ColDel = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.

	for (Rows_in_File(in, vcf_file, &format))
	{
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			in.clear_nf();
			perch::read_meta(in[0]);
			if (str_has(in[0],"##INFO=<ID="+wrWhat+",")) continue;
			if ( h_dels.empty() && str_startsw(in[0],"##INFO=<ID=BayesDel")) ColDel=-1;
			if (!h_dels.empty() && str_startsw(in[0],"##INFO=<ID="+h_dels+",")) ColDel=-1;
			if (!header_not_read) { program.main_data().clear(); continue; }
			print_container(in.contents(),program.outf,' ',true);
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			if (AddInf) program.outf<<"##INFO=<ID="<<wrWhat<<",Number=1,Type=Float,Description=\"vSEG co-segregation analysis result.\">"<<endl;
			if (!header_not_read) { program.main_data().clear(); continue; }
			lns<<showl<<"BayesDel filter cutoff value is "<<perch::FltDel<<flush_logger;
			header_not_read=false;
			format.clear_field_nums();
			FldInf.clear();
			FldFmt.clear();
			FldFlt.clear();
			FldXAF.clear();
			FldChr.clear();
			FldPos.clear();
			FldRef.clear();
			FldAlt.clear();
			FldSym.clear();
			FldFun.clear();
			FldGrp.clear();
			FldIdx.clear();
			FldRes.clear();
			for (int i=0;i<in.NumFields();++i)
			{
				if (exist_element(h_grID,in[i]))	FldGrp.push_back(i+1);
				if (in[i]==perch::h_MxAF)			FldXAF.push_back(i+1);
				if (in[i]=="INFO")					FldInf.push_back(i+1);
				if (in[i]=="FORMAT")				FldFmt.push_back(i+1);
				if (in[i]=="FILTER")				FldFlt.push_back(i+1);
				if (in[i]==perch::h_symb)			FldSym.push_back(i+1);
				if (in[i]==perch::h_func)			FldFun.push_back(i+1);
				if ( h_dels.empty() && str_startsw(in[i],"BayesDel"))	{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for BayesDel"); }
				if (!h_dels.empty() && in[i]==h_dels)					{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for "+h_dels); }
				if (in[i]=="Chr"   && FldChr.no_input())	FldChr.push_back(i+1); if (in[i]=="#CHROM")		{	FldChr.clear(); FldChr.push_back(i+1); }
				if (in[i]=="Start" && FldPos.no_input())	FldPos.push_back(i+1); if (in[i]=="POS")		{	FldPos.clear(); FldPos.push_back(i+1); }
				if (in[i]=="Ref"   && FldRef.no_input())	FldRef.push_back(i+1); if (in[i]=="REF")		{	FldRef.clear(); FldRef.push_back(i+1); }
				if (in[i]=="Alt"   && FldAlt.no_input())	FldAlt.push_back(i+1); if (in[i]=="ALT")		{	FldAlt.clear(); FldAlt.push_back(i+1); }
			}
			if (AddInf && FldInf.no_input()) exit_error("The INFO column is missing.");
			if (Wt_VQS && FldInf.no_input()) exit_error("The INFO  column is missing.");
			if (Wt_DEL && ColDel==-2) exit_error("The BayesDel annotation is missing.");
			if (perch::LFonly && FldSym.no_input()) exit_error("The Gene Symbol column is missing.");
			if (perch::LFonly && FldFun.no_input()) exit_error("The Function column is missing.");
			if (perch::filSAF && FldInf.no_input()) exit_error("The INFO column is missing.");
			if (perch::filFAF && FldInf.no_input()) exit_error("The INFO column is missing.");
			if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing."); // FldChr is required for autosomal / sexlink linkage analysis
			if (FldPos.no_input()) exit_error("The POS/Start column is missing."); // FldPos is required for autosomal / sexlink linkage analysis
			if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
			if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
			FldIdx = FldChr + FldPos + FldRef + FldAlt;
			if (perch::VarCla)	FldGrp = FldIdx;
			else		{ if (FldGrp.no_input()) exit_error("The groupID column is missing."); }
			if (AddInf)	FldRes=FldInf;
			else	  { FldRes.push_back(in.NumFields()+1); in.contents().push_back(wrWhat); format.set_field_nums(FldRes,"",tfile_format::Expand); }

			int found_sid=0;
			for (int i=0;i<in.NumFields();i++) if (exist_element(fld4ind,in[i])) { fld4ind[in[i]]=i; ++found_sid; }
			lns<<showl<<found_sid<<" sequenced samples in Genotype File "<<vcf_file<<flush_logger;
			for (each_element(L2SID,it)) part2.push_back(fld4ind[*it]);
			
			if (cap_pg)
			{
				set<string> SeqPed;
				for (size_t i=0;i<L2SID.size();++i) if (fld4ind[L2SID[i]]) SeqPed.insert(L2FID[i]);
				if (perch::penetr[1]==perch::penetr[2])	max_pg=SeqPed.size();
				else									max_pg=SeqPed.size()*2;
			}

			print_container(in.contents(),program.outf,DLMTR,true);
			program.main_data().clear();
			continue;
		}
		if (header_not_read) exit_error("Header lines missing.");

		// read FORMAT
		genepi::gtp_par gpar;
		if (!FldFmt.no_input()) gpar.read(in[FldFmt[0]]);
		
		// check file
		if (in[FldAlt[0]].find(',')!=std::string::npos) exit_error("The input file has not been split by alternative alleles.");

		// group & index
		string this_group, this_index;
		FldGrp.contents_to_a_string(in.contents(),this_group,DLMTR);
		FldIdx.contents_to_a_string(in.contents(),this_index,DLMTR);
		if (this_group!=prev_group)
		{
			if (exist_element(old_groups,this_group)) elog.add(warning2, this_group); // exit_error("The Genotype File is not sorted by groups.");
			else old_groups.insert(this_group);
			if (!prev_group.empty())
			{
				write_results(AddInf,FldRes,wrWhat,MisStr,false);
				program.main_data().clear_ExceptTheLastRow();
			}
			prev_group=this_group;
			prev_index.clear();
			variant_lnkLH1.clear();
			variant_lnkLH0.clear();
			variant_lnkRes.clear();
			variant_weight.clear();
			variant_BayesD.clear();
			variant_filter.clear();
		}
		else
		{
			if (exist_element(prev_index,this_index))
			{
				variant_lnkLH1.push_back(0);
				variant_lnkLH0.push_back(0);
				variant_lnkRes.push_back(0);
				variant_weight.push_back(0);
				variant_BayesD.push_back(0);
				variant_filter.push_back("SKIPPED:DupVariant");
				elog.add(warning1);
				continue;
			}
		}
		prev_index.insert(this_index);

		// determine whether is_snv
		string& ref = in[FldRef[0]];
		string& alt = in[FldAlt[0]];
		bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
		bool is_lof = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"LoF")||str_has(in[FldFun[0]],"NMD") );
		bool is_vks = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"knClinSig=1") );
		bool is_DNg = ( FldFun.no_input() ? false : perch::is_DomNeg(in[FldFun[0]]) );
		bool is_cds = ( FldFun.no_input() ? false : perch::is_coding(in[FldFun[0]]) );

		vector<string> INFO;
		if (!FldInf.no_input())
		{
			if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
			bool info_modified=false;
			for (vector<string>::iterator it = INFO.begin(); it != INFO.end(); it++)
				if (str_startsw(*it,wrWhat+"=")) { INFO.erase(it); info_modified=true; break; }
			if (info_modified)
			{
				in[FldInf[0]] = str_of_container(INFO,';');
				if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
			}
		}

		// calculate well cover flag
		// skip by chr_num
		int chr_num = genepi::read_chr_num(in[FldChr[0]]);
		int bp = -1; try { bp = boost::lexical_cast<int>(in[FldPos[0]]); } catch (...) { exit_error("Failed to read "+in[FldPos[0]]+" as basepairs."); }

		string GeneSymbol;
		if (!FldSym.no_input()) GeneSymbol=in[FldSym[0]];
		if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
		if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }

		// MaxAF
		double MaxAF = std::numeric_limits<double>::signaling_NaN();
		if (!FldXAF.no_input())	read_val(in[FldXAF[0]],MaxAF);
		else					MaxAF=get_value(INFO,perch::h_MxAF);
		for (auto &id : perch::h_afID)
		{
			double x = get_value(INFO,id);
			if ( fabs(x-0.5)<fabs(MaxAF-0.5) || std::isnan(MaxAF)) MaxAF=x;
		}
		
		// Calculate MAF and AsIs (MA is alt or ref)
		double af = (force_defMAF ? defMAF : ( std::isnan(MaxAF) ? defMAF : MaxAF ));
		if (af==0) af=defMAF;
		if (af==1) af=1-defMAF;
		bool AsIs = (af<=0.5);
		double MAF = ( AsIs ? af : 1-af );

		// prepare num_car, num_gtp, proband
		bool has_missing = false;				// has missing value in at least one sequence sample
		map<string,int>		num_car;			// num_car[FamID] = number of carriers
		map<string,int>		num_gtp;			// fam_pos[FamID] = number of genotyped
		map<string,int>		proband = fid2prb;	// proband[FamID] = line number of the proband (the carrier with the most severe liability class (highest RR))
		map<string,double>	prob_rr;			// proband[FamID] = proband relative risk
		for (unsigned i=0;i<part1.size();++i)
		{
			int	carrier = 0; // carrier of the variant
			int	nonmiss = 0; // successfully genotyped
			if (part2[i])
			{
				// geno allele is either 0 or 1. It's not a problem as Non-ALT alleles are always treated as not the allele of interest.
				char g = GTP::read(in[part2[i]], chr_num, bp, L2SEX[i], gpar);
				if (!AsIs) GTP::swap_allele(g);
				if		(!GTP::usable(g))	has_missing=true;
				else if (GTP::num_alt(g))	{ ++carrier; ++nonmiss; }
				else						{			 ++nonmiss; }
				if (carrier) // previously "carrier && L2AFF[i]", less powerful and not correct
				{
					if (!exist_element(fid2prb,L2FID[i]))
					{
						if ((L2LIA[i]+1)>int(relative_risk.size())) exit_error("A liability class code exceeds the number of liability classes.");
						double this_rr=relative_risk[L2LIA[i]];
						if (!exist_element(proband,L2FID[i]))	{ proband[L2FID[i]]=i; prob_rr[L2FID[i]]=this_rr; }
						else if (this_rr>prob_rr[L2FID[i]])		{ proband[L2FID[i]]=i; prob_rr[L2FID[i]]=this_rr; }
					}
					++num_car[L2FID[i]];
				}
				if (nonmiss)
					++num_gtp[L2FID[i]];
			}
			if (perch::VarCla)
			{
				if (exist_element(proband,L2FID[i]) && proband[L2FID[i]]==(int)i)
					if (!carrier) exit_error("The designated proband is not successfully genotyped or not a carrier.");
			}
		}
		if (perch::VarCla)
		{
			for (auto &f:proband)
			{
				if (!exist_element(fid2prb,f.first))
					lns<<showl<<"selected individual "<<L2SID[f.second]<<" as a proband for pedigree "<<f.first<<flush_logger;
			}
		}
		
		// write prefile_V
		int num_ind=0;
		stringstream prefile_V;
		for (unsigned i=0;i<part1.size();++i)
		{
			if (num_car[L2FID[i]]<min_AnyCar) continue;
			if (num_gtp[L2FID[i]]<min_GtpPer) continue;
			++num_ind;
			prefile_V << part1[i];
			if (part2[i])
			{
				// geno allele is either 0 or 1. It's not a problem as Non-ALT alleles are always treated as not the allele of interest.
				char g = GTP::read(in[part2[i]], chr_num, bp, L2SEX[i], gpar);
				if (!AsIs) GTP::swap_allele(g);
				prefile_V << '\t' << GTP::to_linkage(g) << '\n';
			}
			else
			{
				prefile_V << "\t0 0\n";
			}
		}

		// write prefile_P
		stringstream prefile_P;
		for (unsigned i=0;i<part1.size();++i)
		{
			if (num_car[L2FID[i]]<min_AnyCar) continue;
			if (num_gtp[L2FID[i]]<min_GtpPer) continue;
			prefile_P << part1[i];
			if (part2[i] && proband[L2FID[i]]==(int)i)
			{
				// geno allele is either 0 or 1. It's not a problem as Non-ALT alleles are always treated as not the allele of interest.
				char g = GTP::read(in[part2[i]], chr_num, bp, L2SEX[i], gpar);
				if (!AsIs) GTP::swap_allele(g);
				prefile_P << '\t' << GTP::to_linkage(g) << '\n';
			}
			else
			{
				prefile_P << "\t0 0\n";
			}
		}
		
		bool passed_filters = true;
		string applied_filters;
		
		// skip by region
		if (!perch::within_covered_region(chr_num,bp)) { passed_filters=false; applied_filters+="Region,"; }

		// skip by var type
		if (SNonly && !is_snv)	{ passed_filters=false; applied_filters+="NotSNV,"; }
		if (perch::LFonly && is_lof && perch::is_LoFtol(GeneSymbol))	{ passed_filters=false; applied_filters+="LoF-tolerant,"; }
		if (perch::LFonly && !is_lof)	{ passed_filters=false; applied_filters+="NotLoF,"; }
		if (perch::CDonly && !is_cds)	{ passed_filters=false; applied_filters+="NotCDS,"; }
		if (perch::DomNeg && !is_DNg)	{ passed_filters=false; applied_filters+="NotDomNeg,"; }

		// skip by VICTOR_QC
		if (!FldInf.no_input())
		{
			string vQC = get_string(INFO,"VICTOR_QC");
			if (!vQC.empty() && vQC!="PASS") { passed_filters=false; applied_filters+="vQC="+vQC+","; }
		}
		
		// skip by VQSLOD
		double VQSLOD = get_value(INFO,"VQSLOD");
		/*if (std::isnan(VQSLOD))
		{
			if (perch::VQSnan)	{ passed_filters=false; applied_filters+="noVQSLOD,"; }
			if (perch::HFnoVQ || perch::hardft)
			{
				double infoQD = get_value(INFO,"QD"); // use double so that it is nan if not exist
				double infoMQ = get_value(INFO,"MQ"); // nan automatically return false in comparison
				double infoFS = get_value(INFO,"FS"); // so no need to check isnan()
				double infoHS = get_value(INFO,"HaplotypeScore");
				double infoMR = get_value(INFO,"MQRankSum");
				double infoRP = get_value(INFO,"ReadPosRankSum");
				if (perch::FiltQD) { if (infoQD<perch::FiltQD) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltMQ) { if (infoMQ<perch::FiltMQ) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltFS) { if (infoFS>perch::FiltFS) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltHS) { if (infoHS>perch::FiltHS) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltMR) { if (infoMR<perch::FiltMR) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltRP) { if (infoRP<perch::FiltRP) { passed_filters=false; applied_filters+="HardFilter,"; } }
			}
		}
		else
		{
			if (is_snv) { if (VQSLOD<perch::VQSsnv) { passed_filters=false; applied_filters+="VQSLOD,"; } }
			else		{ if (VQSLOD<perch::VQSidl) { passed_filters=false; applied_filters+="VQSLOD,"; } }
			if (perch::hardft)
			{
				double infoQD = get_value(INFO,"QD"); // use double so that it is nan if not exist
				double infoMQ = get_value(INFO,"MQ"); // nan automatically return false in comparison
				double infoFS = get_value(INFO,"FS"); // so no need to check isnan()
				double infoHS = get_value(INFO,"HaplotypeScore");
				double infoMR = get_value(INFO,"MQRankSum");
				double infoRP = get_value(INFO,"ReadPosRankSum");
				if (perch::FiltQD) { if (infoQD<perch::FiltQD) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltMQ) { if (infoMQ<perch::FiltMQ) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltFS) { if (infoFS>perch::FiltFS) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltHS) { if (infoHS>perch::FiltHS) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltMR) { if (infoMR<perch::FiltMR) { passed_filters=false; applied_filters+="HardFilter,"; } }
				if (perch::FiltRP) { if (infoRP<perch::FiltRP) { passed_filters=false; applied_filters+="HardFilter,"; } }
			}
		} //*/
		
		/*/ skip by FILTER
		if (!FldFlt.no_input())
			if (!perch::filflt.empty() && !exist_element(perch::filflt,in[FldFlt[0]])) { passed_filters=false; applied_filters+="FILTER,"; } //*/

		// skip by BayesDel
		double BayesDel = std::numeric_limits<double>::signaling_NaN();
		if		(ColDel==-1)	BayesDel = get_value_sw(INFO,"BayesDel");
		else if (ColDel>=0)	read_val(in[ColDel],BayesDel);
		if (perch::filter_AnnAF(BayesDel,GeneSymbol,is_lof,is_vks)) { passed_filters=false; applied_filters+="BayesDel,"; }
		
		// skip by cap
		if (max_pg)
		{
			if (variant_BayesD.size()>=max_pg) { passed_filters=false; applied_filters+="VarPerGene,"; }
			if (!variant_BayesD.empty() && !std::isnan(BayesDel))
			{
				if (std::isnan(variant_BayesD.back())) exit_error("in sorting Genotype File by BayesDel, NA should be treated as the minimum value.");
				else if (variant_BayesD.back()<BayesDel) exit_error("the Genotype File is not sorted by BayesDel. See "+this_group);
			}
		}

		// skip by MAF
		if (perch::filXAF)
		{
			if (MaxAF<=0.5 &&    MaxAF >perch::filXAF) { passed_filters=false; applied_filters+="MaxAF,"; }
			if (MaxAF >0.5 && (1-MaxAF)>perch::filXAF) { passed_filters=false; applied_filters+="MaxAF,"; }
		}
		if (perch::filSAF)
		{
			double SplAF = get_value(INFO,"SplAF");
			if (std::isnan(SplAF)) exit_error("cannot apply --filt-SplAF because SplAF is not calculated by vQC");
			if (SplAF>0.5) SplAF=1-SplAF;
			if (SplAF>perch::filSAF) { passed_filters=false; applied_filters+="SplAF,"; }
		}
		if (perch::filFAF)
		{
			double FdrAF = get_value(INFO,"AF_Founder");
			double FdrAN = get_value(INFO,"AN_Founder");
			if (!std::isnan(FdrAF) && !std::isnan(FdrAN))
			{
				if (FdrAF>0.5) FdrAF=1-FdrAF;
				if (FdrAF>perch::filFAF && FdrAN>=perch::filFAFminAN) { passed_filters=false; applied_filters+="AF_Founder,"; }
			}
		}

		// skip by missing rate
		if (perch::MisCut!=1)
		{
			// try { double MsgCSs=get_value(INFO,"MissingInCs"); if (MsgCSs>perch::MisCut) passed_filters=false; } catch (...) {}
			// try { double MsgCTs=get_value(INFO,"MissingInCt"); if (MsgCTs>perch::MisCut) passed_filters=false; } catch (...) {}
			try { double MsgAll=get_value(INFO,"MissingRate"); if (MsgAll>perch::MisCut) { passed_filters=false; applied_filters+="MissingRate,"; } } catch (...) {}
		}
		
		// parameters and default results
		double th(0),al(1),lnkRes(0),lnkLH1(0),lnkLH0(0); // recombination rate, % of linkage, LOD, log10 Likelihood H1, log10 Likelihood H0
		int SexLnk = (genepi::is_chrX(chr_num) && !genepi::is_autoOrPAR(chr_num,bp)); // chrX
		if ( (!has_missing || !SkipPV_mis) && num_ind && passed_filters && (genepi::is_autosomal(chr_num) || genepi::is_chrX(chr_num)) )
		{
			stringstream pedfile_V, ipedfileV, speedfilV;			
			int unknown_errors=0;
			int MUT_LOCUS = 0;
			if (mut_m || mut_f) MUT_LOCUS=1;

			// write datafile0
			stringstream datafile0;
			datafile0<< "2 0 "<<SexLnk<<" 5\n"<<MUT_LOCUS<<" "<<mut_m<<" "<<mut_f<<" 0\n1 2\n1 2\n"<<1-MAF<<' '<<MAF<<'\n';
			if (SexLnk) {
				if (liab_chrX.empty())	datafile0<<"1\n"<<perch::penetr[0]<<" "<<perch::penetr[1]<<" "<<perch::penetr[2]<<"\n"<<perch::penetr[0]<<" "<<perch::penetr[2]<<"\n";
				else					datafile0<<liab_chrX; }
			else {
				if (liability.empty())	datafile0<<"1\n"<<perch::penetr[0]<<" "<<perch::penetr[1]<<" "<<perch::penetr[2]<<"\n";
				else					datafile0<<liability; }
			datafile0<<"3 2\n"<<1-MAF<<'\t'<<MAF<< "\n0 0\n0\n1 1 0.5\n";

			// write datafile1
			stringstream datafile1;
			datafile1<< "2 0 "<<SexLnk<<" 5\n"<<MUT_LOCUS<<" "<<mut_m<<" "<<mut_f<<" 1\n1 2\n1 2\n";
			if (SexLnk) {
				if (liab_chrX.empty())	datafile1<<"1\n"<<perch::penetr[0]<<" "<<perch::penetr[1]<<" "<<perch::penetr[2]<<"\n"<<perch::penetr[0]<<" "<<perch::penetr[2]<<"\n";
				else					datafile1<<liab_chrX; }
			else {
				if (liability.empty())	datafile1<<"1\n"<<perch::penetr[0]<<" "<<perch::penetr[1]<<" "<<perch::penetr[2]<<"\n";
				else					datafile1<<liability; }
			datafile1<<"3 2\n"<<1-MAF<<" 0 0 "<<MAF<< "\n0 0\n0\n1 1 0.5\n";
			
			// prepare V
			if (perch::_Debug)
			{
				fstream lf;
				openfile_or_exit(lf,"vSEG.log",ios::out|ios::app);
				lf << "## Data for " << this_group << " : " << this_index << endl;
				lf << "## datafile0 " << endl << datafile0.str();
				lf << "## datafile1 " << endl << datafile1.str();
				lf << "## prefile_V " << endl << prefile_V.str();
				lf << "## prefile_P " << endl << prefile_P.str();
				closefile(lf);
			}
			rewind(datafile0);
			makeped_program(prefile_V,pedfile_V,true);
			unknown_program(datafile0,pedfile_V,ipedfileV,speedfilV,unknown_errors);
			if (unknown_errors)
			{
				lnkRes = std::numeric_limits<double>::signaling_NaN();
			}
			else // do linkage
			{
				{
					double LL_0_5_V, LL_1_0_V;
					th=0.5; rewind(datafile0); rewind(ipedfileV); rewind(speedfilV); mlink_calculate(datafile0,ipedfileV,speedfilV,th,al,LL_0_5_V);
					th=0.0; rewind(datafile1); rewind(ipedfileV); rewind(speedfilV); mlink_calculate(datafile1,ipedfileV,speedfilV,th,al,LL_1_0_V);
					lnkLH1 = LL_1_0_V;
					lnkLH0 = LL_0_5_V;
					if (perch::_Debug) cerr << LL_1_0_V << ' ' << LL_0_5_V << endl;
				}
				
				if (AscertCorr)
				{
					stringstream pedfile_P, ipedfileP, speedfilP;
					rewind(datafile1);
					makeped_program(prefile_P,pedfile_P,true);
					unknown_program(datafile1,pedfile_P,ipedfileP,speedfilP,unknown_errors);
					double LL_0_5_P, LL_1_0_P;
					th=0.5; rewind(datafile0); rewind(ipedfileP); rewind(speedfilP); mlink_calculate(datafile0,ipedfileP,speedfilP,th,al,LL_0_5_P);
					th=0.0; rewind(datafile1); rewind(ipedfileP); rewind(speedfilP); mlink_calculate(datafile1,ipedfileP,speedfilP,th,al,LL_1_0_P);
					lnkLH1 -= LL_1_0_P;
					lnkLH0 -= LL_0_5_P;
					if (perch::_Debug) cerr << LL_1_0_P << ' ' << LL_0_5_P << endl;
				}
				lnkRes = lnkLH1 - lnkLH0; // lnkRes = (LL_1_0_V - LL_0_5_V) + (LL_0_5_P - LL_1_0_P);
			}
		}
		else
		{
			if (has_missing && SkipPV_mis) applied_filters+="MissingNess,";
			if (num_ind==0) applied_filters+="NoData,";
			if (!genepi::is_autosomal(chr_num) && !genepi::is_chrX(chr_num)) applied_filters+="Chr,";
		}
		if (str_endsw(applied_filters,",")) applied_filters.pop_back();
		
		// weight
		double L10BF = 0; // BF in log10 scale
		if (Wt_VQS) { double v = (VQSLOD<0?VQSLOD:0);	if (!std::isnan(v)) L10BF+=v; }
		if (Wt_DEL) { double v = BayesDel;				if (!std::isnan(v)) L10BF+=v; }
		double weight = posterior_given_log10BF(L10BF);

		variant_lnkLH1.push_back(lnkLH1);
		variant_lnkLH0.push_back(lnkLH0);
		variant_lnkRes.push_back(lnkRes);
		variant_weight.push_back(weight);
		variant_BayesD.push_back(BayesDel);
		variant_filter.push_back(applied_filters);
		
		if (perch::_Debug) cerr << num_ind <<'\t'<< MAF <<'\t'<< weight <<'\t'<< lnkRes << endl;
	}
	write_results(AddInf,FldRes,wrWhat,MisStr,true);

	return 0;
}
