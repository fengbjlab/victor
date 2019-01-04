/*
 If all deleteriousness scores are NA, vDEL still tries to annotate based on the function:
 For coding variants other than missense variants (exclude synonymous), vDEL assigns a DELlbf equal to the worst case missense variants from the same gene.
 For missense variants or substitution variants in a non-coding gene, DELlbf is the median value of the gene.
 
 Removed option descriptions:
 --col-symb STR    In input, header of the Gene Symbol      col. {_Default_col_symb}
 --col-func STR    In input, header of the Func Consequence col. {_Default_col_func}
 --col-1 STRs      In input, header of the first column {_Default_col_one}
 --wt FILE         get weights        fr FILE(Method weight) or wt=1  [internal]
 --bf FILE         get Bayes factor   fr FILE(Method File(Score  BF)) [internal]$
 --ws FILE         get worst score    fr FILE(Method File(Symbol WS)) [internal]#
 --ct FILE         get CADD-to-DELlbf fr FILE(CADD DELlbf)             [internal]
 --out-name STR    in output, field header or INFO variable name is STR [DELlbf]
 --add-info        in output, add DELlbf to the INFO field [add a new field]
 --convert         in output, convert scores to log-scale Bayes factors
 --update          in output, update missing scores / non-missense genes
 --ms STR          Output string for missing values [NA]
 --add-af B        add MaxAF to BayesDel {_Default_add_af}
 --ensemble STR    component score set
 --MaxAF-pr D      Integrate MaxAF into the annotated score with a proportion of D {_Default_MaxAF_pr}
 --MaxAF-no0 B     Change MaxAF=0 to MaxAF=missing, otherwise change MaxAF=missing to MaxAF=0 {_Default_MaxAF_no0}
 --CADD-db F       The CADD whole_genome_SNVs.tsv.gz file full path {_Default_cadd_db}
 --med-for STRs    If the functional consequence has any of the STRs, assign the median score among all possible missense variants within the gene. This applies when the combined BayesDel cannot be calculated.  {_Default_ms_for}
 $  If used --bf, weights become all 1. So better use --wt together.
 #  For Method it's the worst score; unrecognized method is DELlbf.
 */

#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_program.hpp>
#include "vDEL.hpp"
#include "victor_par.hpp"

using namespace std;

class DeleteriousScore {
public:
	int			column; // -2:unknown -1:INFO 0+:other_column
	Interpolate	curve;
	DeleteriousScore():column(-2) {}
};

double get_alternative(const string& method,
					   const string& symbol,
					   const map<string, map<string,double> >& ref,
					   double default_value)
{
	map<string, map<string,double> >::const_iterator it1 = ref.find(method);
	if (it1==ref.end())	return default_value;
	map<string,double>::const_iterator it2 = it1->second.find(symbol);
	if (it2==it1->second.end())	return default_value;
	return it2->second;
}

// output lines sorted by del & maf within symbol
// symbol="" has a special meaning: print everything in memory and add nothing.
// this function is not thread-safe due to static variables, but vDEL doesn't need multi-threading
void add_line_by_gene(const string& line, double del, double maf, const string& symbol)
{
	static map<double, multimap<double,string>, std::greater<double> > stored_lines;
	static string stored_symbol;
	if (stored_symbol!=symbol)
	{
		for (auto &d:stored_lines)
			for (auto &l:d.second) program.outf<<l.second<<endl;
		stored_lines.clear();
		stored_symbol=symbol;
	}
	if (std::isnan(del)) del=-INFINITY;
	if (std::isnan(maf)) maf=1;
	if (!symbol.empty()) stored_lines[del].insert(make_pair(maf,line));
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	string				wrWhat="BayesDel";
	vector<string>		inputs;			// input VCF or annotated file
	tfile_format		format;			// format of input
	string				bfFile;			// Bayes factor
	string				wtFile;			// weight
	string				MisStr="00";	// missing value (each: non-numeric; combined: all non-numeric), use 0 instead of NA, good for sorting
	string				WrMsTo;			// check missing CADD, PROVEAN, fathmmMKL and write variants to S.for_CADD S.for_PROV S.for_fMKL, respectively
	bool				PRVloc=true;	// write files to run PROVEAN locally
	bool				AddInf=false;	// add a field in INFO, otherwise add a column
	bool				UpdDel=false;	// update the deleterious score (MisStr and non-missense)
	bool				ChgDel=false;	// change the deleterious score to log-scale bf (for training)
	bool				NoFunc=false;	// no function-aware
	bool				NoSort=false;	// no sorting output lines
	bool				wr_res=true;	// write result
	bool				af0toM=false;	// change MaxAF=0 to nan, otherwise change MaxAF=nan to 0.
	bool				add_af=false;
	double				m0_cut=0.6;		// Use 0 to disable model 1,2 unless model 0 is completely empty. Use 1 to disable model 0 and enable model 1,2.
	string				wt_set="nsfp33a_noAF";// nsfp23 nsfp23_noMAF nsfp26 nsfp26_noMAF
	
	// Below are strings to identify variant types by str_has(). Case-insensitive (will change to lower later). Input should not have any "non-frameshift" because it matches to both.
	// Default WstFor values work for ANNOVAR, but for vAnnGene it should be --wst-for=LoF. To use scSNV for splice-altering prediction, insert LoF to let vDel read it.
	vector<string>		WstFor={"LoF","NMD"}; // No "transcription","translation" because they must be LoF; removed "StartCodon" "initiator" "frameshift" "ncRNA_InDel"
	vector<string>		MedFor={"stopgain","stop-gain","nonsense","splic","TrEff","nonsyn","non-syn","missense","inframe","in-frame","nonframe","non-frame","stoplos","stop-los"}; // rm "miRNA_Bind","ncRNA_sub","ncRNA_SNV"
	
	program.read_arguments(argc,argv);
	perch::read_arguments();
	format .set_delimiters("\t");
	format .forbid_option(SKIP_NOTES);
	format .read_arguments(program.arg());
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--bf"))				ReadArg(program.arg(),argi,bfFile);
		else if (str_startsw(program.arg()[argi],"--wt"))				ReadArg(program.arg(),argi,wtFile);
		else if (str_startsw(program.arg()[argi],"--ms"))				ReadArg(program.arg(),argi,MisStr);
		else if (str_startsw(program.arg()[argi],"--wst-for"))		{	ReadSet(program.arg(),argi,WstFor); if (!WstFor.empty()) NoFunc=false; }
		else if (str_startsw(program.arg()[argi],"--med-for"))		{	ReadSet(program.arg(),argi,MedFor); if (!MedFor.empty()) NoFunc=false; }
		else if (str_startsw(program.arg()[argi],"--no-func"))			ReadArg(program.arg(),argi,NoFunc);
		else if (str_startsw(program.arg()[argi],"--out-name"))			ReadArg(program.arg(),argi,wrWhat);
		else if (str_startsw(program.arg()[argi],"--add-info"))			ReadArg(program.arg(),argi,AddInf);
		else if (str_startsw(program.arg()[argi],"--convert"))			ReadArg(program.arg(),argi,ChgDel);
		else if (str_startsw(program.arg()[argi],"--update"))			ReadArg(program.arg(),argi,UpdDel);
		else if (str_startsw(program.arg()[argi],"--ensemble"))			ReadArg(program.arg(),argi,wt_set);
		else if (str_startsw(program.arg()[argi],"--check-ms"))			ReadArg(program.arg(),argi,WrMsTo);
		else if	(str_startsw(program.arg()[argi],"--local-provean"))	ReadArg(program.arg(),argi,PRVloc);
		else if (str_startsw(program.arg()[argi],"--no-sort"))			ReadArg(program.arg(),argi,NoSort);
		else if (str_startsw(program.arg()[argi],"--wr-res"))			ReadArg(program.arg(),argi,wr_res);
		else if (str_startsw(program.arg()[argi],"--m0-cut"))			ReadArg(program.arg(),argi,m0_cut);
		else if (str_startsw(program.arg()[argi],"--MaxAF-no0"))		ReadArg(program.arg(),argi,af0toM);
		else if	(str_startsw(program.arg()[argi],"--MaxAF-pr"))			ReadArg(program.arg(),argi,AF_prp);
		else if	(str_startsw(program.arg()[argi],"--add-af"))			ReadArg(program.arg(),argi,add_af);
		else if (str_startsw(program.arg()[argi],"-"))					exit_error("unknown option "+program.arg()[argi]);
		else add_to_container(inputs,program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_AddInfo",str_YesOrNo(AddInf));
	program.help_text_var("_Default_ws_for",str_of_container(WstFor,',',false));
	program.help_text_var("_Default_ms_for",str_of_container(MedFor,',',false));
	program.help_text_var("_Default_ensemble",wt_set);
	program.help_text_var("_Default_no_sort",str_YesOrNo(NoSort));
	program.help_text_var("_Default_wr_res",str_YesOrNo(wr_res));
	program.help_text_var("_Default_MaxAF_no0",str_YesOrNo(af0toM));
	program.help_text_var("_Default_loc_prov",str_YesOrNo(PRVloc));
	program.help_text_var("_Default_add_af",str_YesOrNo(add_af));
	program.help_text_var("_Default_check_ms",WrMsTo);
	program.help_text_var("_Default_MaxAF_pr",ftos(AF_prp));
	program.help_text_var("_Candidate_ensemble","nsfp33a nsfp33a_noAF nsfp26 nsfp26_noAF nsfp23 nsfp23_noAF");
	program.manual += "nsfp23 components: "+str_of_map_first(ds_crv_LJB23,' ',true);
	program.manual += "nsfp26 components: "+str_of_map_first(ds_crv_LJB26,' ',true);
	program.manual += "nsfp33a components: "+str_of_map_first(ds_crv_nsfp33a,' ',true);
	perch::check_arguments();
	
	// check errors
	if (NoFunc) { WstFor.clear(); MedFor.clear(); }
	if (!wr_res) NoSort=true;
	if (UpdDel && ChgDel) exit_error("Use either --update or --convert, but not both.");
	if (UpdDel && ChgDel && MisStr=="00") exit_error("With either --update or --convert, it's very likely you don't want --ms to be 00 by default.");
	for (auto &x:WstFor) if (x.empty()) exit_error("--wst-for elements cannot be empty");
	for (auto &x:MedFor) if (x.empty()) exit_error("--med-for elements cannot be empty");
	for (auto &x:WstFor) for (auto &y:WstFor) if (x!=y && str_has(x,y)) exit_error("Redundant --wst-for: "+x+" "+y);
	for (auto &x:MedFor) for (auto &y:MedFor) if (x!=y && str_has(x,y)) exit_error("Redundant --med-for: "+x+" "+y);
	for (auto &x:WstFor) for (auto &y:MedFor) if (str_has(y,x) || str_has(x,y)) exit_error("Overlap between --wst-for and --med-for: "+x+" "+y);
	for (auto &x:WstFor) boost::to_lower(x);
	for (auto &x:MedFor) boost::to_lower(x);
	boost::replace_all(wt_set, "noMAF", "noAF");
	
	bool write_missing_del = !WrMsTo.empty();
	std::fstream ms_PROV; // ms_CADD, ms_fMKL;
	if (write_missing_del)
	{
		openfile_or_exit(ms_PROV,WrMsTo+".for_PROV",ios::out);
		// openfile_or_exit(ms_CADD,WrMsTo+".for_CADD",ios::out);
		// openfile_or_exit(ms_fMKL,WrMsTo+".for_fMKL",ios::out);
	}
	Interpolate	MaxAF_curve;
	MaxAF_curve.reference_data=ebf_nsfp33a_MaxAF;
	
	map<string, double> m0_wtA, m0_wtB, m1_wtA, m1_wtB, m2_wtA, m2_wtB;
	map<string, map<double,double> > m0_crv, m1_crv, m2_crv;
	double m12_cut = std::numeric_limits<double>::signaling_NaN();
	if		(wt_set=="nsfp23")		 { m0_wtA=ds_wts_LJB23;			m0_crv=ds_crv_LJB23;			m0_wtB=ds_wts_LJB23_noMAF;	}
	else if (wt_set=="nsfp23_noAF")  { m0_wtA=ds_wts_LJB23_noMAF;	m0_crv=ds_crv_LJB23_noMAF;		m0_wtB=ds_wts_LJB23_noMAF;	}
	else if	(wt_set=="nsfp26")		 { m0_wtA=ds_wts_LJB26;			m0_crv=ds_crv_LJB26;			m0_wtB=ds_wts_LJB26_noMAF;	}
	else if (wt_set=="nsfp26_noAF")  { m0_wtA=ds_wts_LJB26_noMAF;	m0_crv=ds_crv_LJB26_noMAF;		m0_wtB=ds_wts_LJB26_noMAF;	}
	else if	(wt_set=="nsfp33a")		 { m0_wtA=ds_wts_nsfp33a;		m0_crv=ds_crv_nsfp33a;			m0_wtB=ds_wts_nsfp33a_noMAF;}
	else if (wt_set=="nsfp33a_noAF") { m0_wtA=ds_wts_nsfp33a_noMAF;	m0_crv=ds_crv_nsfp33a_noMAF;	m0_wtB=ds_wts_nsfp33a_noMAF;}
	else exit_error("wrong --ensemble argument \""+wt_set+"\"");

	if		(wt_set=="nsfp23")		 { m1_wtA=ds_wts_LJB23_noMAF;	m1_crv=ds_crv_LJB23_noMAF;		m1_wtB=ds_wts_LJB23_noMAF;	} // nsfp23 has no alt model, so I use _noMAF instead.
	else if (wt_set=="nsfp23_noAF")  { m1_wtA=ds_wts_LJB23_noMAF;	m1_crv=ds_crv_LJB23_noMAF;		m1_wtB=ds_wts_LJB23_noMAF;	}
	else if	(wt_set=="nsfp26")		 { m1_wtA=m1_wts_LJB26;			m1_crv=m1_crv_LJB26;			m1_wtB=m1_wts_LJB26;		} // nsfp26 has one version of alt model, no MAF.
	else if (wt_set=="nsfp26_noAF")  { m1_wtA=m1_wts_LJB26;			m1_crv=m1_crv_LJB26;			m1_wtB=m1_wts_LJB26;		}
	else if	(wt_set=="nsfp33a")		 { m1_wtA=m1_wts_nsfp33a;		m1_crv=m1_crv_nsfp33a;			m1_wtB=m1_wts_nsfp33a_noMAF;}
	else if (wt_set=="nsfp33a_noAF") { m1_wtA=m1_wts_nsfp33a_noMAF; m1_crv=m1_crv_nsfp33a_noMAF;	m1_wtB=m1_wts_nsfp33a_noMAF;}
	else exit_error("wrong --ensemble argument \""+wt_set+"\"");

	if		(wt_set=="nsfp23")		 { m2_wtA=m2_wts_LJB26;			m2_crv=m2_crv_LJB26;			m2_wtB=m2_wts_LJB26;			m12_cut=m12_cut_LJB26;			}
	else if (wt_set=="nsfp23_noAF")  { m2_wtA=m2_wts_LJB26;			m2_crv=m2_crv_LJB26;			m2_wtB=m2_wts_LJB26;			m12_cut=m12_cut_LJB26;			}
	else if	(wt_set=="nsfp26")		 { m2_wtA=m2_wts_LJB26;			m2_crv=m2_crv_LJB26;			m2_wtB=m2_wts_LJB26;			m12_cut=m12_cut_LJB26;			}
	else if (wt_set=="nsfp26_noAF")  { m2_wtA=m2_wts_LJB26;			m2_crv=m2_crv_LJB26;			m2_wtB=m2_wts_LJB26;			m12_cut=m12_cut_LJB26;			}
	else if	(wt_set=="nsfp33a")		 { m2_wtA=m2_wts_nsfp33a;		m2_crv=m2_crv_nsfp33a;			m2_wtB=m2_wts_nsfp33a_noMAF;	m12_cut=m12_cut_nsfp33a;		}
	else if (wt_set=="nsfp33a_noAF") { m2_wtA=m2_wts_nsfp33a_noMAF;	m2_crv=m2_crv_nsfp33a_noMAF;	m2_wtB=m2_wts_nsfp33a_noMAF;	m12_cut=m12_cut_nsfp33a_noMAF;	}
	else exit_error("wrong --ensemble argument \""+wt_set+"\"");

	// read deleterious score => Bayes factor (not log10) curve
	map<string,DeleteriousScore> ds_dta;			// For combined BayesDel calculation
	map<string,Values<double> > ds_val;				// For combined BayesDel checking each individual predictor
	DeleteriousScore CADD, PROVEAN, fMKL_NC;// For Bayes factor calculation from other predictors
	if (!bfFile.empty())
	{
		m0_wtA.clear();
		for (Rows_in_File(in,bfFile,2))
		{
			DeleteriousScore ds;
			ds.curve.setup(in[1]);
			ds_dta[in[0]]=ds;
			m0_wtA[in[0]]=1;
		}
		m0_wtB = m0_wtA;
		wt_set = substr_before_find(bfFile,".");
		m0_crv.clear();
		m1_crv.clear();
		m2_crv.clear();
		m1_wtA.clear();
		m1_wtB.clear();
		m2_wtA.clear();
		m2_wtB.clear();
	}
	else
	{
		for (auto& i:m0_crv)
		{
			if (!exist_element(m0_wtA,i.first) && !exist_element(m0_wtB,i.first)) continue;
			if (m0_wtA[i.first]==0 && m0_wtB[i.first]==0) continue;
			DeleteriousScore ds;
			ds.curve.reference_data=i.second;
			ds_dta[i.first] = ds;
		}
		if (str_startsw(wt_set,"nsfp23"))
		{
			CADD.curve.reference_data=ebf_LJB23_CADD;
			PROVEAN.curve.reference_data=ebf_nsfp33a_PROVEAN;
			fMKL_NC.curve.reference_data=ebf_nsfp33a_fthmMKL_NC;
		}
		else if (str_startsw(wt_set,"nsfp26"))
		{
			CADD.curve.reference_data=ebf_LJB26_CADD;
			PROVEAN.curve.reference_data=ebf_nsfp33a_PROVEAN;
			fMKL_NC.curve.reference_data=ebf_nsfp33a_fthmMKL_NC;
		}
		else if (str_startsw(wt_set,"nsfp33a"))
		{
			CADD.curve.reference_data=ebf_nsfp33a_CADD;
			PROVEAN.curve.reference_data=ebf_nsfp33a_PROVEAN;
			fMKL_NC.curve.reference_data=ebf_nsfp33a_fthmMKL_NC;
		}
		else exit_error("wrong --ensemble argument "+wt_set);
	}
	for (auto& ds:ds_dta) if (ds.second.curve.reference_data.empty()) exit_error("Curve data is empty for "+ds.first);
	
	// For combined BayesDel calculation, read weight for each individual predictor
	if (!wtFile.empty())
	{
		m0_wtA.clear();
		for (Rows_in_File(in,wtFile,2))
		{
			try { m0_wtA[in[0]] = boost::lexical_cast<double>(in[1]); } catch (...) {}
		}
		set<string> missing;
		for (auto & i:ds_dta)
		{
			if (!exist_element(m0_wtA,i.first)) { missing.insert(i.first); m0_wtA[i.first]=0; }
		}
		if (!missing.empty())
		{
			lns<<showl<<str_of_container(missing,',')<<" are missing in "<<wtFile<<"; their weights are set to 0."<<flush_logger;
		}
		m0_wtB = m0_wtA;
	}
	
	// alternative models when there's missing
	map<string,DeleteriousScore> m1_dta, m2_dta;
	for (auto& i:m1_crv)
	{
		if (!exist_element(m1_wtA,i.first) && !exist_element(m1_wtB,i.first)) continue;
		if (m1_wtA[i.first]==0 && m1_wtB[i.first]==0) continue;
		DeleteriousScore ds;
		ds.curve.reference_data=i.second;
		m1_dta[i.first] = ds;
	}
	for (auto& i:m2_crv)
	{
		if (!exist_element(m2_wtA,i.first) && !exist_element(m2_wtB,i.first)) continue;
		if (m2_wtA[i.first]==0 && m2_wtB[i.first]==0) continue;
		DeleteriousScore ds;
		ds.curve.reference_data=i.second;
		m2_dta[i.first] = ds;
	}
	for (auto& ds:m1_dta) if (ds.second.curve.reference_data.empty()) exit_error("Curve data is empty for "+ds.first);
	for (auto& ds:m2_dta) if (ds.second.curve.reference_data.empty()) exit_error("Curve data is empty for "+ds.first);
	
	// read deleterious score statatistics (MAX MEDIAN) for each gene symbol
	map<string, map<string,double> > WorstCase; // WorstCase[method][symbol]
	map<string, map<string,double> > MedianVal; // MedianVal[method][symbol]
	if ((!WstFor.empty() || !MedFor.empty()) && bfFile.empty())
	{
		vector<string> methods = { wrWhat };
		for (auto &method:methods)
		{
			Values<double> acc_max, acc_med;
			for (Rows_in_File(in2,perch::find_file(method+"_stat.txt"),3)) // previously "@GDB_"+method+".txt"
			{
				if (in2[0]=="Transcript" || in2[0]=="Gene" || in2[0]=="ID") continue;
				double mx = 0; try { mx=boost::lexical_cast<double>(in2[1]); } catch (...) { exit_error("Error reading @GDB_"+method+".txt: failed to read "+in2[1]+" as a number."); }
				double md = 0; try { md=boost::lexical_cast<double>(in2[2]); } catch (...) { exit_error("Error reading @GDB_"+method+".txt.txt"+": failed to read "+in2[2]+" as a number."); }
				WorstCase[method][in2[0]] = mx; acc_max.push_back(mx);
				MedianVal[method][in2[0]] = md; acc_med.push_back(md);
			}
			WorstCase[method]["MeanValue"]  =acc_max.get(STAT::MEAN);
			WorstCase[method]["MedianValue"]=acc_max.get(STAT::MEDIAN);
			MedianVal[method]["MeanValue"]  =acc_med.get(STAT::MEAN);
			MedianVal[method]["MedianValue"]=acc_med.get(STAT::MEDIAN);
		}
	}
	
	int warning1 = elog.get_token("variants calculated a missing value for BayesDel.");
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldSym(false,true);	// field numb for GeneSymbol
	field_numbers	FldFun(false,true);	// field numb for FuncConseq
	field_numbers	FldDet(false,true);	// field numb for FuncDetail
	field_numbers	FldRes(false,true);	// field numb for result
	field_numbers	FldXAF(false,true);	// field numb for MaxAF
	int	ColDel = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
	int AnnGen=0; // 0: no vAnnGene in INFO; 1: vAnnGene in INFO
	int num_MaxAF_INFO = 0;
	format.set_option(SKIP_NOTES,false);
	bool row1_not_read = true;
	int by_model_0 = elog.get_token("variants calcualted by model 0");
	int by_model_1 = elog.get_token("variants calcualted by model 1");
	int by_model_2 = elog.get_token("variants calcualted by model 2");
	int by_PROVEAN = elog.get_token("variants calcualted by PROVEAN");
	int by_fMKL_NC = elog.get_token("variants calcualted by fathmmMKL_NC");
	int by_CADD = elog.get_token("variants calcualted by CADD");
	int by_worst = elog.get_token("variants calcualted by the higest score among all possible missense variants");
	int by_median = elog.get_token("variants calcualted by the median score among all possible missense variants");
	for (Rows_in_File(in, inputs, &format))
	{
		// read comments and header
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			for (auto &ds:ds_dta) if (str_startsw(in[0],"##INFO=<ID="+ds.first+",")) ds.second.column=-1;
			for (auto &ds:m1_dta) if (str_startsw(in[0],"##INFO=<ID="+ds.first+",")) ds.second.column=-1;
			for (auto &ds:m2_dta) if (str_startsw(in[0],"##INFO=<ID="+ds.first+",")) ds.second.column=-1;
			if (str_startsw(in[0],"##INFO=<ID="+perch::i_func+",")) AnnGen=1;
			if (str_startsw(in[0],"##INFO=<ID=fthmMKL_NC,"))	fMKL_NC.column=-1;
			if (str_startsw(in[0],"##INFO=<ID=PROVEAN,"))		PROVEAN.column=-1;
			if (str_startsw(in[0],"##INFO=<ID=CADD,"))			CADD.column=-1;
			if (str_startsw(in[0],"##INFO=<ID="+wrWhat+"_"+wt_set+",")) { ColDel=-1; AddInf=true; }
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_MxAF+","))	++num_MaxAF_INFO;
			else { for (auto &id : perch::h_afID) if (str_startsw(in[0],"##INFO=<ID="+id+","))	++num_MaxAF_INFO; }
			if (str_startsw(in[0],"##INFO=<ID=Included_MaxAF_in_"+wrWhat+",") && add_af) exit_error("--add-af is applied twice");
			print_container(in.contents(),program.outf,' ',true);
			in.clear_nf();
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			if (AddInf && ColDel!=-1)
				program.outf<<"##INFO=<ID="+wrWhat+"_"+wt_set+",Number=1,Type=Float,Description=\"BayesDel computed by vDEL.\">"<<endl;
			if (add_af || !str_has(wt_set,"_noAF"))
				program.outf<<"##INFO=<ID=Included_MaxAF_in_"+wrWhat+",Number=0,Type=Flag,Description=\"BayesDel computed with --add-af.\">"<<endl;
			row1_not_read=false;
			format.clear_field_nums();
			FldChr.clear();
			FldPos.clear();
			FldRef.clear();
			FldAlt.clear();
			FldInf.clear();
			FldSym.clear();
			FldFun.clear();
			FldDet.clear();
			FldRes.clear();
			FldXAF.clear();
			for (int i=0;i<in.NumFields();++i)
			{
				if (in[i]=="INFO")	FldInf.push_back(i+1);
				if (in[i]==perch::h_MxAF)	FldXAF.push_back(i+1); else { for (auto &id : perch::h_afID) if (in[i]==id) FldXAF.push_back(i+1); }
				if (exist_element(ds_dta,in[i]))	{ if (ds_dta[in[i]].column==-2) ds_dta[in[i]].column=i; else if (ds_dta[in[i]].column>=0) exit_error("duplicated columns "+in[i]); }
				if (exist_element(m1_dta,in[i]))	{ if (m1_dta[in[i]].column==-2) m1_dta[in[i]].column=i; else if (m1_dta[in[i]].column>=0) exit_error("duplicated columns "+in[i]); }
				if (exist_element(m2_dta,in[i]))	{ if (m2_dta[in[i]].column==-2) m2_dta[in[i]].column=i; else if (m2_dta[in[i]].column>=0) exit_error("duplicated columns "+in[i]); }
				if (in[i]=="CADD")					{ if (CADD.column==-2) CADD.column=i; else if (CADD.column>=0) exit_error("duplicated columns "+in[i]); }
				if (in[i]=="PROVEAN")				{ if (PROVEAN.column==-2) PROVEAN.column=i; else if (PROVEAN.column>=0) exit_error("duplicated columns "+in[i]); }
				if (in[i]=="fthmMKL_NC")			{ if (fMKL_NC.column==-2) fMKL_NC.column=i; else if (fMKL_NC.column>=0) exit_error("duplicated columns "+in[i]); }
				if (in[i]==wrWhat+"_"+wt_set)		{ if (ColDel==-2) { ColDel=i; AddInf=false; } else exit_error("multiple annotations for "+wrWhat+"_"+wt_set); }
				if (in[i]==perch::h_symb)	FldSym.push_back(i+1);
				if (in[i]==perch::h_func)	FldFun.push_back(i+1);
				if (in[i]==perch::h_fdet)	FldDet.push_back(i+1);
				string in_i_ = boost::to_upper_copy(in[i]);
				if (exist_element(perch::h_col1,in_i_))												{ FldChr.clear(); FldChr.push_back(i+1); }
				if (in_i_=="START" && FldPos.no_input()) FldPos.push_back(i+1); if (in_i_=="POS")	{ FldPos.clear(); FldPos.push_back(i+1); }
				if (in_i_=="REF")	{ FldRef.clear(); FldRef.push_back(i+1); }
				if (in_i_=="ALT")	{ FldAlt.clear(); FldAlt.push_back(i+1); }
			}
			if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
			if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
			if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
			if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
			if (AddInf && FldInf.no_input()) exit_error("The INFO column is missing.");
			if (add_af && num_MaxAF_INFO==0 && FldXAF.no_input()) exit_error("MaxAF annotation is required by --add-af but is missing");
			if (!WstFor.empty() || !MedFor.empty())
			{
				if ((FldSym.no_input()||FldFun.no_input()||FldDet.no_input()) && !AnnGen)
					exit_error("Functional consequence annotation is missing.");
			}
			if (ColDel==-2)
			{
				for (auto &ds:ds_dta) if (ds.second.column==-2) exit_error(ds.first+" annotation is missing.");
				for (auto &ds:m1_dta) if (ds.second.column==-2) exit_error(ds.first+" annotation is missing.");
				for (auto &ds:m2_dta) if (ds.second.column==-2) exit_error(ds.first+" annotation is missing.");
			}
			if (FldSym.no_input()) NoSort=true;
			if (wr_res)
			{
				if (AddInf)	FldRes=FldInf;
				else
				{
					if (ColDel>=0) { FldRes.push_back(ColDel+1); }
					else { FldRes.push_back(in.NumFields()+1); in.contents().push_back(wrWhat+"_"+wt_set); format.set_field_nums(FldRes,"",tfile_format::Expand); }
				}
			}
			print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
		if (row1_not_read) exit_error("Header line missing.");

		int	chr_num;	if (!genepi::read_chr_num(in[FldChr[0]],chr_num))	exit_error("Failed to read "+in[FldChr[0]]+" as a chromosome.");
		int	bp;			if (!read_val_ge(in[FldPos[0]],bp,1))				exit_error("Failed to read "+in[FldPos[0]]+" as a 1-based bp.");
		string& ref = in[FldRef[0]];
		string& alt = in[FldAlt[0]];
		bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
		string this_index = in[FldChr[0]]+'\t'+in[FldPos[0]]+'\t'+in[FldRef[0]]+'\t'+in[FldAlt[0]];
		
		vector<string> INFO;
		if (!FldInf.no_input()) { if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";")); }

		// get variant type. ann_ws and ann_ms should be exclusive.
		string GeneSymbol;
		string FuncConseq;
		string FuncDetail;
		if (!FldSym.no_input() && !FldFun.no_input() && !FldDet.no_input())
		{
			GeneSymbol=in[FldSym[0]]; if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
			if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
			FuncConseq=to_lower_copy(in[FldFun[0]]);
			if (PRVloc) FuncDetail=in[FldDet[0]];
		}
		else if (AnnGen)
		{
			string func_ann = get_string(INFO,perch::i_func);
			if (!func_ann.empty()) // exit_error(perch::i_func+" not in INFO or empty"); // removed to allow intergenic variants
			{
				vector<string> func_vec;
				boost::split(func_vec,func_ann,boost::is_any_of(","));
				if (func_vec.size()<3 || (func_vec.size()%3)) exit_error(perch::i_func+" vector size wrong");
				GeneSymbol=func_vec[0]; if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
				if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
				FuncConseq=to_lower_copy(func_vec[1]);
				if (PRVloc) FuncDetail=func_vec[2];
			}
		}
		string Func_tx_id;
		string HGVS_p_dot;
		if (!FuncDetail.empty())
		{
			string HGVS_nom=FuncDetail;
			if (str_has(HGVS_nom,"[")) HGVS_nom=substr_before_find(HGVS_nom,"[");
			vector<string> detail_vec;
			boost::split(detail_vec,HGVS_nom,boost::is_any_of(":"));
			if (detail_vec.size()==3 && str_startsw(detail_vec[2],"p."))
			{
				Func_tx_id=detail_vec[0];
				HGVS_p_dot=detail_vec[2].substr(2);
			}
		}
		bool for_ann_by_provean = (!is_snv && (str_has(FuncConseq,"missense")||str_has(FuncConseq,"inframe")));
		//bool for_ensemble_Bayes = ( is_snv && (str_has(FuncConseq,"missense")||str_has(FuncConseq,"stop")||str_has(FuncConseq,"splice")||str_has(FuncConseq,"translat")) );
		
		bool NonMissense_NonSynon = false; for (auto &x:WstFor) if (str_has(FuncConseq,x)) NonMissense_NonSynon=true; // complete loss-of-function variant
		bool MissenseSubstitution = false; for (auto &x:MedFor) if (str_has(FuncConseq,x)) MissenseSubstitution=true; // partial loss-of-function variant
		if (NonMissense_NonSynon) MissenseSubstitution=false;
		bool ann_ws = (NonMissense_NonSynon && !WorstCase.empty());
		bool ann_ms = (MissenseSubstitution || (MedFor.size()==1 && MedFor[0]=="others")) && !MedianVal.empty();
		
		// MaxAF
		double MaxAF = get_value(INFO,perch::h_MxAF);
		for (size_t i=0;i<FldXAF.size();++i)
		{
			double x = std::numeric_limits<double>::signaling_NaN();
			read_val(in[FldXAF[i]],x);
			if ( !std::isnan(x) && (std::isnan(MaxAF) || fabs(x-0.5)<fabs(MaxAF-0.5)) ) MaxAF=x;
		}
		for (auto &id : perch::h_afID)
		{
			double x = get_value(INFO,id);
			if ( !std::isnan(x) && (std::isnan(MaxAF) || fabs(x-0.5)<fabs(MaxAF-0.5)) ) MaxAF=x;
		}
		if (af0toM) { if (MaxAF==0) MaxAF=std::numeric_limits<double>::signaling_NaN(); }
		else		{ if (std::isnan(MaxAF)) MaxAF=0; }
		double MaxAF_lod = 0;
		if (add_af && AF_prp && !std::isnan(MaxAF)) MaxAF_lod = log10(MaxAF_curve.solve(MaxAF));

		// get PROVEAN score in log10(BF) scale
		double PROVlbf = std::numeric_limits<double>::signaling_NaN();
		if (PROVEAN.column>=-2)
		{
			double ReadScore = std::numeric_limits<double>::signaling_NaN();
			if		(PROVEAN.column>=0)		ReadStr(in[PROVEAN.column],ReadScore,0);
			else if (PROVEAN.column==-1)	ReadScore = get_value(INFO,"PROVEAN");
			if (!std::isnan(ReadScore))	PROVlbf = log10(PROVEAN.curve.solve(ReadScore));
			else if (write_missing_del && for_ann_by_provean)
			{
				if (PRVloc && !Func_tx_id.empty() && !HGVS_p_dot.empty())
				{
					ms_PROV << in[FldChr[0]] << ',' << in[FldPos[0]] << ',' << ref << ',' << alt << '\t' << Func_tx_id << '\t' << HGVS_p_dot << endl;
				}
				else
					ms_PROV << in[FldChr[0]] << ',' << in[FldPos[0]] << ',' << ref << ',' << alt << endl;
			}
		}
		
		// cal combined DELlbf is isnan
		double DELlbf = std::numeric_limits<double>::signaling_NaN();
		if (ColDel!=-2)
		{
			if		(ColDel==-1)	DELlbf = get_value(INFO,wrWhat+"_"+wt_set);
			else if (ColDel>=0)	read_val(in[ColDel],DELlbf);
		}
		if (std::isnan(DELlbf))
		{
			double fMNClbf = std::numeric_limits<double>::signaling_NaN();
			double CADDlbf = std::numeric_limits<double>::signaling_NaN();
			if (CADD.column>=-1)
			{
				double ReadScore = std::numeric_limits<double>::signaling_NaN();
				if		(CADD.column>=0)	ReadStr(in[CADD.column],ReadScore,0);
				else if (CADD.column==-1)	ReadScore = get_value(INFO,"CADD");
				if (!std::isnan(ReadScore)) CADDlbf = log10(CADD.curve.solve(ReadScore));
				// else if (write_missing_del && !is_snv) ms_CADD << in[FldChr[0]] << DLMTR << in[FldPos[0]] << DLMTR << '.' << DLMTR << ref << DLMTR << alt << endl;
			}
			if (fMKL_NC.column>=-1)
			{
				double ReadScore = std::numeric_limits<double>::signaling_NaN();
				if		(fMKL_NC.column>=0)		ReadStr(in[fMKL_NC.column],ReadScore,0);
				else if (fMKL_NC.column==-1)	ReadScore = get_value(INFO,"fthmMKL_NC");
				if (!std::isnan(ReadScore))	fMNClbf = log10(fMKL_NC.curve.solve(ReadScore));
				// else if (write_missing_del && is_snv)	ms_fMKL << in[FldChr[0]] << ',' << in[FldPos[0]] << ',' << ref << ',' << alt << endl;
			}
			
			if (for_ann_by_provean && !std::isnan(PROVlbf))	// for coding in-frame indels or missense indels
			{
				if (MaxAF_lod)	DELlbf = AF_prp*MaxAF_lod + (1-AF_prp)*PROVlbf;
				else			DELlbf = PROVlbf;
				elog.add(by_PROVEAN);
			}
			else
			{
				map<string, double>& m0_wts = (std::isnan(MaxAF) ? m0_wtB : m0_wtA);
				map<string, double>& m1_wts = (std::isnan(MaxAF) ? m1_wtB : m1_wtA);
				map<string, double>& m2_wts = (std::isnan(MaxAF) ? m2_wtB : m2_wtA);
				int num_ms = 0;				// number of missing components
				double sum_wt = 0;			// sum of weights of the observed predictors
				map<string,double> ds_bf;	// deleteriousness score _ Bayes factor
				for (auto& ds:ds_dta)
				{
					double r = std::numeric_limits<double>::signaling_NaN();
					if (ds.first==perch::h_MxAF)
						r=MaxAF;
					else
					{
						if		(ds.second.column>=0)	read_val(in[ds.second.column],r);
						else if (ds.second.column==-1)	r=get_value(INFO,ds.first);
					}
					if (!std::isnan(r))
					{
						if (m0_wts[ds.first])
						{
							ds_bf[ds.first] = log10(ds.second.curve.solve(r)); // won't be nan
							sum_wt += m0_wts[ds.first];
							ds_val[ds.first].push_back(r);
						}
					}
					else
					{
						if (ann_ws) r = get_alternative(ds.first,GeneSymbol,WorstCase,r);
						if (ann_ms) r = get_alternative(ds.first,GeneSymbol,MedianVal,r);
						++num_ms;
					}
					if (UpdDel) in[ds.second.column] = ftos(r							   ,MisStr);
					if (ChgDel) in[ds.second.column] = ftos(log10(ds.second.curve.solve(r)),MisStr);
				}
				
				if (sum_wt && sum_wt>m0_cut)
				{
					double wt_sum=0;	// weighted sum
					for (auto& i:ds_bf) wt_sum += m0_wts[i.first] * i.second;
					DELlbf = wt_sum/sum_wt;// for Missense caused by SNV
					elog.add(by_model_0);
				}
				else
				{
					sum_wt = 0;
					ds_bf.clear();
					for (auto& ds:m1_dta)
					{
						double r = std::numeric_limits<double>::signaling_NaN();
						if		(ds.second.column>=0)	read_val(in[ds.second.column],r);
						else if (ds.second.column==-1)	r=get_value(INFO,ds.first);
						if (!std::isnan(r))
						{
							if (m1_wts[ds.first])
							{
								ds_bf[ds.first] = log10(ds.second.curve.solve(r)); // won't be nan
								sum_wt += m1_wts[ds.first];
							}
						}
					}
					if (sum_wt > m12_cut)
					{
						double wt_sum=0;	// weighted sum
						for (auto& i:ds_bf) wt_sum += m1_wts[i.first] * i.second;
						DELlbf = wt_sum/sum_wt;// for Missense caused by SNV
						elog.add(by_model_1);
					}
					else
					{
						sum_wt = 0;
						ds_bf.clear();
						for (auto& ds:m2_dta)
						{
							double r = std::numeric_limits<double>::signaling_NaN();
							if		(ds.second.column>=0)	read_val(in[ds.second.column],r);
							else if (ds.second.column==-1)	r=get_value(INFO,ds.first);
							if (!std::isnan(r))
							{
								if (m2_wts[ds.first])
								{
									ds_bf[ds.first] = log10(ds.second.curve.solve(r)); // won't be nan
									sum_wt += m2_wts[ds.first];
								}
							}
						}
						if (sum_wt)
						{
							double wt_sum=0;	// weighted sum
							for (auto& i:ds_bf) wt_sum += m2_wts[i.first] * i.second;
							DELlbf = wt_sum/sum_wt;// for Missense caused by SNV
							elog.add(by_model_2);
						}
					}
				}																									// for coding SNVs
			}
			if (std::isnan(DELlbf) && !std::isnan(CADDlbf) && !is_snv)	{ DELlbf=CADDlbf; elog.add(by_CADD);	 }	// for non-coding InDels
			if (std::isnan(DELlbf) && !std::isnan(fMNClbf) &&  is_snv)	{ DELlbf=fMNClbf; elog.add(by_fMKL_NC); }	// for non-coding SNVs
			if (std::isnan(DELlbf) && !std::isnan(fMNClbf) && !is_snv)	{ DELlbf=fMNClbf; elog.add(by_fMKL_NC); }	// for others
			if (std::isnan(DELlbf) && !std::isnan(CADDlbf) &&  is_snv)	{ DELlbf=CADDlbf; elog.add(by_CADD);	 }	// for others
		}
		
		// modify BayesDel
		if (for_ann_by_provean && !std::isnan(PROVlbf))
		{
			if (MaxAF_lod)	DELlbf = AF_prp*MaxAF_lod+(1-AF_prp)*PROVlbf;
			else			DELlbf = PROVlbf;
			elog.add(by_PROVEAN);
		}
		else if	(ann_ws)					  {	elog.add(by_worst);  DELlbf = get_alternative(wrWhat,GeneSymbol,WorstCase,DELlbf); if (MaxAF_lod) DELlbf=AF_prp*MaxAF_lod+(1-AF_prp)*DELlbf; }
		else if (ann_ms && std::isnan(DELlbf)){	elog.add(by_median); DELlbf = get_alternative(wrWhat,GeneSymbol,MedianVal,DELlbf); if (MaxAF_lod) DELlbf=AF_prp*MaxAF_lod+(1-AF_prp)*DELlbf; }
		else if (MaxAF_lod) { DELlbf=AF_prp*MaxAF_lod+(1-AF_prp)*DELlbf; }
		
		// write output
		if (std::isnan(DELlbf)) elog.add(warning1);
		string result_str = ftos(DELlbf,MisStr);
		if (wr_res)
		{
			if (AddInf)
			{
				bool exist_BDEL=false;
				for (auto &i:INFO)
					if (str_startsw(i,wrWhat+"_"+wt_set+"=")) { i=wrWhat+"_"+wt_set+"="+result_str; exist_BDEL=true; }
				if (!exist_BDEL) INFO.push_back(wrWhat+"_"+wt_set+"="+result_str);
				in[FldInf[0]]=str_of_container(INFO,';');
			}
			else
				in[FldRes[0]]=result_str;
		}
		if (NoSort)
		{
			print_container(in.contents(),program.outf,DLMTR,true);
		}
		else
		{
			double MAF = (MaxAF<=0.5 ? MaxAF : 1-MaxAF) ;
			add_line_by_gene(str_of_container(in.contents(),DLMTR),DELlbf,MAF,(GeneSymbol.empty()?this_index:GeneSymbol));
		}
	}
	if (!NoSort) add_line_by_gene("",-INFINITY,1,"");
	if (write_missing_del)
	{
		closefile(ms_PROV);
		// closefile(ms_CADD);
		// closefile(ms_fMKL);
	}
	
	// check DEL. compressed: MA FATHMM. Flipped: SIFT (20%<0.5) MT (40%<0.5) FATHMM (92%<0.5) LRT (11%<0.5)
	if (ds_val["MA"].get(STAT::N)>1000 && str_startsw(wt_set,"nsfp23") && !UpdDel && !ChgDel)
	{
		if (ds_val["MA"]    .get(STAT::MIN)<0 || ds_val["MA"]    .get(STAT::MAX)>1) exit_error("The MA score is not converted.");
		if (ds_val["FATHMM"].get(STAT::MIN)<0 || ds_val["FATHMM"].get(STAT::MAX)>1) exit_error("The FATHMM score is not converted.");
		if (ds_val["FATHMM"].num_le(0.5)/ds_val["FATHMM"].get(STAT::N) <0.5) lns<<showw<<"<50% of the variants have a FATHMM score <=0.5. FATHMM might not be converted." << flush_logger;
		if (ds_val["SIFT"]  .num_le(0.5)/ds_val["SIFT"]  .get(STAT::N) >0.5) lns<<showw<<">50% of the variants have a SIFT   score <=0.5. SIFT might not be converted." << flush_logger;
		if (ds_val["MT"]    .num_le(0.5)/ds_val["MT"]    .get(STAT::N) >0.5) lns<<showw<<">50% of the variants have a MT     score <=0.5. MT might not be converted." << flush_logger;
		if (ds_val["LRT"]   .num_le(0.5)/ds_val["LRT"]   .get(STAT::N) >0.5) lns<<showw<<">50% of the variants have a LRT    score <=0.5. LRT might not be converted." << flush_logger;
	} //*/
	return 0;
}

/*
To debug, 2.tmp and 3.tmp should have a same results
 sed 's/Class1/0/;s/Class2/0/;s/Class3/./;s/Class4/1/;s/Class5/1/' b1b2_ann.xls | sed 1d | MAKE_BF2 | df.cut -f 37,6,8,10,16,18,20,22,24,26,28,30 > 1.tmp
 NaiveBayes each 1.tmp --weights 0.0389957,0.00897645,0.0379296,0.0801522,0.0655452,0.0880422,0.20379,0.174036,0.170067,0.0404486,0.0920167 > 2.tmp
 sed 's/Class1/0/;s/Class2/0/;s/Class3/./;s/Class4/1/;s/Class5/1/' b1b2_ann.xls | df.insert --titles INFO --strg '' | vDEL > 3.tmp
*/
