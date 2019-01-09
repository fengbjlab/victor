/*
 Removed options:
 --groups STRs     In input, header of group columns {_Default_groups}
 --no-zero B       Exclude the groups that have zero scores for all components except biological relevance and BayesDel {_Default_no_zero}
 --norm STAT/BOOL  Normalize each component before adding up {no. If yes, default stat is MED_ABS}
 --except-biol B   Do not normalize biological relevant scores {No}
 */
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_genepi.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	vector<string>	inputs;					// input files
	vector<string>	aLBF_f;					// filenames for additional log10 Bayes factors
	set<string>		aLBF_h={"gene","trait","gene\\trait","trait\\gene","symbol\\trait","trait\\symbol","symbol","gene_symbol","genesymbol"}; // header of column 1 of aLBF_f
	bool			AddSEG=true;			// add component co-segregation
	bool			AddHLR=true;			// add component haplotype likelihood ratio
	bool			AddGLR=true;			// add component composite likelihood ratio
	bool			AddGBA=true;			// add component gene relevance
	bool			HasSEG=false;			// has h_SEGb
	bool			HasHLR=false;			// has h_HLRb
	bool			HasGLR=false;			// has h_GLRb
	bool			HasMLP=false;			// has h_mLgP
	bool			HasPCR=false;			// has h_Pcrt
	bool			Has_OR=false;			// has h_o_r_
	bool			HasGBA=false;			// has gene relevance
	bool			NegGBA=true;			// add negative GBA score
	bool			NoZero=true;			// Some variants have 0 HLR & SEG but high GBA & DEL. They may rank high but should not be counted.
	bool			No0seg=false;			// exclude variant that vSEG log10 Bayes factor=0
	bool			No0hlr=false;			// exclude variant that vAAA log10 Bayes factor=0
	bool			NoNcGn=false;			// exclude ncRNA genes.
	bool			NoNrmB=false;			// no normalizing GBA
	double			NoRVIS=0;				// no top % of genes with high RVIS scores (tolerance). BRCA2 is 99.59%, no wonder it show up in my UGP. 171 genes 99%, 849 95%, 1703 90%.
	double			priorp=0.5;				// prior probability for variant classification, 0.5 is non-informative
	string			norm_t;					// normalization type, empty=no_normalization
	string			MisStr="NA";			// missing value
	string			biol_f;					// GeneSymbol->GBAlbf file
	string			biol_t;					// GBA trait name, if biol_f has multiple traits then find this trait
	string			h_grNM;					// header of the groupID name (col1:col2:col3)
	string			h_dels;					// header of BayesDel defined by user
	string			plot_t;					// manhattan plot
	string			plot_f;					// --plot-to file
	set<string>		h_grID={perch::h_symb};	// header of the groupID column (col1 col2 col3)
	string			DGD_in = "dgd_Hsa_all_v71.tsv";
	string			RVISfn = "RVIS";
	bool			do_nothing = false;
	bool			outsym = false;			// for GSA, output each gene instead of gene set name
	string			pergen;
	vector<string>	inGSet;
	set<string>		det_fn = { perch::h_Adet };

	// handle program options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(program.arg()[argi]=="--no-seg")	AddSEG=false;
		else if (program.arg()[argi]=="--no-glr")	AddGLR=false;
		else if (program.arg()[argi]=="--no-hlr")	AddHLR=false;
		else if (program.arg()[argi]=="--no-biol")	AddGBA=false;
		else if (str_startsw(program.arg()[argi],"--no-zero"))		ReadArg(program.arg(),argi,NoZero);
		else if (str_startsw(program.arg()[argi],"--no-0seg"))		ReadArg(program.arg(),argi,No0seg);
		else if (str_startsw(program.arg()[argi],"--no-0hlr"))		ReadArg(program.arg(),argi,No0hlr);
		else if (str_startsw(program.arg()[argi],"--add-seg"))		ReadArg(program.arg(),argi,AddSEG);
		else if (str_startsw(program.arg()[argi],"--add-glr"))		ReadArg(program.arg(),argi,AddGLR);
		else if (str_startsw(program.arg()[argi],"--add-hlr"))		ReadArg(program.arg(),argi,AddHLR);
		else if (str_startsw(program.arg()[argi],"--add-biol"))		ReadArg(program.arg(),argi,AddGBA);
		else if (str_startsw(program.arg()[argi],"--neg-biol"))		ReadArg(program.arg(),argi,NegGBA);
		else if (str_startsw(program.arg()[argi],"--dgd"))			ReadSet(program.arg(),argi,DGD_in);
		else if (str_startsw(program.arg()[argi],"--lbf"))			ReadSet(program.arg(),argi,aLBF_f);
		else if (str_startsw(program.arg()[argi],"--norm"))			ReadArg(program.arg(),argi,norm_t);
		else if (str_startsw(program.arg()[argi],"--plot-to"))		ReadArg(program.arg(),argi,plot_f);
		else if	(str_startsw(program.arg()[argi],"--group"))		ReadSet(program.arg(),argi,h_grID);
		else if (str_startsw(program.arg()[argi],"--biol"))			ReadArg(program.arg(),argi,biol_f);
		else if (str_startsw(program.arg()[argi],"--trait"))		ReadArg(program.arg(),argi,biol_t);
		else if (str_startsw(program.arg()[argi],"--prior"))		ReadArg(program.arg(),argi,priorp);
		else if (str_startsw(program.arg()[argi],"--no-ncRNA"))		ReadArg(program.arg(),argi,NoNcGn);
		else if (str_startsw(program.arg()[argi],"--no-rvis"))		ReadArg(program.arg(),argi,NoRVIS);
		else if (str_startsw(program.arg()[argi],"--do-nothing"))	ReadArg(program.arg(),argi,do_nothing);
		else if (str_startsw(program.arg()[argi],"--per-gene"))		ReadArg(program.arg(),argi,pergen);
		else if (str_startsw(program.arg()[argi],"--gene-set"))		ReadSet(program.arg(),argi,inGSet);
		else if (str_startsw(program.arg()[argi],"--out-fn"))		ReadSet(program.arg(),argi,det_fn);
		else if (str_startsw(program.arg()[argi],"--show-genes"))	ReadArg(program.arg(),argi,outsym);
		//else if (str_startsw(program.arg()[argi],"--except-biol"))	ReadArg(program.arg(),argi,NoNrmB);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else add_to_container(inputs,program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_groups",str_of_container(h_grID,',',false));
	program.help_text_var("_Default_lbf",str_of_container(aLBF_f,',',false));
	program.help_text_var("_Default_plot_to",plot_f);
	program.help_text_var("_Default_biol",biol_f);
	program.help_text_var("_Default_dgd",DGD_in);
	program.help_text_var("_Default_rvis",ftos(NoRVIS));
	program.help_text_var("_Default_neg_biol",str_YesOrNo(NegGBA));
	program.help_text_var("_Default_add_seg",str_YesOrNo(AddSEG));
	program.help_text_var("_Default_add_glr",str_YesOrNo(AddGLR));
	program.help_text_var("_Default_add_hlr",str_YesOrNo(AddHLR));
	program.help_text_var("_Default_add_biol",str_YesOrNo(AddGBA));
	program.help_text_var("_Default_no_zero",str_YesOrNo(NoZero));
	program.help_text_var("_Default_no_0seg",str_YesOrNo(No0seg));
	program.help_text_var("_Default_no_0hlr",str_YesOrNo(No0hlr));
	program.help_text_var("_Default_no_ncRNA",str_YesOrNo(NoNcGn));
	program.help_text_var("_Default_except_biol",str_YesOrNo(NoNrmB));
	program.help_text_var("_Default_do_nothing",str_YesOrNo(do_nothing));
	program.help_text_var("_Default_out_gene",str_YesOrNo(outsym));
	program.help_text_var("_Default_out_fn",str_of_container(det_fn,',',false));
	program.help_text_var("_Default_per_gene",pergen);
	program.help_text_var("_Default_gene_set",str_of_container(inGSet,',',false));
	perch::check_arguments();

	if (do_nothing)
	{
		tfile_format format;
		format.set_delimiters("\t");
		format.set_option(SKIP_NOTES,false);
		for (Lines_in_File(in, inputs, &format)) program.outf << in[0] << endl;
		return 0;
	}

	// check errors
	if (!plot_f.empty())
	{
		plot_t = substr_after_rfind(plot_f,".");
		if (plot_t!="pdf") exit_error(plot_t+" is not a supported figure file format. Please use pdf.");
		string gnuplot_version = exec(perch::gnuplot+" --version | cut -d' ' -f 2",true);
		if (gnuplot_version.empty()) exit_error("cannot execute the gnuplot program.");
		int n1 = extract_int (gnuplot_version);
		extract_char(gnuplot_version);
		int n2 = extract_int (gnuplot_version);
		if (n1<4 || (n1==4 && n2<6)) exit_error("gnuplot version ("+exec(perch::gnuplot+" --version | xargs echo -n",true)+") is too old.");
	}
	if (perch::VarCla)
	{
		perch::h_FINb="PostPrb";
		lns<<showl<<"The --group option will be ignored for --vc."<<flush_logger;
	}
	for (auto &f:inGSet) f=perch::find_file(f);
	
	// about normalization
	// STAT use MED_POS is better than MAD because
	// 1) it's more meaningful to compute deviation from 0 rather than median;
	// 2) use positive values only will focus on the front instead of the whole, which may improve partial AUC
	STAT::StatType	n_stat = STAT::MED_POS;	// divided by statistics
	bool			toNorm = false;			// normalize before adding. It's good for full AUC but bad for partial AUC. The latter is more important.
	if (perch::VarCla)
	{
		toNorm=false;
	}
	else
	{
		try
		{
			toNorm=IsYes(norm_t);
		}
		catch (...)
		{
			if (!norm_t.empty()) { toNorm=true; n_stat=to_StatType(norm_t); }
			else toNorm=false;
		}
	}
	
	// read gene db
	set<string> ncRNA;
	if (NoNcGn)
	{
		genepi::read_genes();
		for (auto &g:genepi::gene_byTranscript) for (auto &t:g.second) if (t.second.cdsStart==t.second.cdsEnd) ncRNA.insert(t.second.name2);
	}

	// read GBAlbf
	map<string,double> SymbToGBA; // SymbToGBA[symbol]=GBAlbf
	if (!biol_f.empty())
	{
		int num_columns=0; for (Rows_in_File(in,biol_f,1)) { num_columns=in.NumFields(); break; }
		if (num_columns>2)
		{
			vector<string> GeneSymbols;
			for (Rows_in_File(in,biol_f,1))
			{
				if (in.RowNumber()==0) {	GeneSymbols=in.contents(); continue; }
				if (!biol_t.empty() && in[0]!=biol_t) continue;
				for (int i=1;i<in.NumFields();++i)
				{
					double v = 0;
					try { v=boost::lexical_cast<double>(in[i]); }
					catch (...) { exit_error("Error reading "+biol_f+": failed to read "+in[i]+" as a number."); }
					if (!std::isnan(v)) SymbToGBA[GeneSymbols[i]] = v;
				}
				break;
			}
		}
		else if (num_columns==2)
		{
			for (Rows_in_File(in,biol_f,2))
			{
				if (exist_element(aLBF_h,to_lower_copy(in[0]))) continue;
				double v = 0;
				try { v=boost::lexical_cast<double>(in[1]); }
				catch (...) { exit_error("Error reading "+biol_f+": failed to read "+in[1]+" as a number."); }
				if (!std::isnan(v)) SymbToGBA[in[0]] = v;
			}
		}
		else exit_error(biol_f+" has the wrong number of columns.");
		if (SymbToGBA.empty()) exit_error("Nothing was read from "+biol_f);
		HasGBA=true;
		AddGBA=true;
	}
	else
	{
		HasGBA=false;
		AddGBA=false;
	}

	// read DGD
	if (!DGD_in.empty())
	{
		//int summary1=elog.get_token("duplicated genes with unknown GBA scores have been filled in.");
		map<string, vector<string> > dgd; // dgd[group]=vector<gene_symbol>
		tfile_format DGDfmt;
		DGDfmt.set_delimiters("\t");
		DGDfmt.set_titlelines(1);
		string fn = perch::find_file(DGD_in);
		for (Rows_in_File(in,fn,&DGDfmt))
			dgd[in[1]].push_back(in[7]);
		for (auto &gene_group:dgd)
		{
			Values<double> gba;
			vector<string> unk;
			for (auto &gene:gene_group.second)
				if (exist_element(SymbToGBA,gene))	gba.push_back(SymbToGBA[gene]);
				else								unk.push_back(gene);
			if (gba.empty() || unk.empty()) continue;
			for (auto &gene:unk)
			{
				SymbToGBA[gene] = gba.get(STAT::MEAN);
				//elog.add(summary1,gene);
			}
		}
	}
	
	// read RVIS
	map<string,double> RVISdb;
	for (Rows_in_File(in,perch::find_file(RVISfn),5))
	{
		if (in[0]=="GENE") continue;
		double v; if (!read_val(in[4],v)) continue; // NA
		RVISdb[in[0]]=v;
	}
	
	// read additional LBF score
	vector< map<string,double> > aLBF_s; // aLBF_s[#][group_name_with_colon_as_field_delimiter] = additional LBF score
	#define X_it(j) size_t j=0;j<aLBF_s.size();++j
	for (auto &f:aLBF_f)
	{
		aLBF_s.push_back(map<string,double>());
		map<string,double>& this_map = aLBF_s.back();
		for (Rows_in_File(in,f,2))
		{
			if (exist_element(aLBF_h,to_lower_copy(in[0]))) continue;
			double v = 0;
			try { v=boost::lexical_cast<double>(in[1]); }
			catch (...) { exit_error("Error reading "+f+": failed to read "+in[1]+" as a number."); }
			if (!std::isnan(v)) this_map[in[0]] = v;
		}
	}
	
	// read per gene score file and gene set DB. They are used for gene set analysis.
	map<string,double> pg_score;
	if (!pergen.empty())
	{
		for (Rows_in_File(in,pergen,4))
		{
			if (in[0]==perch::h_symb) continue;
			string GeneSymbol=in[0];
			if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
			if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
			try { pg_score[GeneSymbol] = extract_double(in[3]); } catch (...) {}
		}
	}
	map<string, set<string> > all_sets;
	if (!inGSet.empty())
	{
		tfile_format gs_format;
		gs_format.set_delimiters("\t");
		gs_format.forbid_nf_rpt();
		for (Rows_in_File(gs,inGSet,&gs_format))
		{
			string& GeneSetName=gs[1];
			set<string> genes;
			for (int i=2;i<gs.NumFields();++i) genes.insert(gs[i]);
			all_sets[GeneSetName]=genes;
		}
	}
	
	// data. groups[#], vec...[#], lbf...[#] and ...Xtr[][#] have the same coordinate #. Val...[#] is for statistics, so the coordinate doesn't matter.
	set<string> groups;
	map<string,double> lbfDEL;
	map<string,double> lbfSEG;
	map<string,double> lbfHLR;
	map<string,double> lbfGLR;
	map<string,double> lbfGBA;
	map<string,double> resMLP;
	map<string,double> resPCR;
	map<string,double> res_OR;
	vector< map<string,double> > l_X(aLBF_s.size());
	map<string,string> resXtr;
	
	// info about genes
	map<string, string>				geneChr; // Gene Symbol => chromosome string
	map<string, Values<long> >		genePos; // Gene Symbol => position in basepair
	map<string, Values<double> >	geneFin; // Gene Symbol => final score
	
	// info about gropus
	map<string, vector<string> >	grpSym;
	map<string, vector<string> >	grpChr;
	map<string, Values<long> >		grpPos;
	map<string, Values<double> >	grpGBAval;
	map<string, set<string> >		grpGBASym;
	map<string, set<string> >		grp2gen; // group => gene symbol
	map<string, string >			grp2fun; // group => functional consequence, for VarCla only

	// read
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldSym(false,true);	// field numb for GeneSymbol
	field_numbers	FldFun(false,true);	// field numb for FuncConseq
	field_numbers	FldDet(false,true);	// field numb for FuncDetail
	field_numbers	FldGrp(false,true);	// field numb for group ID
	field_numbers	FldXAF(false,true);	// field numb for MaxAF
	int	ColSEG=-2;	// column number for vSEG analysis results as log10 Bayes factors
	int	ColHLR=-2;	// column number for vAAA analysis results as log10 Bayes factors (default)
	int	ColGLR=-2;	// column number for vAAA analysis results as log10 Bayes factors (--glr)
	int	ColMLP=-2;	// column number for vAAA analysis results as minus log1- p (--out-mlp)
	int	ColPCR=-2;	// column number for vAAA analysis results as corrected p (--permute)
	int	Col_OR=-2;	// column number for vAAA analysis results as odds ratios (--out-or)
	int	ColDel=-2;	// column number for BayesDel. -2 = not annotated; -1 = in INFO; 0+ = column.
	int	ColXtr=-2;	// column number for Xtr. -2 = not annotated; -1 = in INFO; 0+ = column.
	tfile_format	format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);
	string		prev_group;	// previous group ID
	set<string> old_groups; // to check whether rows are sorted by groups
	int AnnGen=0; // 0: no vAnnGene in INFO; 1: vAnnGene in INFO
	bool	row1_not_read = true;
	bool	collapse_out = false;
	bool	details_out = false;
	set<string>	collapse_head = { "Gene","No.Var","cs_neg","cs_het","cs_hom" };
	field_numbers details_fld(false,true);
	int warning1=elog.get_token("duplicated groups");
	for (Rows_in_File(in, inputs, &format))
	{
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			in.clear_nf();
			perch::read_meta(in[0]);
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_SEGb+",")) ColSEG=-1;
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_HLRb+",")) ColHLR=-1;
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_GLRb+",")) ColGLR=-1;
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_mLgP+",")) ColMLP=-1;
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_Pcrt+",")) ColPCR=-1;
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_o_r_+",")) Col_OR=-1;
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_Axtr+",")) ColXtr=-1;
			if ( h_dels.empty() && str_startsw(in[0],"##INFO=<ID=BayesDel")) ColDel=-1;
			if (!h_dels.empty() && str_startsw(in[0],"##INFO=<ID="+h_dels+",")) ColDel=-1;
			if (str_has(in[0],"##INFO=<ID="+perch::i_func+",")) AnnGen=1;
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			row1_not_read=false;
			if (exist_any(collapse_head, in.contents()))
			{
				collapse_out=true;
				print_container(in.contents(),program.outf,DLMTR,true);
				continue;
			}
			format.clear_field_nums();
			FldChr.clear();
			FldPos.clear();
			FldRef.clear();
			FldAlt.clear();
			FldInf.clear();
			FldSym.clear();
			FldFun.clear();
			FldDet.clear();
			FldGrp.clear();
			FldXAF.clear();
			for (int i=0;i<in.NumFields();++i)
			{
				if (in[i]==perch::h_SEGb)	{ ColSEG=i; }
				if (in[i]==perch::h_HLRb)	{ ColHLR=i; }
				if (in[i]==perch::h_GLRb)	{ ColGLR=i; }
				if (in[i]==perch::h_mLgP)	{ ColMLP=i; }
				if (in[i]==perch::h_Pcrt)	{ ColPCR=i; }
				if (in[i]==perch::h_o_r_)	{ Col_OR=i; }
				if (in[i]==perch::h_Axtr)	{ ColXtr=i; }
				if (in[i]=="INFO")					FldInf.push_back(i+1);
				if (in[i]==perch::h_symb)			FldSym.push_back(i+1);
				if (in[i]==perch::h_func)			FldFun.push_back(i+1);
				if (in[i]==perch::h_fdet)			FldDet.push_back(i+1);
				if (in[i]==perch::h_MxAF)			FldXAF.push_back(i+1);
				if (exist_element(h_grID,in[i]))	FldGrp.push_back(i+1);
				if ( h_dels.empty() && str_startsw(in[i],"BayesDel"))	{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for BayesDel"); }
				if (!h_dels.empty() && in[i]==h_dels)					{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for "+h_dels); }
				if (in[i]=="Chr"   && FldChr.no_input())	FldChr.push_back(i+1); if (in[i]=="#CHROM")		{	FldChr.clear(); FldChr.push_back(i+1); }
				if (in[i]=="Start" && FldPos.no_input())	FldPos.push_back(i+1); if (in[i]=="POS")		{	FldPos.clear(); FldPos.push_back(i+1); }
				if (in[i]=="Ref"   && FldRef.no_input())	FldRef.push_back(i+1); if (in[i]=="REF")		{	FldRef.clear(); FldRef.push_back(i+1); }
				if (in[i]=="Alt"   && FldAlt.no_input())	FldAlt.push_back(i+1); if (in[i]=="ALT")		{	FldAlt.clear(); FldAlt.push_back(i+1); }
				if (exist_element(det_fn,in[i])) { details_out=true; details_fld.push_back(i+1); }
			}
			if (ColSEG!=-2) HasSEG=true;
			if (ColHLR!=-2) HasHLR=true;
			if (ColGLR!=-2) HasGLR=true;
			if (ColMLP!=-2 && !perch::VarCla) HasMLP=true;
			if (ColPCR!=-2 && !perch::VarCla) HasPCR=true;
			if (Col_OR!=-2 && !perch::VarCla) Has_OR=true;
			if (HasMLP && HasHLR) exit_error("contain both "+perch::h_mLgP+" and "+perch::h_HLRb+", don't know which one to use.");
			if (HasMLP && HasPCR) exit_error("contain both "+perch::h_mLgP+" and "+perch::h_Pcrt+", don't know which one to use.");
			if (HasMLP && Has_OR) exit_error("contain both "+perch::h_mLgP+" and "+perch::h_o_r_+", don't know which one to use.");
			if (HasPCR && Has_OR) exit_error("contain both "+perch::h_Pcrt+" and "+perch::h_o_r_+", don't know which one to use.");
			if (FldChr.no_input())						exit_error("The #CHROM/Chr column is missing.");
			if (FldPos.no_input())						exit_error("The POS/Start column is missing.");
			if (FldSym.no_input() && !AnnGen)			exit_error("The "+perch::h_symb+" annotation is missing.");
			if (perch::VarCla && FldRef.no_input())		exit_error("The REF column is missing.");
			if (perch::VarCla && FldAlt.no_input())		exit_error("The ALT column is missing.");
			if (perch::VarCla) FldGrp = FldChr + FldPos + FldRef + FldAlt + FldSym;
			if (FldGrp.no_input())						exit_error("The groupID column is missing.");
			FldGrp.contents_to_a_string(in.contents(),h_grNM,':');
			if (details_out)
			{
				program.outf<<"#CHROM:POS:REF:ALT"<<DLMTR<<"GBAlbf"<<DLMTR<<perch::h_MxAF<<DLMTR<<(h_dels.empty()?"BayesDel":h_dels)<<DLMTR<<perch::h_symb<<DLMTR<<perch::h_func<<DLMTR<<perch::h_fdet<<DLMTR;
				in.write_r(program.outf,details_fld,true);
				continue;
			}
			continue;
		}
		if (row1_not_read) exit_error("Header lines missing.");
		
		// write
		if (collapse_out)
		{
			print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}

		vector<string> INFO;
		if (!FldInf.no_input()) { if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";")); }

		string GeneSymbol;
		string FuncConseq;
		string FuncDetail;
		if (AnnGen)
		{
			string func_ann = get_string(INFO,perch::i_func);
			if (func_ann.empty()) exit_error(perch::i_func+" not in INFO or empty");
			vector<string> func_vec;
			boost::split(func_vec,func_ann,boost::is_any_of(","));
			if (func_vec.size()<3) exit_error(perch::i_func+" vector size wrong");
			GeneSymbol=func_vec[0];
			FuncConseq=func_vec[1];
			FuncDetail=func_vec[2];
		}
		if (!FldSym.no_input()) GeneSymbol=in[FldSym[0]];
		if (!FldFun.no_input()) FuncConseq=in[FldFun[0]];
		if (!FldDet.no_input()) FuncDetail=in[FldDet[0]];
		if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
		if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
		
		int bp = -1; try { bp = boost::lexical_cast<int>(in[FldPos[0]]); } catch (...) { exit_error("Failed to read "+in[FldPos[0]]+" as a position in basepairs."); }
		if (NoNcGn && exist_element(ncRNA,GeneSymbol)) continue;
		if (NoRVIS && exist_element(RVISdb,GeneSymbol) && RVISdb[GeneSymbol]>=NoRVIS) continue;
		geneChr[GeneSymbol] = in[FldChr[0]];
		genePos[GeneSymbol].push_back(bp);
		
		// group
		string current_group;
		FldGrp.contents_to_a_string(in.contents(),current_group,':');
		grp2gen[current_group].insert(GeneSymbol);
		grp2fun[current_group]=FuncConseq;
		if (current_group==prev_group)
		{
			// if (perch::VarCla) exit_error("For variant classification, grouping column(s) should uniquely identify each variant.");
			// removed this line because I also use this for variant ranking, and one variant could be related to two overlapping genes.
		}
		else
		{
			if (exist_element(old_groups,current_group))
			{
				if (perch::VarCla) ; // exit_error("For variant classification, grouping column(s) should uniquely identify each variant."); // same reason as above
				else elog.add(warning1,current_group); // exit_error("The Genotype File is not sorted by groups.");
			}
			old_groups.insert(current_group);
			prev_group=current_group;
		}
		
		double BayesDel = std::numeric_limits<double>::signaling_NaN();
		if		(ColDel==-1)	BayesDel = get_value_sw(INFO,"BayesDel");
		else if (ColDel>=0)	read_val(in[ColDel],BayesDel);

		double MaxAF = std::numeric_limits<double>::signaling_NaN();
		if (!FldXAF.no_input())	read_val(in[FldXAF[0]],MaxAF);
		else					MaxAF=get_value(INFO,perch::h_MxAF);
		if (std::isnan(MaxAF)) MaxAF=0;
		
		// write
		if (details_out)
		{
			program.outf<<in[FldChr[0]]<<":"<<in[FldPos[0]]<<":"<<in[FldRef[0]]<<":"<<in[FldAlt[0]]<<DLMTR;
			if (exist_element(SymbToGBA,GeneSymbol)) program.outf<<SymbToGBA[GeneSymbol]; else program.outf<<".";
			program.outf<<DLMTR<<MaxAF<<DLMTR<<BayesDel<<DLMTR<<GeneSymbol<<DLMTR<<FuncConseq<<DLMTR<<FuncDetail<<DLMTR;
			in.write_r(program.outf,details_fld,true);
			continue;
		}

		// add data, robust to nan or missing
		groups.insert(current_group);
		if (grpSym[current_group].empty() || grpSym[current_group].back()!=GeneSymbol)
		{
			grpSym[current_group].push_back(GeneSymbol);
			grpChr[current_group].push_back(in[FldChr[0]]);
		}
		grpPos[current_group].push_back(bp);
		if (HasGBA)
		{
			if (!exist_element(grpGBASym[current_group],GeneSymbol) && exist_element(SymbToGBA,GeneSymbol))
			{
				double v=SymbToGBA[GeneSymbol];
				grpGBASym[current_group].insert(GeneSymbol);
				grpGBAval[current_group].push_back(v);
			}
			if (grpGBAval[current_group].empty())	lbfGBA[current_group]=0;
			else 									lbfGBA[current_group]=grpGBAval[current_group].get(STAT::MAX);
		}
		for (X_it(j)){double v=aLBF_s[j][current_group]; double v0=l_X[j][current_group]; if (!std::isnan(v)&&(v>v0||v0==0)) l_X[j][current_group]=v;} // higest score among this group (domain/gene/pathway)
		//if (HasGBA) { double v=SymbToGBA[GeneSymbol];    double v0=lbfGBA[current_group]; if (!std::isnan(v)&&(v>v0||v0==0)) lbfGBA[current_group]=v;} // higest score among this group (domain/gene/pathway)
		if (perch::VarCla) { double v=BayesDel;	if (!std::isnan(v)) lbfDEL[current_group]=v; }
		if (HasSEG) { double v=0; if (perch::read_variable(in,ColSEG,INFO,perch::h_SEGb,v))	lbfSEG[current_group]=v; }
		if (HasHLR) { double v=0; if (perch::read_variable(in,ColHLR,INFO,perch::h_HLRb,v))	lbfHLR[current_group]=v; }
		if (HasGLR) { double v=0; if (perch::read_variable(in,ColGLR,INFO,perch::h_GLRb,v))	lbfGLR[current_group]=v; }
		if (HasMLP) { double v=0; if (perch::read_variable(in,ColMLP,INFO,perch::h_mLgP,v))	resMLP[current_group]=v; }
		if (HasPCR) { double v=0; if (perch::read_variable(in,ColPCR,INFO,perch::h_Pcrt,v))	resPCR[current_group]=v; }
		if (Has_OR) { double v=0; if (perch::read_variable(in,Col_OR,INFO,perch::h_o_r_,v))	res_OR[current_group]=v; }
		if (ColXtr==-1) { string v=get_string(INFO,perch::h_Axtr); if (!v.empty()) resXtr[current_group]=v; }
		
	}
	if (collapse_out || details_out) return 0;
	
	// normalize. Previously removed normalizing GBA because we don't want to make its weight >1 no matter what.
	double wt_DEL=1, wt_SEG=1, wt_HLR=1, wt_GLR=1, wt_GBA=1;
	vector<double> w_X(aLBF_s.size(),1);
	if (toNorm)
	{
		Values<double> ValGBA;
		Values<double> ValOth;
		for (auto &i:groups)
		{
			double gba=0, oth=0;
			for (X_it(j))		if (exist_element(l_X[j],i)) oth += l_X[j][i];
			if (HasSEG && AddSEG && exist_element(lbfSEG,i)) oth += lbfSEG[i];
			if (HasHLR && AddHLR && exist_element(lbfHLR,i)) oth += lbfHLR[i];
			if (HasGLR && AddGLR && exist_element(lbfGLR,i)) oth += lbfGLR[i];
			if (HasGBA && AddGBA && exist_element(lbfGBA,i)) gba  = lbfGBA[i];
			ValGBA.push_back(gba);
			ValOth.push_back(oth);
		}
		if (!ValGBA.empty() && !ValOth.empty())
		{
			double StaGBA = ValGBA.get(n_stat);
			double StaOth = ValOth.get(n_stat);
			if (!std::isnan(StaGBA) && !std::isnan(StaOth) && StaGBA!=0 && StaOth!=0)
			{
				double num_wt = l_X.size()+1;
				if (HasSEG && AddSEG) ++num_wt;
				if (HasHLR && AddHLR) ++num_wt;
				if (HasGLR && AddGLR) ++num_wt;
				double sum_wt = 1/StaGBA + 1/StaOth;
				double gba_wt = num_wt/StaGBA/sum_wt;
				double oth_wt = (num_wt-gba_wt)/(num_wt-1);
				if (gba_wt < 1)	// The purpose of normalization is to down-weight GBA if sample size is small, so that GBA does not swamp the others.
				{				// But if gba_wt is >=1 (more than its share), then either n_stat is not right or the normalization is not necessary.
					w_X.assign(aLBF_s.size(),oth_wt);
					wt_SEG = oth_wt;
					wt_HLR = oth_wt;
					wt_GLR = oth_wt;
					wt_GBA = gba_wt;
					lns<<showl << "Normalized --biol so that it won't swamp the other components if your sample size is small."<<flush_logger;
				}
				else
				{
					lns<<showl<<"Normalization of --biol was not done because the calculated weight for GBA is >=1."<<flush_logger;
				}
			}
		}
		// cerr<<ValOth.get(n_stat)<<' '<<ValGBA.get(n_stat)<<' '<<wt_HLR<<' '<<wt_GBA<<endl;
		
	/*	vector< Values<double> > v_X(aLBF_s.size());
				for (X_it(j)) { for (auto &x:l_X[j]) v_X[j].push_back(x.second); }
		Values<double> ValDEL;	for (auto &x:lbfDEL) ValDEL.push_back(x.second);
		Values<double> ValSEG;	for (auto &x:lbfSEG) ValSEG.push_back(x.second);
		Values<double> ValHLR;	for (auto &x:lbfHLR) ValHLR.push_back(x.second);
		Values<double> ValGLR;	for (auto &x:lbfGLR) ValGLR.push_back(x.second);
		Values<double> ValGBA;	for (auto &x:lbfGBA) ValGBA.push_back(x.second);
		double num_wt=0;
		double sum_wt=0;
		for (X_it(j))	if (!std::isnan(v_X[j].get(n_stat))) if (v_X[j].get(n_stat)!=0) { sum_wt+=1/v_X[j].get(n_stat); ++num_wt; }
						if (!std::isnan(ValSEG.get(n_stat))) if (ValSEG.get(n_stat)!=0) { sum_wt+=1/ValSEG.get(n_stat); ++num_wt; }
						if (!std::isnan(ValHLR.get(n_stat))) if (ValHLR.get(n_stat)!=0) { sum_wt+=1/ValHLR.get(n_stat); ++num_wt; }
						if (!std::isnan(ValGLR.get(n_stat))) if (ValGLR.get(n_stat)!=0) { sum_wt+=1/ValGLR.get(n_stat); ++num_wt; }
		if (!NoNrmB) {	if (!std::isnan(ValGBA.get(n_stat))) if (ValGBA.get(n_stat)!=0) { sum_wt+=1/ValGBA.get(n_stat); ++num_wt; } }

		for (X_it(j))	if (!std::isnan(v_X[j].get(n_stat))) if (v_X[j].get(n_stat)!=0) w_X[j] = num_wt/v_X[j].get(n_stat)/sum_wt;
						if (!std::isnan(ValSEG.get(n_stat))) if (ValSEG.get(n_stat)!=0) wt_SEG = num_wt/ValSEG.get(n_stat)/sum_wt;
						if (!std::isnan(ValHLR.get(n_stat))) if (ValHLR.get(n_stat)!=0) wt_HLR = num_wt/ValHLR.get(n_stat)/sum_wt;
						if (!std::isnan(ValGLR.get(n_stat))) if (ValGLR.get(n_stat)!=0) wt_GLR = num_wt/ValGLR.get(n_stat)/sum_wt;
		if (!NoNrmB) {	if (!std::isnan(ValGBA.get(n_stat))) if (ValGBA.get(n_stat)!=0) wt_GBA = num_wt/ValGBA.get(n_stat)/sum_wt; }
		// cerr<<ValSEG.get(n_stat)<<' '<<ValHLR.get(n_stat)<<' '<<ValGBA.get(n_stat)<<' '<<wt_SEG<<' '<<wt_HLR<<' '<<wt_GBA<<endl;
	 */
	}
	
	// remove negative GBA. Do it after normalization so that median of GBA is rarely zero.
	if (!NegGBA)
		for (auto &g:lbfGBA)
			if (g.second<0) g.second=0;
	
	// calculate lbfFIN
	map<string,double>	lbfFIN; // lbfFIN[group] = final log10 Bayes factor
	map<string,int>		NmComp; // number of components except GBA
	bool HasFIN;
	if		(HasPCR)
	{
		for (auto &i:groups) if (exist_element(resPCR,i)) lbfFIN[i] = -log10(resPCR[i]);
		AddSEG = false;
		AddHLR = false;
		AddGLR = false;
		AddGBA = false;
		HasFIN = false;
	}
	else if (Has_OR)
	{
		for (auto &i:groups) if (exist_element(res_OR,i)) lbfFIN[i] = res_OR[i];
		AddSEG = false;
		AddHLR = false;
		AddGLR = false;
		AddGBA = false;
		HasFIN = false;
	}
	else if (HasMLP && !HasSEG)
	{
		for (auto &i:groups) if (exist_element(resMLP,i)) lbfFIN[i] = resMLP[i];
		AddSEG = false;
		AddHLR = false;
		AddGLR = false;
		AddGBA = false;
		HasFIN = false;
	}
	else if (HasMLP && HasSEG)
	{
		for (auto &i:groups)
		{
			double SEGpv=0; if (lbfSEG[i]!=0) { SEGpv=genepi::LOD_to_P(lbfSEG[i]); ++NmComp[i]; }
			double AAApv=0; if (exist_element(resMLP,i)) { AAApv=pow(10,-resMLP[i]); ++NmComp[i]; }
			if (SEGpv && AAApv) { double x2 = -2*(log(SEGpv)+log(AAApv)); lbfFIN[i] = -log10(cdf_chisq_q(x2,4)); }
			else if (SEGpv) { lbfFIN[i] = -log10(SEGpv); }
			else if (AAApv) { lbfFIN[i] = -log10(AAApv); }
			else 			{ lbfFIN[i] = std::numeric_limits<double>::signaling_NaN(); }
		}
		HasFIN = true;
		perch::h_FINb="FINmlp";
	}
	else // !HasPCR && !Has_OR && !HasMLP
	{
		if (perch::VarCla)	  for (auto &i:groups) { if (lbfDEL[i]!=0) { lbfFIN[i]+=lbfDEL[i]*wt_DEL;              } }
		if (HasSEG && AddSEG) for (auto &i:groups) { if (lbfSEG[i]!=0) { lbfFIN[i]+=lbfSEG[i]*wt_SEG; ++NmComp[i]; } }
		if (HasHLR && AddHLR) for (auto &i:groups) { if (lbfHLR[i]!=0) { lbfFIN[i]+=lbfHLR[i]*wt_HLR; ++NmComp[i]; } }
		if (HasGLR && AddGLR) for (auto &i:groups) { if (lbfGLR[i]!=0) { lbfFIN[i]+=lbfGLR[i]*wt_GLR; ++NmComp[i]; } }
		if (HasGBA && AddGBA) for (auto &i:groups) { if (lbfGBA[i]!=0) { lbfFIN[i]+=lbfGBA[i]*wt_GBA;              } }
		for (X_it(j))		  for (auto &i:groups) { if (l_X[j][i]!=0) { lbfFIN[i]+=l_X[j][i]*w_X[j]; ++NmComp[i]; } }
		HasFIN = true;
	}

	// output
	program.outf << h_grNM << DLMTR << "#CHROM" << DLMTR << "POS";
	if (HasFIN)				program.outf << DLMTR << perch::h_FINb;
	if (HasMLP)				program.outf << DLMTR << perch::h_mLgP;
	if (HasPCR)				program.outf << DLMTR << perch::h_Pcrt;
	if (Has_OR)				program.outf << DLMTR << perch::h_o_r_;
	if (HasSEG && AddSEG)	program.outf << DLMTR << perch::h_SEGb;
	if (HasHLR && AddHLR)	program.outf << DLMTR << perch::h_HLRb;
	if (HasGLR && AddGLR)	program.outf << DLMTR << perch::h_GLRb;
	if (HasGBA)				program.outf << DLMTR << "GBAlbf";
	if (perch::VarCla)		program.outf << DLMTR << "BayesDel";
	if (perch::VarCla)		program.outf << DLMTR << perch::h_func;
	for (auto &f:aLBF_f)	program.outf << DLMTR << f;
	if (!resXtr.empty())	program.outf << DLMTR << perch::h_Axtr;
	program.outf << endl;
	
	multimap< double,string,std::greater<double> > sorted_output;
	for (auto &i:groups)
	{
		if (NoZero)
		{
			if (HasFIN && NmComp[i]==0) continue;
			if (HasMLP && resMLP[i]==0) continue;
			if (HasPCR && resPCR[i]==0) continue;
			if (Has_OR && res_OR[i]==0) continue;
		}
		if (No0seg && lbfSEG[i]==0) continue;
		if (No0hlr && lbfHLR[i]==0) continue;
		string chr_str="NA", pos_str="NA", grp_str=i;
		if (exist_element(grpSym,i))
		{
			if (grpSym[i].size()==1)
			{
				int mean = (*grpPos[i].begin() + *grpPos[i].rbegin())/2;
				pos_str=itos(mean);
				int chrNum = genepi::read_chr_num(*grpChr[i].begin()); if (!chrNum) exit_error("Can't read chr "+*grpChr[i].begin());
				chr_str = genepi::cytoband(chrNum,*grpPos[i].begin(),*grpPos[i].rbegin());
			}
			else
			{
				if (outsym)
				{
					chr_str = str_of_container(grpChr[i],",");
					grp_str = str_of_container(grpSym[i],",");
					if (!inGSet.empty())
					{
						vector<string> gs_names;
						for (auto &gs:all_sets)
						{
							if (exist_all(gs.second,grpSym[i])) gs_names.push_back(gs.first);
						}
						if (!gs_names.empty())
						{
							grp_str+='(';
							grp_str+=str_of_container(gs_names,";"); // possible delimiters are ; \ | %
							// GeneSet_bi_fam.gmt      name has SPACE # ( ) . / : < _
							// GeneSet_MSigDB_lt50.gmt name has - . / : _
							// GeneSet_hgnc_fam.gmt    name has ' ( ) + , - . / : _
							grp_str+=')';
						}
					}
				}
				if (!pergen.empty() && !perch::VarCla)
				{
					double overall_score = std::numeric_limits<double>::signaling_NaN();
					if (HasFIN)	overall_score=lbfFIN[i];
					if (HasMLP)	overall_score=resMLP[i];
					if (HasPCR)	overall_score=resPCR[i];
					if (Has_OR)	overall_score=res_OR[i];
					int gt_num=0;
					int tot_num=0;
					if (!std::isnan(overall_score))
					{
						for (auto &symbol:grpSym[i])
							if (exist_element(pg_score,symbol))
							{
								tot_num+=1;
								if (overall_score>pg_score[symbol]) ++gt_num;
							}
					}
					double pc_incr = (tot_num ? (double)gt_num/tot_num : std::numeric_limits<double>::signaling_NaN());
					if 		(tot_num==0)	pos_str=".";
					else if (pc_incr==1)	pos_str="+++++:"+itos(gt_num)+"/"+itos(tot_num);
					else if	(pc_incr>=0.8)	pos_str="++++_:"+itos(gt_num)+"/"+itos(tot_num);
					else if	(pc_incr>=0.6)	pos_str="+++__:"+itos(gt_num)+"/"+itos(tot_num);
					else if	(pc_incr>=0.4)	pos_str="++___:"+itos(gt_num)+"/"+itos(tot_num);
					else if	(pc_incr>=0.2)	pos_str="+____:"+itos(gt_num)+"/"+itos(tot_num);
					else 					pos_str="_____:"+itos(gt_num)+"/"+itos(tot_num);
				}
			}
		}
		stringstream ss;
		ss << grp_str << DLMTR << chr_str << DLMTR << pos_str;
		// below previously used ftos_MaxWidthNoTiny() for the sort command to work. now since I sort the output, I nolonger need it.
		if (HasFIN)				ss << DLMTR << ( perch::VarCla ? posterior_given_log10BF(lbfFIN[i],priorp) : lbfFIN[i] );
		if (HasMLP)				ss << DLMTR << ( resMLP[i] );
		if (HasPCR)				ss << DLMTR << ( resPCR[i] );
		if (Has_OR)				ss << DLMTR << ( res_OR[i] );
		if (HasSEG && AddSEG)	ss << DLMTR << ( lbfSEG[i] );
		if (HasHLR && AddHLR)	ss << DLMTR << ( lbfHLR[i] );
		if (HasGLR && AddGLR)	ss << DLMTR << ( lbfGLR[i] );
		if (HasGBA)				ss << DLMTR << ( lbfGBA[i] );
		if (perch::VarCla)		ss << DLMTR << ( lbfDEL[i] );
		if (perch::VarCla)		ss << DLMTR << ( exist_element(grp2fun,i)?grp2fun[i]:MisStr );
		for (X_it(j))			ss << DLMTR << ( l_X[j][i] );
		if (!resXtr.empty())	ss << DLMTR << ( exist_element(resXtr,i)?resXtr[i]:MisStr );
		sorted_output.insert(make_pair(lbfFIN[i],ss.str()));
	}
	for (auto &i:sorted_output)
		program.outf << i.second << endl;
	
	// out manhattan plot of lbfFIN
	if (!plot_f.empty())
	{
		double y_max=0, y_min=0;
		for (auto &i:groups)
		{
			if (NoZero && NmComp[i]==0) continue;
			for (auto &g:grp2gen[i])
			{
				geneFin[g].push_back(lbfFIN[i]);
				if (lbfFIN[i]>y_max) y_max=lbfFIN[i];
				if (lbfFIN[i]<y_min) y_min=lbfFIN[i];
			}
		}
		y_max = ceil (y_max/10)*10;
		y_min = floor(y_min/10)*10;
		y_max = std::max(y_max,-y_min);
		y_min = std::min(-y_max,y_min);
		if (HasPCR||HasMLP) y_min=0;
		
		boost::iostreams::filtering_ostream outf;
		if (!openOutFile( outf,"perch.plot.txt" )) exit_cannotOpen("perch.plot.txt");
		for (auto &g:geneFin)
		{
			int chrNum = genepi::read_chr_num(geneChr[g.first]);
			if (!chrNum) exit_error("Failed to read "+geneChr[g.first]+" as a chromosome.");
			outf << genepi::chrlen_bp_cumulative(chrNum)+genePos[g.first].get(STAT::MEAN) << '\t' << geneFin[g.first].get(STAT::MAX) << '\t' << chrNum%2+2 << endl;
		}
		closefile(outf);
		
		if (!openOutFile( outf,"perch.plot.cmd" )) exit_cannotOpen("perch.plot.cmd");
		outf << "# usage: "<<perch::gnuplot<<" perch.plot.cmd" << endl;
		outf << "set terminal "<<plot_t<<" size 11,8.5" << endl;
		outf << "set output \""<<plot_f<<"\"" << endl;
		outf << "set multiplot" << endl;
		outf << "set bmargin 2" << endl;
		outf << "set tmargin 1.5" << endl;
		outf << "unset xtics" << endl;
		outf << "plot [0:"<<(double)genepi::total_genome_len()<<"] ["<<y_min<<":"<<y_max<<"] \"perch.plot.txt\" using 1:2:3 with points pointtype 7 pointsize 0.5 linecolor variable tit \"\"" << endl;
		closefile(outf);
		
		exec(perch::gnuplot+" perch.plot.cmd",false);
		exec("rm perch.plot.cmd",false);
		exec("rm perch.plot.txt",false);
	}
	
	return 0;
}
