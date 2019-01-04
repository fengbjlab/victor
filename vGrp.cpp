#include <boost/range/adaptor/reversed.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_genepi.hpp>
#include "victor_par.hpp"

using namespace std;
using boost::iostreams::filtering_ostream;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	vector<string>	inputs;
	vector<string>	inGSet;
	string			inCllp;
	set<string>		excl_G;
	string			excl_F;
	set<string>		incl_G;
	string			incl_F;
	set<string>		invo_G;
	string			invo_F;
	bool 			epacts=false;
	bool 			SNonly=false;	// analysis restricted to SNVs
	bool			toRest=true;
//	int				FltMAC=0;		// removed FltMAC because it's hard to decide whether the cutoff makes sense without knowing the TotSpl
	string			InfoAC="AC";	// INFO field for AC
	string			InfoAN="AN";	// INFO field for AN
	string			h_dels;			// header of BayesDel defined by user
	int				minSet=2;
	int				maxSet=50;
	string			incGBA;
	string			invGBA;
	double			gbaCut=0.302426; // 0.302426 is the median of known genes by IC2. So the power to include a 2-gene set with at least 1 known gene is 75%.

	// handle program options
	perch::MisCut=1;
	perch::VQSsnv=-INFINITY;
	perch::VQSidl=-INFINITY;
	perch::filflt.clear();
	program.enable_option("--nt");
	program.enable_option("--prefix");
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--epacts"))			ReadArg(program.arg(),argi,epacts);
		else if (str_startsw(program.arg()[argi],"--snv-only"))			ReadArg(program.arg(),argi,SNonly);
//		else if (str_startsw(program.arg()[argi],"--filt-mac"))			ReadArg(program.arg(),argi,FltMAC);
		else if (str_startsw(program.arg()[argi],"--info-ac"))			ReadArg(program.arg(),argi,InfoAC);
		else if (str_startsw(program.arg()[argi],"--info-an"))			ReadArg(program.arg(),argi,InfoAN);
		else if (str_startsw(program.arg()[argi],"--gene-set"))			ReadSet(program.arg(),argi,inGSet);
		else if (str_startsw(program.arg()[argi],"--collapse"))			ReadArg(program.arg(),argi,inCllp);
		else if (str_startsw(program.arg()[argi],"--min-set"))			ReadArg(program.arg(),argi,minSet);
		else if (str_startsw(program.arg()[argi],"--max-set"))			ReadArg(program.arg(),argi,maxSet);
		else if (str_startsw(program.arg()[argi],"--excl-genes"))		ReadSet(program.arg(),argi,excl_G);
		else if (str_startsw(program.arg()[argi],"--excl-gene-file"))	ReadArg(program.arg(),argi,excl_F);
		else if (str_startsw(program.arg()[argi],"--incl-genes"))		ReadSet(program.arg(),argi,incl_G);
		else if (str_startsw(program.arg()[argi],"--incl-gene-file"))	ReadArg(program.arg(),argi,incl_F);
		else if (str_startsw(program.arg()[argi],"--invo-genes"))		ReadSet(program.arg(),argi,invo_G);
		else if (str_startsw(program.arg()[argi],"--invo-gene-file"))	ReadArg(program.arg(),argi,invo_F);
		else if	(str_startsw(program.arg()[argi],"--restore"))			ReadArg(program.arg(),argi,toRest);
		else if (str_startsw(program.arg()[argi],"--incl-gba-file"))	ReadArg(program.arg(),argi,incGBA);
		else if (str_startsw(program.arg()[argi],"--invo-gba-file"))	ReadArg(program.arg(),argi,invGBA);
		else if (str_startsw(program.arg()[argi],"--invo-gba-cut"))		ReadArg(program.arg(),argi,gbaCut);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else add_to_container(inputs,program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_epacts",str_YesOrNo(epacts));
	program.help_text_var("_Default_snv_only",str_YesOrNo(SNonly));
//	program.help_text_var("_Default_filt_mac",itos(FltMAC));
	program.help_text_var("_Default_info_ac",InfoAC);
	program.help_text_var("_Default_info_an",InfoAN);
	program.help_text_var("_Default_gene_set",str_of_container(inGSet,',',false));
	program.help_text_var("_Default_min_set",itos(minSet));
	program.help_text_var("_Default_max_set",itos(maxSet));
	program.help_text_var("_Default_excl_genes",str_of_container(excl_G,',',false));
	program.help_text_var("_Default_excl_gene_file",excl_F);
	program.help_text_var("_Default_incl_genes",str_of_container(incl_G,',',false));
	program.help_text_var("_Default_incl_gene_file",incl_F);
	program.help_text_var("_Default_invo_genes",str_of_container(invo_G,',',false));
	program.help_text_var("_Default_invo_gene_file",invo_F);
	program.help_text_var("_Default_incl_gba_file",incGBA);
	program.help_text_var("_Default_invo_gba_file",invGBA);
	program.help_text_var("_Default_collapse",inCllp);
	program.help_text_var("_Default_gba_cut",ftos(gbaCut));
	program.help_text_var("_Default_restore",str_YesOrNo(toRest));
	perch::check_arguments();
	
	// check errors
	for (auto &f:inGSet) f=perch::find_file(f);
	if (!inCllp.empty()) inCllp=perch::find_file(inCllp);
	if (!inGSet.empty() && inCllp.empty()) exit_error("--gene-set requires --collapse to be set");
	if (!excl_F.empty()) for (Rows_in_File(in,excl_F,1)) excl_G.insert(in[0]);
	if (!incl_F.empty()) for (Rows_in_File(in,incl_F,1)) incl_G.insert(in[0]);
	if (!invo_F.empty()) for (Rows_in_File(in,invo_F,1)) invo_G.insert(in[0]);
	if (!incGBA.empty()) for (Rows_in_File(in,incGBA,2)) { if (as_double_or_nan(in[1])>=gbaCut) incl_G.insert(in[0]); }
	if (!invGBA.empty()) for (Rows_in_File(in,invGBA,2)) { if (as_double_or_nan(in[1])>=gbaCut) invo_G.insert(in[0]); }
	if (program.prefix().empty()) exit_error("--prefix is not set");
	
	// read gene list
	set<string> GwPV; // genes with a gene-wise p-value, which means there're variants in data
	if (!inCllp.empty())
	{
		tfile_format inCllp_format;
		inCllp_format.set_delimiters("\t");
		inCllp_format.comment_sw()="##";
		int gene_column=-1;
		for (Rows_in_File(in,inCllp,&inCllp_format))
		{
			bool is_header=false;
			for (int i=0;i<in.NumFields();++i)
				if (in[i]==perch::h_symb) { gene_column=i; is_header=true; break; }
			if (is_header) continue;
			if (gene_column==-1) exit_error("column "+perch::h_symb+" not found in "+inCllp);
			bool is_skipped=false;
			for (int i=0;i<in.NumFields();++i)
				if (str_has(in[i],"SKIPPED")&&i!=gene_column) { is_skipped=true; }
			if (is_skipped) continue;
			string GeneSymbol=in[gene_column];
			if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
			if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"("))
			{
				GeneSymbol=substr_after_find(GeneSymbol,"(");
				GeneSymbol.pop_back();
			}
			GwPV.insert(GeneSymbol);
		}
		if (GwPV.empty()) exit_error("nothing read from "+inCllp);
		lns << showl << "Read " << GwPV.size() << " genes from " << inCllp << flush_logger;
	}
	
	// read variant list
	typedef tuple<int,int,string,string> 				loc_t;	// chr,pos,ref,alt
	typedef tuple<double,double,int,int,string,string> 	del_t;	// del,xaf,chr,pos,ref,alt
	map<string, map< loc_t,string > > 				data_by_loc;// data_by_loc[gene][loc]=line
	map<string, map< del_t,string > > 				data_by_del;// data_by_loc[gene][del]=line
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldXAF(false,true);	// field numb for MaxAF
	field_numbers	FldSym(false,true);	// field numb for "Gene Symbol"
	field_numbers	FldFun(false,true);	// field numb for "Function Type"
	field_numbers	FldDet(false,true);	// field numb for FuncDetail
	field_numbers 	FldOut(false,true); // field numb for output
	FldOut.push_back(1,5);
	int	ColDel = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
	tfile_format format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);
	vector<string> header;
	int err_chr=elog.get_token("variants removed due to wrong #CHROM.");
	int err_pos=elog.get_token("variants removed due to wrong POS.");
	int fil_reg=elog.get_token("variants filtered out by region.");
	int fil_xaf=elog.get_token("variants filtered out by MaxAF.");
	int fil_saf=elog.get_token("variants filtered out by SplAF.");
	int fil_faf=elog.get_token("variants filtered out by FdrAF.");
	int fil_del=elog.get_token("variants filtered out by BayesDel.");
	int fil_mis=elog.get_token("variants filtered out by missing rate.");
	int fil_snv=elog.get_token("variants filtered out by type (SNV).");
	int fil_lof=elog.get_token("variants filtered out by LoF.");
	int fil_DNg=elog.get_token("variants filtered out by DomNeg.");
	int fil_cds=elog.get_token("variants filtered out by CDS.");
	int var_pas=elog.get_token("variants pass filtering.");
	stringstream output_header;
	for (Rows_in_File(in, inputs, &format))
	{
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#') // meta lines
		{
			if (in.FileNumber()==0)
			{
				perch::read_meta(in[0]);
				if ( h_dels.empty() && str_startsw(in[0],"##INFO=<ID=BayesDel")) ColDel=-1;
				if (!h_dels.empty() && str_startsw(in[0],"##INFO=<ID="+h_dels+",")) ColDel=-1;
				if (!epacts) print_container(in.contents(),output_header,' ',true);
			}
			in.clear_nf();
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			lns << showl << "Read " << in.FileName()<< flush_logger;
			if (in.FileNumber()==0)
			{
				header=in.contents();
				
				// get columns
				FldChr.clear();
				FldPos.clear();
				FldRef.clear();
				FldAlt.clear();
				FldInf.clear();
				FldXAF.clear();
				FldSym.clear();
				FldFun.clear();
				FldDet.clear();
				for (int i=0;i<in.NumFields();++i)
				{
					if (in[i]=="#CHROM")		{	if (FldChr.no_input()) { FldChr.push_back(i+1); 					   } else exit_error("multiple columns for #CHROM");}
					if (in[i]=="POS")			{	if (FldPos.no_input()) { FldPos.push_back(i+1); 					   } else exit_error("multiple columns for POS"); 	}
					if (in[i]=="REF")			{	if (FldRef.no_input()) { FldRef.push_back(i+1); 					   } else exit_error("multiple columns for REF"); 	}
					if (in[i]=="ALT")			{	if (FldAlt.no_input()) { FldAlt.push_back(i+1); 					   } else exit_error("multiple columns for ALT"); 	}
					if (in[i]=="INFO")			{	if (FldInf.no_input()) { FldInf.push_back(i+1); FldOut.push_back(i+1); } else exit_error("multiple columns for INFO"); 	}
					if (in[i]==perch::h_MxAF)	{	if (FldXAF.no_input()) { FldXAF.push_back(i+1); FldOut.push_back(i+1); } else exit_error("multiple columns for MaxAF"); }
					if (in[i]==perch::h_symb)	{	if (FldSym.no_input()) { FldSym.push_back(i+1); FldOut.push_back(i+1); } else exit_error("multiple columns for Gene Symbol");}
					if (in[i]==perch::h_func)	{	if (FldFun.no_input()) { FldFun.push_back(i+1); FldOut.push_back(i+1); } else exit_error("multiple columns for Consequence Type");}
					if (in[i]==perch::h_fdet)	{	if (FldDet.no_input()) { FldDet.push_back(i+1); FldOut.push_back(i+1); } else exit_error("multiple columns for Consequence Detail");}
					if ( h_dels.empty() && str_startsw(in[i],"BayesDel"))	{ if (ColDel==-2) { ColDel=i; FldOut.push_back(i+1); } else exit_error("multiple columns for BayesDel"); }
					if (!h_dels.empty() && in[i]==h_dels)					{ if (ColDel==-2) { ColDel=i; FldOut.push_back(i+1); } else exit_error("multiple columns for "+h_dels); }
				}
				if (FldChr.no_input()) exit_error("The #CHROM column is missing.");
				if (FldPos.no_input()) exit_error("The POS column is missing.");
				if (FldRef.no_input()) exit_error("The REF column is missing.");
				if (FldAlt.no_input()) exit_error("The ALT column is missing.");
				if (FldSym.no_input()) exit_error("The Gene Symbol column is missing.");
				if (FldFun.no_input()) exit_error("The Functional Consequence column is missing.");

				// write header
				if (!epacts)
				{
					if (!inGSet.empty()) {	in.write_r(output_header,FldOut,false); output_header<<DLMTR<<"GeneSet"<<endl; }
					else					in.write_r(output_header,FldOut,true);
				}
			}
			else
			{
				if (header != in.contents()) exit_error("input VCF files have different headers");
			}
			continue;
		}
		
		// basic info
		int chr_num = genepi::read_chr_num(in[FldChr[0]]); if (chr_num<1) { elog.add(err_chr); continue; }
		int bp=-1; if (!read_val_gt(in[FldPos[0]],bp,0)) { elog.add(err_pos); continue; }
		string& ref = in[FldRef[0]];
		string& alt = in[FldAlt[0]];
		bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
		bool is_lof = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"LoF")||str_has(in[FldFun[0]],"NMD") );
		bool is_vks = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"knClinSig=1") );
		bool is_DNg = ( FldFun.no_input() ? false : perch::is_DomNeg(in[FldFun[0]]) );
		bool is_cds = ( FldFun.no_input() ? false : perch::is_coding(in[FldFun[0]]) );

		// get INFO
		vector<string> INFO;
		if (!FldInf.no_input())
			if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
		
		string GeneSymbol=in[FldSym[0]];
		if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
		if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }

		// skip by chr region
		if (!perch::within_covered_region(chr_num,bp)) { elog.add(fil_reg); continue; }

		// skip by variant type
		if (SNonly && !is_snv) { elog.add(fil_snv); continue; }
		if (perch::LFonly && is_lof && perch::is_LoFtol(GeneSymbol)) { elog.add(fil_lof); continue; }
		if (perch::LFonly && !is_lof) { elog.add(fil_lof); continue; }
		if (perch::CDonly && !is_cds) { elog.add(fil_cds); continue; }
		if (perch::DomNeg && !is_DNg) { elog.add(fil_DNg); continue; }

		// skip by MaxAF
		double MaxAF = std::numeric_limits<double>::signaling_NaN();
		if (!FldXAF.no_input())	read_val(in[FldXAF[0]],MaxAF);
		else					MaxAF=get_value(INFO,perch::h_MxAF);
		if (perch::filXAF)
		{
			double f = ( std::isnan(MaxAF) ? 0 : (MaxAF>0.5 ? 1-MaxAF : MaxAF) );
			if (!perch::rf_XAF && f> perch::filXAF) { elog.add(fil_xaf); continue; }
			if ( perch::rf_XAF && f<=perch::filXAF) { elog.add(fil_xaf); continue; }
		}

		// skip by observed allele frequency
		if (perch::filSAF)
		{
			double SplAF = get_value(INFO,"SplAF");
			if (std::isnan(SplAF)) exit_error("cannot apply --filt-SplAF because SplAF is not calculated by vQC");
			if (SplAF>0.5) SplAF=1-SplAF;
			if (!perch::rf_SAF && SplAF> perch::filSAF) { elog.add(fil_saf); continue; }
			if ( perch::rf_SAF && SplAF<=perch::filSAF) { elog.add(fil_saf); continue; }
		}
		
		// skip by founder allele frequency
		if (perch::filFAF)
		{
			double FdrAF = get_value(INFO,"AF_Founder");
			double FdrAN = get_value(INFO,"AN_Founder");
			if (!std::isnan(FdrAF) && !std::isnan(FdrAN))
			{
				if (FdrAF>0.5) FdrAF=1-FdrAF;
				if (FdrAN>=perch::filFAFminAN)
				{
					if (!perch::rf_FAF && FdrAF> perch::filFAF) { elog.add(fil_faf); continue; }
					if ( perch::rf_FAF && FdrAF<=perch::filFAF) { elog.add(fil_faf); continue; }
				}
			}
		}

		/*/ skip by MAC
		if (FltMAC)
		{
			int INFO_AC=-1; get_int(INFO,InfoAC,INFO_AC);
			int INFO_AN=-1; get_int(INFO,InfoAN,INFO_AN);
			if (INFO_AC<0) exit_error(InfoAC+" sub-field missing or not a valid number (non-negative integer)");
			if (INFO_AN<0) exit_error(InfoAN+" sub-field missing or not a valid number (non-negative integer)");
			int mac = std::min(INFO_AC, INFO_AN-INFO_AC);
			if (mac<FltMAC) { elog.add(fil_mac); continue; }
		} //*/

		// skip by BayesDel
		double BayesDel = std::numeric_limits<double>::signaling_NaN();
		if		(ColDel==-1)	BayesDel = get_value_sw(INFO,"BayesDel");
		else if (ColDel>=0)		read_val(in[ColDel],BayesDel);
		if (perch::filter_AnnAF(BayesDel,GeneSymbol,is_lof,is_vks)) { elog.add(fil_del); continue; }

		// skip by missing rate
		if (perch::MisCut!=1)
		{
			if (perch::Mis_ea)
			{
				double MssCS=get_value(INFO,"MissingCs");
				double MssCT=get_value(INFO,"MissingCt");
				if (MssCS>perch::MisCut || MssCT>perch::MisCut) { elog.add(fil_mis); continue; }
			}
			else
			{
				double MsgAll=get_value(INFO,"MissingRate");
				if (MsgAll>perch::MisCut) { elog.add(fil_mis); continue; }
			}
		}

		// restore original POS,REF,ALT
		if (toRest)
		{
			for (auto& f:boost::adaptors::reverse(INFO))
				if (str_startsw(f,"OriginalIndex="))
				{
					string original_index=f.substr(14);
					vector<string> components;
					boost::split(components,original_index,boost::is_any_of("_"));
					if (components.size()==4)
					{
						in[FldChr[0]]=components[0];
						in[FldPos[0]]=components[1]; if (!read_val_gt(in[FldPos[0]],bp,0)) exit_error("failed to read basepair from "+in[FldPos[0]]);
						in[FldRef[0]]=components[2];
						in[FldAlt[0]]=components[3];
					}
					break;
				}
		}

		// add data
		string line;
		if (epacts)	line=in[FldChr[0]]+":"+in[FldPos[0]]+"_"+in[FldRef[0]]+"/"+in[FldAlt[0]];
		else		FldOut.contents_to_a_string(in.contents(),line,DLMTR);
		loc_t variant = make_tuple(chr_num,bp,ref,alt);
		data_by_loc[GeneSymbol][variant]=line;
		elog.add(var_pas);
	}
	
	// write file
	int output_ptr=0;
	vector<std::shared_ptr<filtering_ostream> >	grpOut;
	for (int i=1;i<=program.nt;++i)
	{
		grpOut.push_back(make_shared<filtering_ostream>());
		if ( !openOutFile(*grpOut.back(),program.prefix()+"."+itos(i)) ) exit_cannotOpen(program.prefix()+"."+itos(i));
		*grpOut.back() << output_header.str();
	}
	if (!inGSet.empty())
	{
		set< set<string> > all_sets;
		int num_gene_sets = elog.get_token("gene sets written");
		int num_lt_min	=   elog.get_token("gene sets excluded because number of genes < minimum allowed");
		int num_gt_max	=   elog.get_token("gene sets excluded because number of genes > maximum allowed");
		int num_dup		=   elog.get_token("gene sets excluded due to duplication");
		int num_not_inv	=   elog.get_token("gene sets excluded due to --invo-xxx");
		tfile_format gs_format;
		gs_format.set_delimiters("\t");
		gs_format.forbid_nf_rpt();
		for (Rows_in_File(gs,inGSet,&gs_format))
		{
			string GeneSetName=gs[0];
			set<string> genes;
			int involved=0;
			for (int i=2;i<gs.NumFields();++i)
				if (exist_element(data_by_loc,gs[i]) && exist_element(GwPV,gs[i]) && !exist_element(excl_G,gs[i]) && (incl_G.empty()||exist_element(incl_G,gs[i])))
				{
					genes.insert(gs[i]);
					if (invo_G.empty()||exist_element(invo_G,gs[i])) involved++;
				}
			if ((int)genes.size()<minSet) { elog.add(num_lt_min); continue; }
			if ((int)genes.size()>maxSet) { elog.add(num_gt_max); continue; }
			if (exist_element(all_sets,genes)) { elog.add(num_dup); continue; }
			if (!involved) { elog.add(num_not_inv); continue; }
			all_sets.insert(genes);
			map< loc_t,string > this_set;
			for (int i=2;i<gs.NumFields();++i)
			{
				if (exist_element(data_by_loc,gs[i]) && !exist_element(excl_G,gs[i]) && (incl_G.empty()||exist_element(incl_G,gs[i])))
				{
					for (auto &l:data_by_loc[gs[i]])
						this_set[l.first]=l.second;
				}
			}
			if (epacts)
			{
				*grpOut[output_ptr]<<GeneSetName;
				for (auto &l:this_set) *grpOut[output_ptr]<<DLMTR<<l.second;
				*grpOut[output_ptr]<<endl;
			}
			else
			{
				for (auto &l:this_set) *grpOut[output_ptr]<<l.second<<DLMTR<<GeneSetName<<endl;
			}
			elog.add(num_gene_sets);
			++output_ptr; if (output_ptr>=program.nt) output_ptr=0;
		}
	}
	else
	{
		int num_genes_wr = elog.get_token("genes written");
		for (auto &g:data_by_loc)
		{
			if (!excl_G.empty() &&  exist_element(excl_G,g.first)) continue;
			if (!incl_G.empty() && !exist_element(incl_G,g.first)) continue;
			if (epacts)
			{
				*grpOut[output_ptr]<<g.first;
				for (auto &l:g.second) *grpOut[output_ptr]<<DLMTR<<l.second;
				*grpOut[output_ptr]<<endl;
			}
			else
			{
				for (auto &l:g.second) *grpOut[output_ptr]<<l.second<<endl;
			}
			elog.add(num_genes_wr);
			++output_ptr; if (output_ptr>=program.nt) output_ptr=0;
		}
	}
	for (auto & out : grpOut) closefile(*out);

	return 0;
}
