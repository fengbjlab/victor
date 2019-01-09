/*
 Caveat: --detail and --ann-af filters look at one variant at a time. But it's possible that similar or more severe variants exceed the cutoff aggregately.
 However, it's hard to solve this problem because we don't know the LD of the variants. So I leave it as is. 
*/

#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	string			var_in;
	string			det_in;
	string			sym_in;
	string			log_fn;
	int				cut_obs=1;
	double			cut_seg=-INFINITY;
	double			cut_aaa=-INFINITY;
	double			cut_gba=-INFINITY;

	// handle program options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--detail"))		ReadArg(program.arg(),argi,det_in);
		else if	(str_startsw(program.arg()[argi],"--log"))			ReadArg(program.arg(),argi,log_fn);
		else if	(str_startsw(program.arg()[argi],"--gene-file"))	ReadArg(program.arg(),argi,sym_in);
		else if	(str_startsw(program.arg()[argi],"--filt-obs"))		ReadArg(program.arg(),argi,cut_obs);
		else if	(str_startsw(program.arg()[argi],"--filt-seg"))		ReadArg(program.arg(),argi,cut_seg);
		else if	(str_startsw(program.arg()[argi],"--filt-aaa"))		ReadArg(program.arg(),argi,cut_aaa);
		else if	(str_startsw(program.arg()[argi],"--filt-gba"))		ReadArg(program.arg(),argi,cut_gba);
		else if (var_in.empty()) var_in=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	program.help_text_var("_Default_det",det_in);
	program.help_text_var("_Default_log",log_fn);
	program.help_text_var("_Default_obs",itos(cut_obs));
	program.help_text_var("_Default_seg",ftos(cut_seg));
	program.help_text_var("_Default_aaa",ftos(cut_aaa));
	program.help_text_var("_Default_gba",ftos(cut_gba));
	perch::check_arguments();

	set<string> incl_gene_symbols;
	if (!sym_in.empty())
	{
		for (Rows_in_File(in,sym_in,1)) incl_gene_symbols.insert(in[0]);
	}
	
	map<string, double> background1; // background1[GeneSymbol]=max_BayesDel
	if (!det_in.empty())
	{
		field_numbers	FldSym(false,true);	// field numb for GeneSymbol
		field_numbers	FldDel(false,true);	// field numb for BayesDel
		field_numbers	FldAAA(false,true);	// field numb for AAA_details
		tfile_format det_format;
		det_format.set_delimiters("\t");
		det_format.set_option(SKIP_NOTES,false);
		for (Rows_in_File(in,det_in,&det_format))
		{
			if (exist_any(perch::h_col1, in.contents()))
			{
				FldSym.clear();
				FldDel.clear();
				FldAAA.clear();
				for (int i=0;i<in.NumFields();++i)
				{
					if (in[i]==perch::h_symb)	FldSym.push_back(i+1);
					if (in[i]=="BayesDel")		FldDel.push_back(i+1);
					if (in[i]=="AAA_details")	FldAAA.push_back(i+1);
				}
				if (FldSym.no_input())	exit_error("The "+perch::h_symb+" column is missing.");
				if (FldDel.no_input())	exit_error("The BayesDel column is missing.");
				if (FldAAA.no_input())	exit_error("The AAA_details column is missing.");
				continue;
			}
			if (FldSym.no_input()) exit_error("header row not read from "+det_in);
			if (str_has(in[FldAAA[0]],"SKIPPED")) continue;
			int num_obs=0;
			vector<string> AAA_details;
			boost::split(AAA_details,in[FldAAA[0]],boost::is_any_of(";"));
			for (auto &x:AAA_details)
			{
				if ((str_startsw(x,"het_ct:")||str_startsw(x,"hom_ct:"))&&!str_endsw(x,":"))
				{
					vector<string> samples;
					x=x.substr(7);
					boost::split(samples,x,boost::is_any_of(","));
					num_obs+=samples.size();
				}
			}
			if (num_obs<cut_obs) continue;

			double BayesDel = std::numeric_limits<double>::signaling_NaN();
			read_val(in[FldDel[0]],BayesDel);
			
			string GeneSymbol;
			if (!FldSym.no_input()) GeneSymbol=in[FldSym[0]];
			if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
			if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
			
			if (!std::isnan(BayesDel))
			{
				if (exist_element(background1,GeneSymbol)) { if (BayesDel>background1[GeneSymbol]) background1[GeneSymbol]=BayesDel; }
				else background1[GeneSymbol]=BayesDel;
			}
		}
	}
	
	{
		boost::iostreams::filtering_ostream logout;			// file, set by main if !log_fn.empty()
		if (!log_fn.empty()) { if (!openOutFile(logout, log_fn)) exit_error("cannot open log file."); }
		field_numbers	FldSym(false,true);	// field numb for GeneSymbol
		field_numbers	FldFun(false,true);	// field numb for FuncConseq
		field_numbers	FldGrp(false,true);	// field numb for #CHROM:POS:REF:ALT:Func_Gene
		field_numbers	FldDel(false,true);	// field numb for BayesDel
		field_numbers	FldAAA(false,true);	// field numb for AAA_details
		field_numbers	FldSEG(false,true);	// field numb for AAA_details
		field_numbers	FldGBA(false,true);	// field numb for AAA_details
		tfile_format var_format;
		var_format.set_delimiters("\t");
		var_format.set_option(SKIP_NOTES,false);
		for (Rows_in_File(in,var_in,&var_format))
		{
			if (exist_any(perch::h_col1, in.contents()))
			{
				FldSym.clear();
				FldFun.clear();
				FldGrp.clear();
				FldDel.clear();
				FldAAA.clear();
				FldSEG.clear();
				FldGBA.clear();
				for (int i=0;i<in.NumFields();++i)
				{
					if (in[i]==perch::h_symb)						FldSym.push_back(i+1);
					if (in[i]==perch::h_func)						FldFun.push_back(i+1);
					if (in[i]=="#CHROM:POS:REF:ALT:"+perch::h_symb)	FldGrp.push_back(i+1);
					if (in[i]=="BayesDel")							FldDel.push_back(i+1);
					if (in[i]=="AAAlbf")							FldAAA.push_back(i+1);
					if (in[i]=="SEGlbf")							FldSEG.push_back(i+1);
					if (in[i]=="GBAlbf")							FldGBA.push_back(i+1);
				}
				if (FldDel.no_input())	exit_error("The BayesDel column is missing.");
				if (FldGrp.no_input())	exit_error("The #CHROM:POS:REF:ALT:"+perch::h_symb+" column is missing.");
				in.write_r(program.outf,true);
				continue;
			}
			if (FldDel.no_input()) exit_error("header row not read from "+var_in);
			
			string GeneSymbol;
			if (!FldSym.no_input()) GeneSymbol=in[FldSym[0]];
			else if (!FldGrp.no_input()) GeneSymbol=substr_after_rfind(in[FldGrp[0]],":");
			if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
			if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }
			if (!incl_gene_symbols.empty()) { if (!exist_element(incl_gene_symbols,GeneSymbol)) continue; }
			
			double BayesDel = std::numeric_limits<double>::signaling_NaN();
			double SEGlbf = std::numeric_limits<double>::signaling_NaN();
			double AAAlbf = std::numeric_limits<double>::signaling_NaN();
			double GBAlbf = std::numeric_limits<double>::signaling_NaN();
			if (!FldDel.no_input()) read_val(in[FldDel[0]],BayesDel);
			if (!FldSEG.no_input()) read_val(in[FldSEG[0]],SEGlbf);
			if (!FldAAA.no_input()) read_val(in[FldAAA[0]],AAAlbf);
			if (!FldGBA.no_input()) read_val(in[FldGBA[0]],GBAlbf);
			bool is_lof = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"LoF")||str_has(in[FldFun[0]],"NMD") );
			bool is_vks = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"knClinSig=1") );
			bool is_DNg = ( FldFun.no_input() ? false : perch::is_DomNeg(in[FldFun[0]]) );
			bool is_cds = ( FldFun.no_input() ? false : perch::is_coding(in[FldFun[0]]) );
			
			string reasons;
			if (perch::LFonly && is_lof && perch::is_LoFtol(GeneSymbol))		reasons+="LoF-tolerant;";
			if (perch::LFonly && !is_lof)										reasons+="NotLoF;";
			if (perch::CDonly && !is_cds)										reasons+="NotCDS;";
			if (perch::DomNeg && !is_DNg)										reasons+="NotDomNeg;";
			if (!std::isinf(cut_seg) && !std::isnan(SEGlbf)	&& SEGlbf<cut_seg)	reasons+="SEGlbf;";
			if (!std::isinf(cut_aaa) && !std::isnan(AAAlbf)	&& AAAlbf<cut_aaa)	reasons+="AAAlbf;";
			if (!std::isinf(cut_gba) && !std::isnan(GBAlbf)	&& GBAlbf<cut_gba)	reasons+="GBAlbf;";
			if (perch::filter_AnnAF(BayesDel,GeneSymbol,is_lof,is_vks))			reasons+="BayesDel;";
			if (exist_element(background1,GeneSymbol) && !std::isnan(BayesDel) && background1[GeneSymbol]>=BayesDel && !is_lof && !is_vks) reasons+="observation;";
			if (!reasons.empty()) { reasons.pop_back(); if (!log_fn.empty()) logout<<in[FldGrp[0]]<<'\t'<<reasons<<endl; continue; }
			
			in.write_r(program.outf,true);
		}
	}
	return 0;
}
