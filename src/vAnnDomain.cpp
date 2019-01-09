#include <tft/libfbj_file.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// parameters
	string			d_File = "@GDB_InterPro";	// annotation file:   "@GDB_CDD"  "@GDB_PFam"  "@GDB_InterPro"
	string			wrWhat;						// annotation header: "CDD"       "PFam"       "InterPro"
	string			inputs;						// input genotype file in VCF format
	string			MisStr = "NA";				// missing value
	char			DelChr = '|';				// delimiter between multiple domains
	bool			AddInf = true;				// add a field in INFO, otherwise add a column
	set<string>		xFtype = {"Intergenic","Intronic","Downstream","Upstream","UTR3","UTR5","miRNAbinding","ncRNA_SNV","ncRNA_InDel"};
	
	// handle program arguments
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"--add-info"))	{ ReadArg(program.arg(),argi,AddInf); }
		else if (str_startsw(program.arg()[argi],"--ann"))		{ ReadArg(program.arg(),argi,d_File); }
		else if (str_startsw(program.arg()[argi],"--wr"))		{ ReadArg(program.arg(),argi,wrWhat); }
		else if (inputs.empty()) inputs=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	program.help_text_var("_Default_wr","automatically decided from filename");
	program.help_text_var("_Default_domain",d_File);
	program.help_text_var("_Default_add_info",str_YesOrNo(AddInf));
	perch::check_arguments();
	
	// prepare
	if (d_File.empty())
	{
		tfile_format format;
		format.set_delimiters("\t");
		format.set_option(SKIP_NOTES,false);
		for (Lines_in_File(in,inputs,&format)) program.outf << in[0] << endl;
		return 0;
	}
	else
	{
		if (wrWhat.empty())
		{
			string short_name = substr_after_rfind(d_File,"/");
			boost::algorithm::replace_all(short_name,"@GDB_","");
			wrWhat=short_name;
		}
		d_File=perch::find_file(d_File);
	}
	
	// read gene db, domain db
	genepi::read_genes();
	typedef std::tuple<int,int,string> domain_t;
	map< string, set<domain_t> > pfam_db;		// pfam_db[Tx]=set<lb,ub,ID>. Use set to remove duplicated records.
	tfile_format pfam_format;
	pfam_format.set_option(SKIP_NOTES,false);
	pfam_format.set_field_nums(4,"lines do not have sufficient columns",tfile_format::Exit);
	bool is_aa=false;
	bool is_cds=false;
	bool is_rna=false;
	bool is_dna=false;
	for (Rows_in_File(in,d_File,&pfam_format))
	{
		if (!in[0].empty()&&in[0][0]=='#')
		{
			string inType = boost::to_upper_copy(in[1]);
			if		(str_startsw(inType,"RNA")) is_rna=true;
			else if (str_startsw(inType,"CDS")) is_cds=true;
			else if (str_startsw(inType,"DNA")) is_dna=true;
			else if (str_startsw(inType,"AA"))  is_aa=true;
			else	exit_error("Input type wrong. Column 2,3 header should start with RNA/CDS/AA/DNA.");
			continue;
		}
		if (!is_aa && !is_cds && !is_rna && !is_dna) exit_error("The header row is missing in "+d_File);
		auto it = genepi::gene_byTranscript.find(in[0]);
		if (it==genepi::gene_byTranscript.end()) continue;
		genepi::gene_info& g=it->second.begin()->second;
		int loc1; if (!read_val(in[1],loc1)) exit_error("Failed to read "+in[1]+" as a position.");
		int loc2; if (!read_val(in[2],loc2)) exit_error("Failed to read "+in[2]+" as a position.");
		string chr;
		int lb,ub;
		if		(is_aa)
		{
			int locA=(loc1-1)*3+1;
			int locB=loc2*3;
			if (!cds_to_genomic_location(g,locA,chr,lb)) { if (perch::_Debug) lns<<showe<<"wrong AA position "<<loc1<<" for "<<in[0]<<flush_logger; else exit_error("wrong position "+itos(loc1)+" for "+in[0]); }
			if (!cds_to_genomic_location(g,locB,chr,ub)) { if (perch::_Debug) lns<<showe<<"wrong AA position "<<loc2<<" for "<<in[0]<<flush_logger; else exit_error("wrong position "+itos(loc2)+" for "+in[0]); }
		}
		else if (is_cds)
		{
			if (!cds_to_genomic_location(g,loc1,chr,lb)) { if (perch::_Debug) lns<<showe<<"wrong CDS position "<<loc1<<" for "<<in[0]<<flush_logger; else exit_error("wrong position "+itos(loc1)+" for "+in[0]); }
			if (!cds_to_genomic_location(g,loc2,chr,ub)) { if (perch::_Debug) lns<<showe<<"wrong CDS position "<<loc2<<" for "<<in[0]<<flush_logger; else exit_error("wrong position "+itos(loc2)+" for "+in[0]); }
		}
		else if (is_rna)
		{
			if (!rna_to_genomic_location(g,loc1,chr,lb)) { if (perch::_Debug) lns<<showe<<"wrong RNA position "<<loc1<<" for "<<in[0]<<flush_logger; else exit_error("wrong position "+itos(loc1)+" for "+in[0]); }
			if (!rna_to_genomic_location(g,loc2,chr,ub)) { if (perch::_Debug) lns<<showe<<"wrong RNA position "<<loc2<<" for "<<in[0]<<flush_logger; else exit_error("wrong position "+itos(loc2)+" for "+in[0]); }
		}
		else
		{
			lb=loc1;
			ub=loc2;
		}
		if (lb>ub) swap(lb,ub);
		pfam_db[in[0]].insert(make_tuple(lb,ub,in[3]));
	}
	
	// read VCF
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldFun(false,true);	// field numb for FuncConseq
	field_numbers	FldSym(false,true);	// field numb for GeneSymbol
	field_numbers	FldRes(false,true);	// field numb for result
	tfile_format	format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);
	bool header_not_read = true;
	bool VCFHeader_not_w = true;
	int AnnGen=0; // 0: no vAnnGene in INFO; 1: vAnnGene in INFO
	for (Rows_in_File(in, inputs, &format))
	{
		// read header
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			in.clear_nf();
			if (str_has(in[0],"##INFO=<ID="+perch::i_func+",")) AnnGen=1;
			if (str_has(in[0],"##INFO=<ID="+wrWhat+",")) VCFHeader_not_w=false;
			if (!header_not_read) continue;
			print_container(in.contents(),program.outf,' ',true);
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			if (AddInf && VCFHeader_not_w)
			{
				program.outf<<"##INFO=<ID="+wrWhat+",Number=1,Type=String,Description=\"Protein domain in "<<d_File<<"\">"<<endl;
				VCFHeader_not_w=false;
			}
			if (!header_not_read) continue;
			header_not_read=false;
			format.clear_field_nums();
			FldChr.clear();
			FldPos.clear();
			FldRef.clear();
			FldAlt.clear();
			FldInf.clear();
			FldFun.clear();
			FldSym.clear();
			FldRes.clear();
			for (int i=0;i<in.NumFields();++i)
			{
				if (in[i]==perch::h_func)	FldFun.push_back(i+1);
				if (in[i]==perch::h_symb)	FldSym.push_back(i+1);
				if (in[i]=="INFO")	FldInf.push_back(i+1);
				if (in[i]=="Chr"   && FldChr.no_input())	FldChr.push_back(i+1); if (in[i]=="#CHROM")		{	FldChr.clear(); FldChr.push_back(i+1); }
				if (in[i]=="Start" && FldPos.no_input())	FldPos.push_back(i+1); if (in[i]=="POS")		{	FldPos.clear(); FldPos.push_back(i+1); }
				if (in[i]=="Ref"   && FldRef.no_input())	FldRef.push_back(i+1); if (in[i]=="REF")		{	FldRef.clear(); FldRef.push_back(i+1); }
				if (in[i]=="Alt"   && FldAlt.no_input())	FldAlt.push_back(i+1); if (in[i]=="ALT")		{	FldAlt.clear(); FldAlt.push_back(i+1); }
			}
			if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
			if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
			if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
			if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");			
			if ((FldSym.no_input()||FldFun.no_input()) && !AnnGen) exit_error("Functional consequence annotation is missing.");
			if (AddInf)	{ if (FldInf.no_input()) exit_error("The INFO column is missing."); FldRes=FldInf; }
			else if (FldRes.no_input()) { FldRes.push_back(in.NumFields()+1); in.contents().push_back(wrWhat); format.set_field_nums(FldRes,"",tfile_format::Expand); }
			
			print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
		if (header_not_read) exit_error("Header lines missing.");

		// check file
		string this_index = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+in[FldRef[0]]+DLMTR+in[FldAlt[0]];
		
		// read INFO
		vector<string> INFO;
		bool info_modified=false;
		if (!FldInf.no_input())
		{
			if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
			for (vector<string>::iterator it = INFO.begin(); it != INFO.end(); it++)
				if (str_startsw(*it,wrWhat+"=")) { INFO.erase(it); info_modified=true; break; }
		}
		if (!AddInf) in[FldRes[0]] = MisStr;
		string SVTYPE=get_string(INFO,"SVTYPE");
		string END = get_string(INFO,"END");
		
		// basic info
		string& ref=in[FldRef[0]];
		string& alt=in[FldAlt[0]];
		int	chr_num;	if (!genepi::read_chr_num(in[FldChr[0]],chr_num))	exit_error("Failed to read "+in[FldChr[0]]+" as a chromosome.");
		int	bp;			if (!read_val_ge(in[FldPos[0]],bp,1))				exit_error("Failed to read "+in[FldPos[0]]+" as a position in basepairs.");
		bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
		bool is_sv = (!SVTYPE.empty() && str_has(in[FldAlt[0]],"<"));
		int b2;
		if (is_sv && (SVTYPE=="INV"||SVTYPE=="DEL"||SVTYPE=="DUP"||SVTYPE=="CNV")) { if (!ReadStr(END,b2,0)) exit_error("No INFO:END for "+this_index); }
		else if (is_snv) b2=bp;
		else b2 = bp+in[FldRef[0]].size()-1;
		
		// annotate protein domain
		string GeneSymbol = "NotAnnotated";
		string FuncConseq;
		if (!FldSym.no_input() && !FldFun.no_input())
		{
			GeneSymbol=in[FldSym[0]]; if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
			FuncConseq=to_lower_copy(in[FldFun[0]]);
		}
		else if (AnnGen)
		{
			string func_ann = get_string(INFO,perch::i_func);
			if (func_ann.empty()) exit_error(perch::i_func+" not in INFO or empty");
			vector<string> func_vec;
			boost::split(func_vec,func_ann,boost::is_any_of(","));
			if (func_vec.size()<3) exit_error(perch::i_func+" vector size wrong");
			GeneSymbol=func_vec[0]; if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
			FuncConseq=to_lower_copy(func_vec[1]);
		}
		
		if (!exist_element(xFtype,FuncConseq))
		{
			set<string> result_set;
			if (exist_element(genepi::gene_bySymbol,GeneSymbol))
			{
				for (auto &i:genepi::gene_bySymbol[GeneSymbol])
				{
					genepi::gene_info& g=i.second;
					if (exist_element(pfam_db,g.name) && overlap_exon(g,chr_num,bp,b2))
					{
						for (auto &d:pfam_db[g.name])
						{
							int lb=std::get<0>(d);
							int ub=std::get<1>(d);
							if (lb<=b2 && ub>=bp) result_set.insert(std::get<2>(d));
						}
					}
				}
			}
			// else result_set = {"NA"}; // gene symbol not found
			
			// print result
			string result = str_of_container(result_set,DelChr);
			if (AddInf)
			{
				INFO.push_back(wrWhat+"="+result);
				info_modified=true;
			}
			else
			{
				in[FldRes[0]] = result;
			}
		}
		
		if (info_modified) in[FldInf[0]]=str_of_container(INFO,';');
		print_container(in.contents(),program.outf,DLMTR,true);
	}
	return 0;
}
