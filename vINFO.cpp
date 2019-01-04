#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	vector<string>	inputs;
	set<string>	header_chr = {"#CHROM","CHROM","#CHR","CHR"};
	set<string>	header_pos = {"POS","START","POSITION"};
	set<string>	header_ref = {"REF"};
	set<string>	header_alt = {"ALT"};
	set<string>		toKeep;
	set<string>		toRemv;
	bool			toRest=false;
	tfile_format	format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);
	
	// handle program options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--keep"))		ReadSet(program.arg(),argi,toKeep);
		else if	(str_startsw(program.arg()[argi],"--remove"))	ReadSet(program.arg(),argi,toRemv);
		else if	(str_startsw(program.arg()[argi],"--restore"))	ReadArg(program.arg(),argi,toRest);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else add_to_container(inputs,program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_keep",str_of_container(toKeep,',',false));
	program.help_text_var("_Default_remove",str_of_container(toRemv,',',false));
	program.help_text_var("_Default_restore",str_YesOrNo(toRest));
	perch::check_arguments();
	
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	bool header_not_read = true;
	for (Rows_in_File(in, inputs, &format))
	{
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#') // VCFHeader
		{
			in.clear_nf();
			if (!toRemv.empty())
			{
				bool skip_this=false;
				for (auto &id:toRemv)
					if (str_startsw(in[0],"##INFO=<ID="+id+",")) { skip_this=true; break; }
				if (skip_this) continue;
			}
			if (!toKeep.empty() && str_startsw(in[0],"##INFO=<ID="))
			{
				bool skip_this=true;
				for (auto &id:toKeep)
					if (str_startsw(in[0],"##INFO=<ID="+id+",")) { skip_this=false; break; }
				if (skip_this) continue;
			}
			print_container(in.contents(),program.outf,' ',true);
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			header_not_read=false;
			format.clear_field_nums();
			FldInf.clear();
			FldChr.clear();
			FldPos.clear();
			FldRef.clear();
			FldAlt.clear();
			for (int i=0;i<in.NumFields();++i)
			{
				string in_i_ = boost::to_upper_copy(in[i]);
				if (in_i_=="INFO")	FldInf.push_back(i+1);
				if (exist_element(header_chr,in_i_))	{ if (FldChr.no_input()) FldChr.push_back(i+1); else exit_error("Multilpe CHR columns."); }
				if (exist_element(header_pos,in_i_))	{ if (FldPos.no_input()) FldPos.push_back(i+1); else exit_error("Multilpe POS columns."); }
				if (exist_element(header_ref,in_i_))	{ if (FldRef.no_input()) FldRef.push_back(i+1); else exit_error("Multilpe REF columns."); }
				if (exist_element(header_alt,in_i_))	{ if (FldAlt.no_input()) FldAlt.push_back(i+1); else exit_error("Multilpe ALT columns."); }
			}
			if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
			if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
			if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
			if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
			if (FldInf.no_input()) exit_error("The INFO column is missing.");
			format.set_field_nums(FldChr,"lines missing the #CHROM column.",tfile_format::Continue);
			format.set_field_nums(FldPos,"lines missing the POS column.",tfile_format::Continue);
			format.set_field_nums(FldRef,"lines missing the REF column.",tfile_format::Continue);
			format.set_field_nums(FldAlt,"lines missing the ALT column.",tfile_format::Continue);
			format.set_field_nums(FldInf,"lines missing the INFO column.",tfile_format::Continue);
			print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
		if (header_not_read) exit_error("Header lines missing.");
		
		vector<string> INFO;
		if (!FldInf.no_input())
			if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
		if (toRest)
		{
			for (auto& f:INFO)
				if (str_startsw(f,"OriginalIndex="))
				{
					string original_index=f.substr(14);
					vector<string> components;
					boost::split(components,original_index,boost::is_any_of("_"));
					if (components.size()==4)
					{
						in[FldChr[0]]=components[0];
						in[FldPos[0]]=components[1];
						in[FldRef[0]]=components[2];
						in[FldAlt[0]]=components[3];
					}
					break;
				}
		}
		if (!toKeep.empty() || !toRemv.empty())
		{
			vector<string> temp;
			for (auto& f:INFO)
			{
				string id = substr_before_find(f,"=");
				if (exist_element(toRemv,id)) continue;
				if (!toKeep.empty() && !exist_element(toKeep,id)) continue;
				temp.push_back(f);
			}
			INFO=temp;
			in[FldInf[0]] = str_of_container(INFO,';');
			if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
		}
		
		print_container(in.contents(),program.outf,DLMTR,true);
	}
	return 0;
}
