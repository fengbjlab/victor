/*
 Removed option descriptions:
 --MaxAF-pr D      Integrate MaxAF into the annotated score with a proportion of D {_Default_MaxAF_pr}
 --af Ss           Combine MaxAF with allele frequencies in these INFO fields to calculate overall MaxAF {_Default_af}
 --add-af B        Add MaxAF to the annotated score (InDel only) {_Default_add_af}
*/

#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_program.hpp>
#include "vDEL.hpp"
#include "victor_par.hpp"

using namespace std;
typedef std::tuple<int,int,string,string> ID_t;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	string			MisVal = ".";
	int				mnCode = 0;
	int				mxCode = 249;
	int				msCode = 255;
	int				extend = 0;
	bool			AddInf = false;
	bool			add_af = false;
	vector<string>	inputs;
	string			anFile;
	string			wrWhat;
	int				per_bp = 0;
	double			minval = std::numeric_limits<double>::signaling_NaN();
	double			stpval = std::numeric_limits<double>::signaling_NaN();
	set<string>	header_chr = {"#CHROM","CHROM","#CHR","CHR"};
	set<string>	header_pos = {"POS","START","POSITION"};
	set<string>	header_ref = {"REF"};
	set<string>	header_alt = {"ALT"};
	bool	replace_old_ann=false;
	string			indels;
	tfile_format	format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);
	
	// handle program options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--ann"))		ReadArg(program.arg(),argi,anFile);
		else if (str_startsw(program.arg()[argi],"-x"))			ReadArg(program.arg(),argi,per_bp);
		else if (str_startsw(program.arg()[argi],"--min"))		ReadArg(program.arg(),argi,minval);
		else if (str_startsw(program.arg()[argi],"--step"))		ReadArg(program.arg(),argi,stpval);
		else if (str_startsw(program.arg()[argi],"--wr"))		ReadArg(program.arg(),argi,wrWhat);
		else if (str_startsw(program.arg()[argi],"--add-info"))	ReadArg(program.arg(),argi,AddInf);
		else if	(str_startsw(program.arg()[argi],"--replace"))	ReadArg(program.arg(),argi,replace_old_ann);
		else if	(str_startsw(program.arg()[argi],"--indel"))	ReadArg(program.arg(),argi,indels);
		else if	(str_startsw(program.arg()[argi],"--missing"))	ReadArg(program.arg(),argi,MisVal);
		else if	(str_startsw(program.arg()[argi],"--ms-code"))	ReadArg(program.arg(),argi,msCode);
		else if	(str_startsw(program.arg()[argi],"--padding"))	ReadArg(program.arg(),argi,extend);
		else if	(str_startsw(program.arg()[argi],"--MaxAF-pr"))	ReadArg(program.arg(),argi,AF_prp);
		else if	(str_startsw(program.arg()[argi],"--add-af"))	ReadArg(program.arg(),argi,add_af);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else add_to_container(inputs,program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_file",anFile);
	program.help_text_var("_Default_wr",wrWhat);
	program.help_text_var("_Default_min",ftos(minval));
	program.help_text_var("_Default_step",ftos(stpval));
	program.help_text_var("_Default_x",itos(per_bp));
	program.help_text_var("_Default_replace",str_YesOrNo(replace_old_ann));
	program.help_text_var("_Default_add_info",str_YesOrNo(AddInf));
	program.help_text_var("_Default_indels",indels);
	program.help_text_var("_Default_missing",MisVal);
	program.help_text_var("_Default_ms_code",itos(msCode));
	program.help_text_var("_Default_MaxAF_pr",ftos(AF_prp));
	program.help_text_var("_Default_add_af",str_YesOrNo(add_af));
	perch::check_arguments();
	
	// check errors
	if (std::isnan(minval)) exit_error("--min is required");
	if (std::isnan(stpval)) exit_error("--step is required");
	if (per_bp==0) exit_error("-x is required");
	if (str_endsw(anFile,"/")) anFile.pop_back();
	if (!anFile.empty() && !DirExists(anFile))
	{
		if (DirExists(perch::DBpath()+anFile)) anFile=perch::DBpath()+anFile;
		else exit_error("cannot find "+anFile);
	}
	if (!indels.empty() && indels!="max" && indels!="min") exit_error("If set, --indel should be max or min.");
	if (indels.empty()) lns<<showw<<"--indel was not set. Therefore, InDels will not be annotated."<<flush_logger;
	boost::to_lower(indels);
	
	// do nothing
	if (anFile.empty())
	{
		for (Lines_in_File(in,inputs,&format)) program.outf << in[0] << endl;
		return 0;
	}

	// prepare
	bool ann_indel_max = (indels=="max");
	bool ann_indel_min = (indels=="min");
	if (wrWhat.empty()) wrWhat=substr_after_rfind(anFile,"/");
	Interpolate	MaxAF_curve;
	MaxAF_curve.reference_data=ebf_nsfp33a_MaxAF;
	
	fstream inpf[50];
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldRes(false,true);	// field numb for result
	field_numbers	FldXAF(false,true);	// field numb for MaxAF
	bool header_not_read = true;
	int num_MaxAF_INFO = 0;
	set<string> exclude_chr;
	for (Rows_in_File(in, inputs, &format))
	{
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#') // VCFHeader
		{
			in.clear_nf();
			if (str_startsw(in[0],"##INFO=<ID="+wrWhat+"_includes_MaxAF,") && add_af) exit_error("--add-af is applied twice");
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_MxAF+","))	++num_MaxAF_INFO;
			else { for (auto &id : perch::h_afID) if (str_startsw(in[0],"##INFO=<ID="+id+","))	++num_MaxAF_INFO; }
			if (str_has(in[0],"##INFO=<ID="+wrWhat+",")) continue;
			print_container(in.contents(),program.outf,' ',true);
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			if (AddInf)
				program.outf<<"##INFO=<ID="+wrWhat+",Number=1,Type=Float,Description=\"Score annotated by vAnnBase\">"<<endl;
			if (add_af || !str_has(anFile,"_noAF"))
				program.outf<<"##INFO=<ID="+wrWhat+"_includes_MaxAF,Number=0,Type=Flag,Description=\""+wrWhat+" computed with --add-af.\">"<<endl;
			header_not_read=false;
			format.clear_field_nums();
			FldInf.clear();
			FldChr.clear();
			FldPos.clear();
			FldRef.clear();
			FldAlt.clear();
			FldRes.clear();
			FldXAF.clear();
			for (int i=0;i<in.NumFields();++i)
			{
				if (in[i]==wrWhat) 			{ if (FldRes.no_input()) FldRes.push_back(i+1); else exit_error("Multilpe "+wrWhat+" columns."); }
				if (in[i]==perch::h_MxAF)	FldXAF.push_back(i+1); else { for (auto &id : perch::h_afID) if (in[i]==id) FldXAF.push_back(i+1); }
				if (in[i]=="INFO")			FldInf.push_back(i+1);
				string in_i_ = boost::to_upper_copy(in[i]);
				if (exist_element(header_chr,in_i_))	{ if (FldChr.no_input()) FldChr.push_back(i+1); else exit_error("Multilpe CHR columns."); }
				if (exist_element(header_pos,in_i_))	{ if (FldPos.no_input()) FldPos.push_back(i+1); else exit_error("Multilpe POS columns."); }
				if (exist_element(header_ref,in_i_))	{ if (FldRef.no_input()) FldRef.push_back(i+1); else exit_error("Multilpe REF columns."); }
				if (exist_element(header_alt,in_i_))	{ if (FldAlt.no_input()) FldAlt.push_back(i+1); else exit_error("Multilpe ALT columns."); }
			}
			if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
			if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
			if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
			if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
			if (AddInf) { if (FldInf.no_input()) exit_error("The INFO column is missing."); }
			else if (FldRes.no_input())	{ FldRes.push_back(in.NumFields()+1); in.contents().push_back(wrWhat); format.set_field_nums(FldRes,"",tfile_format::Expand); }
			if (add_af && num_MaxAF_INFO==0 && FldXAF.no_input()) exit_error("MaxAF annotation is required by --add-af but is missing");
			format.set_field_nums(FldChr,"lines missing the #CHROM column.",tfile_format::Continue);
			format.set_field_nums(FldPos,"lines missing the POS column.",tfile_format::Continue);
			print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
		if (header_not_read) exit_error("Header lines missing.");
		if (exist_element(exclude_chr,in[FldChr[0]]))
		{
			if (!AddInf) in[FldRes[0]] = MisVal;
			print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
			
		string& ref = in[FldRef[0]];
		string& alt = in[FldAlt[0]];
		if (ref.empty()||alt.empty()||ref=="."||alt=="."||ref=="-"||alt=="-")
			exit_error("The REF/ALT doesn't seem to follow VCF format. Culprit: "+in[FldChr[0]]+","+in[FldPos[0]]+","+ref+","+alt);
		bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
		if (is_snv || ann_indel_max || ann_indel_min)
		{
			// check whether already annotated
			vector<string> INFO;
			if (!FldInf.no_input())
				if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));

			double MaxAF_lod=0;
			if (AF_prp && add_af)
			{
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
				if (!std::isnan(MaxAF)) MaxAF_lod = log10(MaxAF_curve.solve(MaxAF));
			}

			bool already_annotated = false;
			if (AddInf)
			{
				for (auto& f:INFO) if (str_startsw(f,wrWhat+"=")) { already_annotated=true; break; }
			}
			else
			{
				if (!in[FldRes[0]].empty() && in[FldRes[0]]!=MisVal) already_annotated=true;
			}
			if (!replace_old_ann && already_annotated)
			{
				print_container(in.contents(),program.outf,DLMTR,true);
				continue;
			}

			// open annotation file
			int	chr_num;	if (!genepi::read_chr_num(in[FldChr[0]],chr_num))	exit_error("Failed to read "+in[FldChr[0]]+" as a chromosome.");
			int	bp;			if (!read_val_ge(in[FldPos[0]],bp,1))				exit_error("Failed to read "+in[FldPos[0]]+" as a position in basepairs.");
			if (chr_num>=50) exit_error("Chr number out of bound");
			string chrName = genepi::convert_chr_num(chr_num);
			if (!inpf[chr_num].is_open())
			{
				if ( !openfile_successfully(inpf[chr_num],anFile+"/"+chrName,ios::in|ios::binary))
				{
					exclude_chr.insert(in[FldChr[0]]);
					if (!AddInf) in[FldRes[0]] = MisVal;
					print_container(in.contents(),program.outf,DLMTR,true);
					continue;
				}
			}
			
			// do annotation
			string result;
			if (is_snv)
			{
				unsigned char read_byte;
				int offset=0;
				if (per_bp==3)
				{
					if (ref[0]=='A') { if (alt[0]=='G') offset=1; else if (alt[0]=='T') offset=2; }
					if (ref[0]=='C') { if (alt[0]=='G') offset=1; else if (alt[0]=='T') offset=2; }
					if (ref[0]=='G') { if (alt[0]=='C') offset=1; else if (alt[0]=='T') offset=2; }
					if (ref[0]=='T') { if (alt[0]=='C') offset=1; else if (alt[0]=='G') offset=2; }
				}
				inpf[chr_num].seekg((bp-1)*per_bp+offset, std::ios::beg);
				inpf[chr_num].read((char*)&read_byte,sizeof(unsigned char));
				if ((int)read_byte==msCode)	result = MisVal;
				else
				{
					double score = read_byte*stpval+minval;
					if (MaxAF_lod) score = AF_prp*MaxAF_lod + (1-AF_prp)*score;
					result = ftos(score);
				}
			}
			else
			{
				unsigned char read_byte;
				int pos1 = bp;
				int pos2 = bp;
				if		(ref.size()==1 && alt.size()>1) ++pos2;
				else if (ref.size()>1) pos2=pos1+ref.size()-1;
				else exit_error("impossible");
				pos1-=extend;
				pos2+=extend;
				unsigned char max_byte=mnCode;
				unsigned char min_byte=mxCode;
				int read_something=0;
				for (int i=pos1; i<=pos2; ++i)
				{
					for (int j=0;j<per_bp;++j)
					{
						inpf[chr_num].seekg((i-1)*per_bp+j, std::ios::beg);
						inpf[chr_num].read((char*)&read_byte,sizeof(unsigned char));
						if ((int)read_byte==msCode) continue;
						if (read_byte>max_byte) max_byte=read_byte;
						if (read_byte<min_byte) min_byte=read_byte;
						++read_something;
					}
				}
				if (read_something)
				{
					double score;
					if		(ann_indel_max) score = (max_byte*stpval+minval);
					else if (ann_indel_min) score = (min_byte*stpval+minval);
					else exit_error("impossible");
					if (MaxAF_lod) score = AF_prp*MaxAF_lod + (1-AF_prp)*score;
					result = ftos(score);
				}
				else
				{
					result = MisVal;
				}
			}
			
			// write results
			if (AddInf)
			{
				bool update=false;
				for (auto& f:INFO)
				{
					if (str_startsw(f,wrWhat+"="))
					{
						f = wrWhat + "=" + result;
						update=true;
						break;
					}
				}
				if (!update)
				{
					if (in[FldInf[0]].empty() || in[FldInf[0]]==".")	in[FldInf[0]]  =       wrWhat + "=" + result;
					else												in[FldInf[0]] += ";" + wrWhat + "=" + result;
				}
				else
					in[FldInf[0]]=str_of_container(INFO,';');
			}
			else
				in[FldRes[0]] = result;
		}
		else
		{
			if (!AddInf) in[FldRes[0]] = MisVal;
		}
		print_container(in.contents(),program.outf,DLMTR,true);
	}
	return 0;
}
