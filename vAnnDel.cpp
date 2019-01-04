/*
Caveat: 
 1. not robust to FILE2 not sorted.
 
 --check-order B   Check whether input is sorted by Chr,Pos,Alt {_Default_check_order}
 --chr STRs        Header string of the Chr column {_Default_chr}
 --pos STRs        Header string of the Pos column {_Default_pos}
 --ref STRs        Header string of the Ref column {_Default_ref}
 --alt STRs        Header string of the Alt column {_Default_alt}
 */
#include <tuple>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_genepi.hpp>
#include "victor_par.hpp"

using namespace std;

// -------------- modify the following for a specific project --------------
// The most important things to modify are all that relate to the keys for comparison.
// Also modify program.trademark in main(), which will be shown in program's help text.
// Both FILE1 and FILE2 need to be sorted by the keys; no checking in running.

vector<int> kf1 = {1,2,3,4};	// key fields file1 (Chr,Pos,Ref,Alt), although initialized they will be replaced anyway
vector<int> kf2 = {1,2,3,4};	// key fields file2 (Chr,Pos,Ref,Alt), although initialized they will be replaced anyway
struct data {
	int				key1;				// Chr, order is 1-22,X,Y,M
	int				key2;				// Pos, an integer as basepair
	string			key3;				// Ref
	string			key4;				// Alt
	field_numbers	OutFields;			// fields for output
	tfile_format	InpFormat;			// input format
	
	// initialize with the smallest impossible value, so that "*this < xxx" always return true
	data():key1(0),key2(0),key3(""),key4(""),OutFields(false,true) { }
	
	// set the keys to the largest  impossible value, so that "xxx < *this" always return true
	void set_max() { key1=INT_MAX; key2=INT_MAX; } // no need for all keys
};

// http://stackoverflow.com/questions/3882467/defining-operator-for-a-struct
// http://stackoverflow.com/questions/11312448/operator-comparing-multiple-fields
// http://stackoverflow.com/questions/6218812/implementing-comparision-operators-via-tuple-and-tie-a-good-idea
inline bool operator <(const data& x, const data& y) {
	return std::tie(x.key1, x.key2, x.key3, x.key4) < std::tie(y.key1, y.key2, y.key3, y.key4);
}

inline int read_keys_from(tabular_file& tfile, vector<int>& kf, data& to) {
	// check validity, don't mess up with "to" yet
	int key1; if (!genepi::read_chr_num(tfile[kf[0]],key1))	return 0;
	int key2; if (!read_val_ge(tfile[kf[1]],key2,1))		return 0;
	// save to "to"
	to.key1=key1;
	to.key2=key2;
	to.key3=tfile[kf[2]];
	to.key4=tfile[kf[3]];
	return 1;
}

// -------------- modify the above codes for a specific project --------------

data dataA; // used to check whether FILE1 is sorted
data data1;
data data2;
vector<string> header2;
bool AddInf = false;
string	MisStr=".";	// index or anything else
int ColInf=-1;
bool remove_old_ann=false;
bool add_ms=false;

inline bool already_annotated(tabular_file& in1) {
	bool annotated=true;
	if (data1.OutFields.empty())
	{
		if (AddInf)
		{
			for (each_element(data2.OutFields,it))
				if (!str_has(in1[ColInf],";"+header2[*it]+"=")) { annotated=false; break; } // not exctly right but OK most of the time
		}
		else
		{
			annotated=false;
		}
	}
	else
	{
		for (size_t i=0;i<data1.OutFields.size();++i)
			if (in1[data1.OutFields[i]].empty() || in1[data1.OutFields[i]]==MisStr) { annotated=false; break; }
	}
	return annotated;
}

// always annotate, replace the old
inline void write_matched(tabular_file& in1, vector<string>& in2, ostream& out) {
	if (data1.OutFields.empty())
	{
		if (AddInf)
		{
			vector<string> INFO;
			if (!in1[ColInf].empty()&&in1[ColInf]!=MisStr) boost::split(INFO,in1[ColInf],boost::is_any_of(";"));
			for (each_element(data2.OutFields,it))
			{
				bool exist=false;
				for (auto &i:INFO)
					if (str_startsw(i,header2[*it]+"=")) { i=header2[*it]+"="+in2[*it]; exist=true; break; }
				if (!exist) INFO.push_back(header2[*it]+"="+in2[*it]);
				in1[ColInf]=str_of_container(INFO,';');
			}
			in1.write_r(out,true);
		}
		else
		{
			in1.write_r(out);
			for (each_element(data2.OutFields,it)) out<<data1.InpFormat.output_del()<<in2[*it];
			out<<endl;
		}
	}
	else
	{
		for (size_t i=0;i<data1.OutFields.size();++i)
			vec_deq_set(in1.contents(),data1.OutFields[i],in2[data2.OutFields[i]]);
		in1.write_r(out,true);
	}
}

// when remove_old_ann=no: if annotation already exist, no change; if not, write str.
inline void write_uniqFILE1(tabular_file& in1, string str, ostream& out, bool force) {
	if (str=="index") str = in1[kf1[0]]+"_"+in1[kf1[1]]+"_"+in1[kf1[2]]+"_"+in1[kf1[3]];
	if (data1.OutFields.empty())
	{
		if (AddInf)
		{
			vector<string> INFO;
			if (!in1[ColInf].empty()&&in1[ColInf]!=MisStr) boost::split(INFO,in1[ColInf],boost::is_any_of(";"));
			for (each_element(data2.OutFields,it))
			{
				for (vector<string>::iterator info = INFO.begin(); info != INFO.end(); info++)
					if (str_startsw(*info,header2[*it]+"=")) { if (remove_old_ann || force) INFO.erase(info); break; }
				if (force) INFO.push_back(header2[*it]+"="+str); // to be consistent with already_annotated(), don't write if str is MisStr.
				in1[ColInf]=str_of_container(INFO,';');
			}
			in1.write_r(out,true);
		}
		else
		{
			in1.write_r(out);
			for (each_element(data2.OutFields,it)) out<<data1.InpFormat.output_del()<<str;
			out<<endl;
		}
	}
	else
	{
		for (size_t i=0;i<data1.OutFields.size();++i)
			if (in1.NumFields()<=data1.OutFields[i] || remove_old_ann || force) vec_deq_set(in1.contents(),data1.OutFields[i],str);
		in1.write_r(out,true);
	}
}

void annotate(tabular_file& in1, string FILE2, bool print_unmatched, const string& MisStr, int& foundMatch, int& fndNoMatch)
{
	string tbxres; // tabix results
	string& chr = in1[kf1[0]];
	string& pos = in1[kf1[1]];
	try { tbxres = exec("tabix "+FILE2+" "+chr+":"+pos+"-"+pos,false); }
	catch (const std::exception& error) { exit_error("tabix "+FILE2+" "+chr+":"+pos+"-"+pos+" failed"); }
	int matched=0;
	if (!tbxres.empty())
	{
		tbxres.pop_back(); // tbxres ends with \n
		vector<string> tbxrow;
		boost::split(tbxrow,tbxres,boost::is_any_of("\n"));
		for (auto &r:tbxrow)
		{
			vector<string> tbxstr;
			boost::split(tbxstr,r,boost::is_any_of("\t "));
			if (tbxstr[kf2[1]]==in1[kf1[1]] && tbxstr[kf2[2]]==in1[kf1[2]] && tbxstr[kf2[3]]==in1[kf1[3]]) // matched. Be careful, tabix may return lines with wrong bp!
			{
				// if (matched) exit_error("duplicated rows in "+FILE2); // No checking is faster. If want to check, remove the "break" below.
				write_matched(in1,tbxstr,program.outf);
				elog.add(foundMatch);
				++matched;
				break;
			}
		}
	}
	if (!matched)
	{
		if (print_unmatched) write_uniqFILE1(in1,MisStr,program.outf,add_ms);
		elog.add(fndNoMatch);
	}
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// arguments
	string	FILE1;
	string	FILE2; // = perch::DBpath()+"nsfp33a.gz";
	vector<string>	OutHeader;		// empty means not set, output everything; not empty means output only one column to a specific column, supports multi-run
	bool	print_headerRow=true;
	bool	print_unmatched=true;
	bool	unordered_input=false;	// Cannot do fast-join, must do tabix. It is also better for small input files.
	bool	check_ordering1=true;
	bool	ann_snv_only=false;
	bool	ann_indel_only=false;
	bool	replace_old_ann=true;
	bool	preload_data=false;
	bool	check_sorting=false;
	
	// vAnnDel-specific parameters
	data1.InpFormat.set_option(SKIP_NOTES,false);
	data2.InpFormat.set_option(SKIP_NOTES,true);
	data1.InpFormat.set_delimiters("\t");
	data2.InpFormat.set_delimiters("\t");
	data1.InpFormat.comment_sw()="##";
	data2.InpFormat.comment_sw()="##";
	data1.InpFormat.set_titlelines(0);
	data2.InpFormat.set_titlelines(1);
	data2.OutFields.push_back(5,-1);
	set<string>	header_chr = {"#CHROM","CHROM","#CHR","CHR"};
	set<string>	header_pos = {"POS","START","POSITION"};
	set<string>	header_ref = {"REF"};
	set<string>	header_alt = {"ALT"};

	// read arguments
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi) { if (program.arg()[argi]=="--fmt1") { program.arg().erase(program.arg().begin()+argi); data1.InpFormat.read_arguments(program.arg(),argi,true,true); break; } }
	for (size_t argi=1;argi<program.arg().size();++argi) { if (program.arg()[argi]=="--fmt2") { program.arg().erase(program.arg().begin()+argi); data2.InpFormat.read_arguments(program.arg(),argi,true,true); break; } }
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--chr"))			ReadSet(program.arg(),argi,header_chr);
		else if	(str_startsw(program.arg()[argi],"--pos"))			ReadSet(program.arg(),argi,header_pos);
		else if	(str_startsw(program.arg()[argi],"--ref"))			ReadSet(program.arg(),argi,header_ref);
		else if	(str_startsw(program.arg()[argi],"--alt"))			ReadSet(program.arg(),argi,header_alt);
		else if	(str_startsw(program.arg()[argi],"--ann"))			ReadArg(program.arg(),argi,FILE2);
		else if	(str_startsw(program.arg()[argi],"--ms"))			ReadArg(program.arg(),argi,MisStr);
		else if	(str_startsw(program.arg()[argi],"-f2"))		{	data2.OutFields.clear(); ReadArg(program.arg(),argi,data2.OutFields); }
		else if	(str_startsw(program.arg()[argi],"-f"))			{	data2.OutFields.clear(); ReadArg(program.arg(),argi,data2.OutFields); }
		else if	(str_startsw(program.arg()[argi],"--wr"))			ReadSet(program.arg(),argi,OutHeader);
		else if	(str_startsw(program.arg()[argi],"--no-indel"))		ReadArg(program.arg(),argi,ann_snv_only);
		else if	(str_startsw(program.arg()[argi],"--no-snv"))		ReadArg(program.arg(),argi,ann_indel_only);
		else if	(str_startsw(program.arg()[argi],"--snv-only"))		ReadArg(program.arg(),argi,ann_snv_only);
		else if	(str_startsw(program.arg()[argi],"--indel-only"))	ReadArg(program.arg(),argi,ann_indel_only);
		else if	(str_startsw(program.arg()[argi],"--use-tabix"))	ReadArg(program.arg(),argi,unordered_input); // previously --unordered
		else if	(str_startsw(program.arg()[argi],"--check-order"))	ReadArg(program.arg(),argi,check_ordering1);
		else if	(str_startsw(program.arg()[argi],"--replace"))		ReadArg(program.arg(),argi,replace_old_ann);
		else if	(str_startsw(program.arg()[argi],"--remove"))		ReadArg(program.arg(),argi,remove_old_ann);
		else if	(str_startsw(program.arg()[argi],"--add-info"))		ReadArg(program.arg(),argi,AddInf);
		else if	(str_startsw(program.arg()[argi],"--preload"))		ReadArg(program.arg(),argi,preload_data);
		else if	(str_startsw(program.arg()[argi],"--add-ms"))		ReadArg(program.arg(),argi,add_ms);
		else if	(str_startsw(program.arg()[argi],"--check"))		ReadArg(program.arg(),argi,check_sorting);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else if (FILE1.empty()) FILE1=program.arg()[argi];
		else exit_error("Excessive argument "+program.arg()[argi]);
	}

	// show help
	program.help_text_var("_Default_file",FILE2);
	program.help_text_var("_Default_chr",str_of_container(header_chr,',',false));
	program.help_text_var("_Default_pos",str_of_container(header_pos,',',false));
	program.help_text_var("_Default_ref",str_of_container(header_ref,',',false));
	program.help_text_var("_Default_alt",str_of_container(header_alt,',',false));
	program.help_text_var("_Default_unorder",str_YesOrNo(unordered_input));
	program.help_text_var("_Default_check_order",str_YesOrNo(check_ordering1));
	program.help_text_var("_Default_skip_indels",str_YesOrNo(ann_snv_only));
	program.help_text_var("_Default_skip_snvs",str_YesOrNo(ann_indel_only));
	program.help_text_var("_Default_replace",str_YesOrNo(replace_old_ann));
	program.help_text_var("_Default_remove",str_YesOrNo(remove_old_ann));
	program.help_text_var("_Default_add_info",str_YesOrNo(AddInf));
	program.help_text_var("_Default_add_ms",str_YesOrNo(add_ms));
	program.help_text_var("_Default_preload",str_YesOrNo(preload_data));
	program.help_text_var("_Default_f1",data1.OutFields.print());
	program.help_text_var("_Default_f2",data2.OutFields.print());
	program.help_text_var("_Default_wr",str_of_container(OutHeader,',',false));
	program.help_text_var("_Default_ms",MisStr);
	perch::check_arguments();
	
	if (check_sorting)
	{
		string indexA;
		bool header_not_read=true;
		for (Rows_in_File(in1,FILE1,&data1.InpFormat))
		{
			if (in1[0].size()>1 && in1[0][0]=='#' && in1[0][1]=='#') continue;
			if (exist_any(perch::h_col1, in1.contents()))
			{
				header_not_read=false;
				kf1.assign(4,-1);
				for (int i=0;i<in1.NumFields();++i)
				{
					string in1_i_ = boost::to_upper_copy(in1[i]);
					if (kf1[0]<0 && exist_element(header_chr,in1_i_)) { kf1[0]=i; }
					if (kf1[1]<0 && exist_element(header_pos,in1_i_)) { kf1[1]=i; }
					if (kf1[2]<0 && exist_element(header_ref,in1_i_)) { kf1[2]=i; }
					if (kf1[3]<0 && exist_element(header_alt,in1_i_)) { kf1[3]=i; }
				}
				if (*std::min_element(kf1.begin(),kf1.end())<0) exit_error(FILE1+" is missing some columns (Chr,Pos,Alt) or the header row.");
				continue;
			}
			if (header_not_read) exit_error("Header lines missing.");
			if (!read_keys_from(in1,kf1,data1)) exit_error("failed to read key from "+FILE1);
			string index1 = in1[kf1[0]]+","+in1[kf1[1]]+","+in1[kf1[2]]+","+in1[kf1[3]];
			if (dataA.key1!=0 && data1<dataA)
			{
				lns << showe << "input file is out of order between the lines " << indexA << " and " << index1 << fatal;
				return 1;
			}
			indexA=index1;
			dataA.key1 = data1.key1;
			dataA.key2 = data1.key2;
			dataA.key3 = data1.key3;
			dataA.key4 = data1.key4;
		}
		return 0;
	}

	if (FILE2.empty())
	{
		for (Lines_in_File(in1,FILE1,&data1.InpFormat)) program.outf << in1[0] << endl;
		return 0;
	}
	else
	{
		FILE2=perch::find_file(FILE2);
	}
	
	// variables
	int skipIndels = elog.get_token("lines in "+FILE1+" skipped (InDels)");
	int skipSNVs   = elog.get_token("lines in "+FILE1+" skipped (SNVs)");
	int Structural = elog.get_token("lines in "+FILE1+" skipped (structural variation)");
	int multiAltAl = elog.get_token("lines in "+FILE1+" skipped (multiple ALT alleles)");
	int fndNoMatch = elog.get_token("lines in "+FILE1+" have no match in "+FILE2);
	int foundMatch = elog.get_token("lines in "+FILE1+" have  a match in "+FILE2);
	int parse_err1 = elog.get_token("lines in "+FILE1+" failed in reading data; they were joined with missing values.");
	int parse_err2 = elog.get_token("lines in "+FILE2+" failed in reading data; they were ignored.");
	data2.InpFormat.set_field_nums(data2.OutFields, "lines in "+FILE2+" lack the out fields and were skipped.",tfile_format::Continue);
	
	// prepare in2
	bool in2_not_ended=false;
	tabular_file in2(FILE2,&data2.InpFormat);
	for (data2.set_max(); in2.read_r(); in2.next()) // read the first line from in2
	{
		if (in2.is_header()) {
			kf2.assign(4,-1);
			int max_i=-1;
			for (int i=0;i<in2.NumFields();++i)
			{
				string in2_i_ = boost::to_upper_copy(in2[i]);
				if (kf2[0]<0 && exist_element(header_chr,in2_i_)) { kf2[0]=i; if (i>max_i) max_i=i; }
				if (kf2[1]<0 && exist_element(header_pos,in2_i_)) { kf2[1]=i; if (i>max_i) max_i=i; }
				if (kf2[2]<0 && exist_element(header_ref,in2_i_)) { kf2[2]=i; if (i>max_i) max_i=i; }
				if (kf2[3]<0 && exist_element(header_alt,in2_i_)) { kf2[3]=i; if (i>max_i) max_i=i; }
			}
			if (*std::min_element(kf2.begin(),kf2.end())<0) exit_error(FILE2+" is missing some columns (Chr,Pos,Ref,Alt) or the header row.");
			data2.InpFormat.set_field_nums(max_i+1,"lines in "+FILE2+" lack the key fields and were skipped.",tfile_format::Continue);
			header2=in2.contents();
			continue; }
		if (read_keys_from(in2,kf2,data2)) { in2_not_ended=true; break; }
		else { elog.add(parse_err2); }
	}
	if (in2_not_ended==false) exit_error("all lines are invalid in "+FILE2);
	if (header2.empty()) exit_error("No header row in "+FILE2); // added 2016-6-21
	if (!OutHeader.empty())
	{
		if (data2.OutFields.size()!=OutHeader.size()) exit_error("the number of --wr headers not matching the number of -f fields");
		for (size_t i=0;i<data2.OutFields.size();++i)
			header2[data2.OutFields[i]]=OutHeader[i];
	}
	else
	{
		for (size_t i=0;i<data2.OutFields.size();++i)
			OutHeader.push_back(header2[data2.OutFields[i]]);
	}
	
	// join
	typedef std::tuple<int,int,string,string> ID_t;
	map<ID_t,string> annotation; // index => MaxAF
	if (preload_data)
	{
		if (data2.OutFields.size()!=1) exit_error("--preload only works for one annotation");
		
		// prepare annotation		
		lns<<showl<<"Read data from "<<FILE2<<flush_logger;
		do
		{
			annotation[std::make_tuple(data2.key1, data2.key2, data2.key3, data2.key4)] = in2[data2.OutFields[0]];
		}
		while (in2.read_r() && read_keys_from(in2,kf2,data2));
	}

	bool header_not_read=true;
	tabular_file in1(FILE1,&data1.InpFormat);
	for (; in1.read_r(); in1.next())
	{
		if (in1[0].size()>1 && in1[0][0]=='#' && in1[0][1]=='#') // meta data
		{
			in1.clear_nf();
			bool to_remove=false;
			for (auto &h:OutHeader)
				if (str_has(in1[0],"##INFO=<ID="+h+","))
				{
					if (!AddInf) exit_error(h+" is already in INFO, but you didn't use --add-info");
					to_remove=true;
					break;
				}
			if (to_remove) continue;
			print_container(in1.contents(),program.outf,' ',true);
			continue;
		}
		if (exist_any(perch::h_col1, in1.contents()))
		{
			header_not_read=false;
			kf1.assign(4,-1);
			int max_i=-1;
			for (int i=0;i<in1.NumFields();++i)
			{
				string in1_i_ = boost::to_upper_copy(in1[i]);
				if (kf1[0]<0 && exist_element(header_chr,in1_i_)) { kf1[0]=i; if (i>max_i) max_i=i; }
				if (kf1[1]<0 && exist_element(header_pos,in1_i_)) { kf1[1]=i; if (i>max_i) max_i=i; }
				if (kf1[2]<0 && exist_element(header_ref,in1_i_)) { kf1[2]=i; if (i>max_i) max_i=i; }
				if (kf1[3]<0 && exist_element(header_alt,in1_i_)) { kf1[3]=i; if (i>max_i) max_i=i; }
				if (in1_i_=="INFO") ColInf=i;
			}
			for (auto &h:OutHeader)
			{
				for (int i=0;i<in1.NumFields();++i)
					if (in1[i]==h) { data1.OutFields.push_back(i+1); break; }
			}
			if (!data1.OutFields.no_input())
			{
				if (data1.OutFields.size()!=data2.OutFields.size())
				{
					exit_error("some but not all columns are missing");
				}
				else
				{
					if (AddInf) exit_error("Column(s) already exist in "+FILE1+", which conflicts with --add-info.");
				}
			}
			if (AddInf && ColInf==-1) exit_error("The INFO column does not exist");
			if (*std::min_element(kf1.begin(),kf1.end())<0) exit_error(FILE1+" is missing some columns (Chr,Pos,Alt) or the header row.");
			data1.InpFormat.set_field_nums(max_i+1,"lines in "+FILE1+" lack the key fields and were skipped.",tfile_format::Continue);
			if (print_headerRow)
			{
				if (AddInf)
				{
					for (auto &h:OutHeader)
						program.outf<<"##INFO=<ID="<<h<<",Number=1,Type=Float,Description=\""<<h<<" annotated by vAnnDel.\">"<<endl;
					in1.write_r(program.outf,true);
				}
				else
					write_matched(in1,header2,program.outf);
			}
			continue;
		}
		if (header_not_read) exit_error("Header lines missing.");
		if		(str_has(in1[kf1[3]],"<"))
		{
			if (print_unmatched) write_uniqFILE1(in1,MisStr,program.outf,add_ms);
			elog.add(Structural);
		}
		else if (str_has(in1[kf1[3]],","))
		{
			if (print_unmatched) write_uniqFILE1(in1,MisStr,program.outf,add_ms);
			elog.add(multiAltAl);
		}
		else if (ann_snv_only && (in1[kf1[2]].size()!=1 || in1[kf1[3]].size()!=1))
		{
			if (print_unmatched) write_uniqFILE1(in1,MisStr,program.outf,add_ms);
			elog.add(skipIndels);
		}
		else if (ann_indel_only && (in1[kf1[2]].size()==1 && in1[kf1[3]].size()==1))
		{
			if (print_unmatched) write_uniqFILE1(in1,MisStr,program.outf,add_ms);
			elog.add(skipSNVs);
		}
		else if (!replace_old_ann && already_annotated(in1))
		{
			print_container(in1.contents(),program.outf,DLMTR,true);
		}
		else if (read_keys_from(in1,kf1,data1))
		{
			if (preload_data)
			{
				map<ID_t,string>::iterator it = annotation.find(std::make_tuple (data1.key1, data1.key2, data1.key3, data1.key4));
				if (it == annotation.end())
				{
					if (print_unmatched) write_uniqFILE1(in1,MisStr,program.outf,add_ms);
					elog.add(fndNoMatch);
				}
				else
				{
					write_uniqFILE1(in1,it->second,program.outf,true);
					elog.add(foundMatch);
				}
			}
			else if (unordered_input)
			{
				annotate(in1, FILE2, print_unmatched, MisStr, foundMatch, fndNoMatch);
			}
			else
			{
				if (check_ordering1)
				{
					if (data1 < dataA)
					{
						// exit_error("Input is unsorted, please sort or use the option --unsorted.");
						annotate(in1, FILE2, print_unmatched, MisStr, foundMatch, fndNoMatch);
						continue;
					}
					dataA.key1 = data1.key1;
					dataA.key2 = data1.key2;
					dataA.key3 = data1.key3;
					dataA.key4 = data1.key4;
				}
				
				if (data2 < data1)
				{
					for (; in2_not_ended; in2.next())
					{
						if (in2.read_r())
						{
							if (read_keys_from(in2,kf2,data2))
							{
								if (data2 < data1) continue;
								break;
							}
							else { elog.add(parse_err2); }
						}
						else
						{
							data2.set_max();
							in2_not_ended=false;
							break;
						}
					}
				}
				if (data1 < data2)	// unmatched
				{
					if (print_unmatched) write_uniqFILE1(in1,MisStr,program.outf,add_ms);
					elog.add(fndNoMatch);
				}
				else				// matched
				{
					write_matched(in1,in2.contents(),program.outf);
					elog.add(foundMatch);
				}
			}
		}
		else	// FILE1 bad keys
		{
			if (print_unmatched) write_uniqFILE1(in1,MisStr,program.outf,add_ms);
			elog.add(parse_err1);
		}
	}

	return 0;
}
