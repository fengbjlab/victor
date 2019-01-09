#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	string			vcf_in;
	string			ped_in;
	string			id_del=":";
	string			DoWhat="perch2plink"; // perch2plink linkage2vcf
	
	// handle program options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--ped"))	ReadArg(program.arg(),argi,ped_in);
		else if	(str_startsw(program.arg()[argi],"--geno"))	ReadArg(program.arg(),argi,vcf_in);
		else if (str_startsw(program.arg()[argi],"--do"))	ReadArg(program.arg(),argi,DoWhat);
		else exit_error("unknown option "+program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_ped",ped_in);
	program.help_text_var("_Default_geno",vcf_in);
	program.help_text_var("_Default_do",DoWhat);
	perch::check_arguments();

	if (DoWhat=="perch2plink")
	{
		// check arguments
		if (ped_in.empty()) exit_error("Please set --ped.");
		if (vcf_in.empty()) exit_error("Please set --vcf.");
		
		// read vcf
		tfile_format format;
		format.set_delimiters("\t");
		format.set_option(SKIP_NOTES,false);
		set<string> sequenced_ID;
		for (Rows_in_File(in, vcf_in, &format))
		{
			if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#') // meta data
			{
				in.clear_nf();
				continue;
			}
			if (exist_any(perch::h_col1, in.contents()))
			{
				for (int i=9;i<in.NumFields();++i)
					sequenced_ID.insert(in[i]);
				break;
			}
			exit_error("cannot find the header line in "+vcf_in);
		}
		
		// read ped
		set<string>	h_sid = { "SeqID", "SequencingID", "SequenceID", "ExomeID", "GenomeID", "ExperimentID", "WES_ID", "WGS_ID", "NGS_ID" };
		map<string,string> idx2sid; // index (PID:IID) to SeqID
		vector< vector<string> > data;
		for (Rows_in_File(in, ped_in, 7))
		{
			// header
			if (in.RowNumber()==0)
			{
				if (!exist_any(h_sid,in.contents())) exit_error(ped_in+" is not in PERCH format.");
				continue;
			}
			
			// SeqID
			string SeqID = in[6];
			bool is_proband = false;
			perch::read_SeqID(SeqID,is_proband);
			if (SeqID.empty()) SeqID="0";
			if (exist_element(perch::rm_ind,SeqID)) SeqID="0";
			if (!exist_element(sequenced_ID,SeqID)) SeqID="0";
			if (SeqID!="0") idx2sid[in[0]+id_del+in[1]]=SeqID;
			
			// sex and aff
			perch::read_sex(in[4]);
			perch::read_aff(in[5]);
			data.push_back(in.contents());
		}
		
		// write pedigree file
		for (auto &in:data)
		{
			string iid = in[0]+id_del+in[1];
			string dad = in[0]+id_del+in[2];
			string mom = in[0]+id_del+in[3];
			if (exist_element(idx2sid,iid)) in[1]=idx2sid[iid];
			if (exist_element(idx2sid,dad)) in[2]=idx2sid[dad];
			if (exist_element(idx2sid,mom)) in[3]=idx2sid[mom];
			print_container_head(in,6,program.outf,DLMTR,true);
		}
	}
	else if (DoWhat=="linkage2vcf")
	{
		// check arguments
		if (ped_in.empty()) exit_error("Please set --ped.");
		
		exit_error("this function is not done yet");
	}
	else exit_error("unknown run mode");
	
	return 0;
}
