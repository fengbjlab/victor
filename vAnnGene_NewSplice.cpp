
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"
#include "vAnnGene_mes5.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// input
	string	mes5in;								// input genotype file in VCF format
	string	vcf_in;
	
	// handle program arguments
	string version_number = "1.2beta";
	string version_string = "Version "+version_number+" Build " + macro_date_to_boost(string(__DATE__));
	program.trademark = "("+version_string+")";
#ifdef PERCH_DISTRIBUTION
	program.enable_option("--no-web");
	program.set_check_url("http://www.fengbj-laboratory.org/download/victor_version","https://www.dropbox.com/s/49snwy5odsmfzgx/victor_version", version_string);
#endif

	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"--mes5"))			ReadArg(program.arg(),argi,mes5in);
		else if (vcf_in.empty()) vcf_in=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	perch::add_help_text_var();
	program.help_text_var("_Default_mes5",mes5in);
	program.check_help_request();
	
	MaxEntScan_score5ss mes5;
	mes5.read(mes5in);
	cout<<mes5.score("AAAAAAAGC")<<endl; // -9.06
	cout<<mes5.score("AAAAAAGCT")<<endl; // -16.44
	cout<<mes5.score("AAAAAGCTT")<<endl; // -19.12
	cout<<mes5.max("AAAAAAAGCTT")<<endl; // -9.06
	return 0;
}
