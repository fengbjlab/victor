// shapeit maintains the number and order variants and samples, which makes the programing easier.
// shapeit does not change genotype, but will fill in missing genotypes. I'd like to keep the missingness.

#include <tft/libfbj_base.hpp>
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

	// other parameters
	// parameters
	bool			nohead=false;
	string			vcf_in;						// genotype file in VCF format
	string			max_in;						// SHAPEIT max output
	
	// handle program arguments
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"--no-head"))		ReadArg(program.arg(),argi,nohead);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else if (max_in.empty()) max_in=program.arg()[argi];
		else if (vcf_in.empty()) vcf_in=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	perch::check_arguments();

	// check errors
	if (max_in.empty()) exit_error("please specify the SHAPEIT max output file");

	// read VCF and do linkage analysis line by line
	tfile_format format_vcf, format_max;
	format_vcf.set_delimiters("\t");
	format_vcf.set_option(SKIP_NOTES,false);
	bool header_not_read = true;
	tabular_file phased(max_in,&format_max);
	for (Rows_in_File(in, vcf_in, &format_vcf))
	{
		// read header
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			if (!header_not_read) { continue; }
			if (!nohead) print_container(in.contents(),program.outf,' ',true);
			in.clear_nf();
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			if (!header_not_read) { continue; }
			header_not_read=false;
			if (!nohead) print_container(in.contents(),program.outf,DLMTR,true);
			continue;
		}
		if (header_not_read) exit_error("Header lines missing.");
		
		if (!phased.read_r()) exit_error("not enough rows in "+max_in);
		if (phased.NumFields()!=(in.NumFields()-9)*2+5) exit_error("wrong number of columns in "+max_in);
		for (int i=9,j=4;i<in.NumFields();++i)
		{
			char a1 = in[i][0];
			string& h1 = phased[++j];
			string& h2 = phased[++j];
			if (a1!='.')	in[i] = h1+"|"+h2 + trim_before_find(in[i],":") ;
			else			in[i][1]='|';
		}
		phased.next();
		print_container(in.contents(),program.outf,DLMTR,true);
	}
	return 0;
}
