/*
 --no-head B       Do not output header lines that begin with a #
*/
#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"

using namespace std;

int output(vector< vector<string> >& previous)
{
	if (previous.empty()) return 0;
	if (previous.size()==1) { print_container(previous[0],program.outf,DLMTR,true); return 0; }
	vector<string> alleles;		// MNP alleles
	map<string,int> coordinate;	// MNP allele -> ID (0+, 0 is REF+REF+..)
	string ref; for (size_t	i=0;i<previous.size();++i) ref+=previous[i][3];
	alleles.push_back(ref); coordinate[ref]=0;
	size_t columns = previous[0].size();
	size_t lines = previous.size();
	vector<int> allele1(columns);
	vector<int> allele2(columns);
	for (size_t i=9;i<columns;++i)
	{
		bool missing=false;
		string h1,h2;
		for (size_t j=0;j<lines;++j)
		{
			if (previous[j][i][0]=='.') { missing=true; break; }
			if (previous[j][i][0]=='0') h1+=previous[j][3]; else if (previous[j][i][0]=='1') h1+=previous[j][4]; else exit_error("unlikely");
			if (previous[j][i][2]=='0') h2+=previous[j][3]; else if (previous[j][i][2]=='1') h2+=previous[j][4]; else exit_error("unlikely");
		}
		if (missing)
		{
			allele1[i]=-1;
			allele2[i]=-1;
		}
		else
		{
			if (!exist_element(coordinate,h1)) { coordinate[h1]=alleles.size(); alleles.push_back(h1); }
			if (!exist_element(coordinate,h2)) { coordinate[h2]=alleles.size(); alleles.push_back(h2); }
			allele1[i]=coordinate[h1];
			allele2[i]=coordinate[h2];
		}
	}
	if (alleles.size()==2) // not robust to there're only REF1-ALT2 and REF2-ALT1 in data because REF1-REF2 is already added without seeing the data
	{
		for (int alt=1;alt<(int)alleles.size();++alt)
		{
			program.outf<<previous[0][0]<<'\t'<<previous[0][1]<<"\t.\t"<<ref<<'\t'<<alleles[alt]<<"\t.\t.\tMNP_by_VANNER\tGT";
			for (size_t i=9;i<columns;++i)
			{
				program.outf<<'\t'; if (allele1[i]==-1) program.outf<<'.'; else if (allele1[i]==0) program.outf<<'0'; else if (allele1[i]==alt) program.outf<<'1'; else program.outf<<'a';
				program.outf<<'|';  if (allele2[i]==-1) program.outf<<'.'; else if (allele2[i]==0) program.outf<<'0'; else if (allele2[i]==alt) program.outf<<'1'; else program.outf<<'a';
			}
			program.outf<<endl;
		}
		return 1;
	}
	else
	{
		for (size_t j=0;j<lines;++j)
			print_container(previous[j],program.outf,DLMTR,true);
	}
	return 0;
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	// parameters
	bool			nohead=false;
	string			vcf_in;						// genotype file in VCF format
	
	// handle program arguments
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"--no-head"))		ReadArg(program.arg(),argi,nohead);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else if (vcf_in.empty()) vcf_in=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	perch::check_arguments();
	
	// read VCF
	vector< vector<string> > previous;
	tfile_format format_vcf;
	format_vcf.set_delimiters("\t");
	format_vcf.set_option(SKIP_NOTES,false);
	bool header_not_read = true;
	int num_mnp=0;
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
		bool phased = (in[9][1]=='|');
		if (previous.empty())
		{
			if (phased) previous.push_back(in.contents());
			else print_container(in.contents(),program.outf,DLMTR,true);
		}
		else
		{
			vector<string>& last = previous.back();
			int last_bp, this_bp;
			read_val(last[1],last_bp);
			read_val(in[1],this_bp);
			if (this_bp==last_bp+(int)last[3].size() && last[0]==in[0] && phased)
			{
				previous.push_back(in.contents());
			}
			else
			{
				num_mnp+=output(previous);
				previous.clear();
				if (phased) previous.push_back(in.contents());
				else print_container(in.contents(),program.outf,DLMTR,true);
			}
		}
	}
	num_mnp+=output(previous);
	lns<<showl<<num_mnp<<" MNPs merged."<<flush_logger;
	return 0;
}
