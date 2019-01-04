#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_math.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	vector<string>	inputs;
	set<string>		must_have;
	size_t			num_genes=250;
	
	// handle program options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--num-genes"))	ReadArg(program.arg(),argi,num_genes);
		else if	(str_startsw(program.arg()[argi],"--must-have"))	ReadSet(program.arg(),argi,must_have);
		else inputs.push_back(program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_num_genes",itos(num_genes));
	program.help_text_var("_Default_must_have",str_of_container(must_have,',',false));
	perch::check_arguments();
	
	// check error
	if (must_have.size()>=num_genes) exit_error("--must-have included sufficient genes for --num-genes already");

	// read
	vector< vector<string> > gene_list;
	for (auto &x:inputs)
	{
		gene_list.push_back(vector<string>());
		vector<string>& this_list=gene_list.back();
		for (Rows_in_File(in,x,1))
		{
			if (in[0]==perch::h_symb) continue;
			this_list.push_back(in[0]);
		}
	}
	
	// run
	vector<string> chosen;
	for (auto &x:must_have) chosen.push_back(x);
	for (size_t i=0;;++i)
	{
		bool done=true;
		for (size_t j=0;j<gene_list.size();++j)
			if (i<gene_list[j].size())
			{
				if (!exist_element(must_have,gene_list[j][i]))
				{
					chosen.push_back(gene_list[j][i]);
					must_have.insert(gene_list[j][i]);
				}
				done=false;
			}
		if (must_have.size()>=num_genes)
		{
			lns<<showl<<"read "<<must_have.size()<<" genes by picking the top "<<i+1<<" from each list."<<flush_logger;
			done=true;
		}
		if (done) break;
	}
	
	// out
	for (auto &x:chosen) program.outf<<x<<endl;
	return 0;
}
