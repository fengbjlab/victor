#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_genepi.hpp>
#include "victor_par.hpp"
#include "vAAA_if.hpp"

namespace AAA_IF
{
	using namespace std;

	map<string, set<string> > panel_ifg;
	set<string> panel_rec;
	map<string, map<int,vector<string> > > partial;

	void read_panels()
	{
		for (Rows_in_File(in,perch::find_file("panel_incidental"),2)) panel_ifg[in[0]].insert(in[1]);
		for (Rows_in_File(in,perch::find_file("panel_recessive"),1)) panel_rec.insert(in[0]);
	}
	
	void to_rpt(const string& symbol, int A, int SeqID, string line, string& output)
	{
		if (exist_element(panel_ifg,symbol))
		{
			line=line+"="+str_of_container(panel_ifg[symbol],';',true);
			if (exist_element(panel_rec,symbol))
			{
				if (A==2) output+=line;
				else if (A==1)
				{
					partial[symbol][SeqID].push_back(line);
					size_t size=partial[symbol][SeqID].size();
					if (size==2) output+=partial[symbol][SeqID][0];
					if (size>=2) output+=line;
				}
			}
			else
			{
				if (A>0) output+=line;
			}
		}
	}
};
