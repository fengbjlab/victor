#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_genepi.hpp>
#include "vAAA_sv.hpp"

namespace AAA_SV
{
	using namespace std;

	bool del_only=true;
	const int ABS_MAX_CP = 5;
	const int MAX_CNV_CP = (del_only?2:ABS_MAX_CP) ;
	
	map<int, map<string, map<string, vector<pair<int,int> > > > >	sv_data_wiSymb[ABS_MAX_CP];// sv_data_wiSymb[copy][chr_num][symbol][SeqID][#]=<start,end>
	map<string, map<int, vector<pair<int,int> > > > 				sv_data_woSymb[ABS_MAX_CP];// sv_data_woSymb[copy][SeqID][chr_num][#]=<start,end>

	// to do
	// 1) don't know whether Start is 0-based in bed
	// 2) handle multiple gene symbols in bed
	// 3) convert chrX chrY copy numbers in bed
	
	void sv_read_wiSymb(const string& filename)
	{
		if (str_endsw(filename,".bed")||str_endsw(filename,".bed.gz"))
		{
			int ColChr=0;
			int ColBeg=1;
			int ColEnd=2;
			int ColSym=14;
			int ColSID=3;
			int ColCpy=5;
			for (Rows_in_File(in,filename,7))
			{
				int	chr_num;	if (!genepi::read_chr_num(in[ColChr],chr_num))	exit_error("Failed to read "+in[ColChr]+" as a chromosome.");
				int	bp_beg;		if (!read_val_ge(in[ColBeg],bp_beg,1))			exit_error("Failed to read "+in[ColBeg]+" as a position in basepairs.");
				int	bp_end;		if (!read_val_ge(in[ColEnd],bp_end,1))			exit_error("Failed to read "+in[ColEnd]+" as a position in basepairs.");
				int	copy;		if (!read_val_ge(in[ColCpy],copy,0))			exit_error("Failed to read "+in[ColCpy]+" as a copy number.");
				if (copy>=MAX_CNV_CP) continue;
				if (del_only&&copy>1) continue;
				sv_data_wiSymb[copy][chr_num][in[ColSym]][in[ColSID]].push_back(make_pair(bp_beg,bp_end));
			}
		}
		else
		{
			int ColChr=0;
			int ColBeg=1;
			int ColEnd=2;
			int ColSym=4;
			int ColSID=5;
			int ColCpy=6;
			for (Rows_in_File(in,filename,7))
			{
				int	chr_num;	if (!genepi::read_chr_num(in[ColChr],chr_num))	exit_error("Failed to read "+in[ColChr]+" as a chromosome.");
				int	bp_beg;		if (!read_val_ge(in[ColBeg],bp_beg,1))			exit_error("Failed to read "+in[ColBeg]+" as a position in basepairs.");
				int	bp_end;		if (!read_val_ge(in[ColEnd],bp_end,1))			exit_error("Failed to read "+in[ColEnd]+" as a position in basepairs.");
				int	copy;		if (!read_val_ge(in[ColCpy],copy,0))			exit_error("Failed to read "+in[ColCpy]+" as a copy number.");
				if (copy>=MAX_CNV_CP) continue;
				if (del_only&&copy>1) continue;
				sv_data_wiSymb[copy][chr_num][in[ColSym]][in[ColSID]].push_back(make_pair(bp_beg,bp_end));
			}
		}
	}
	
	void sv_read_wiSymb(const vector<string>& filenames)
	{
		for (auto &f:filenames) sv_read_wiSymb(f);
	}

	void sv_read_woSymb(const string& filename)
	{
		if (str_endsw(filename,".bed")||str_endsw(filename,".bed.gz"))
		{
			int ColChr=0;
			int ColBeg=1;
			int ColEnd=2;
			int ColSID=3;
			int ColCpy=5;
			for (Rows_in_File(in,filename,7))
			{
				int	chr_num;	if (!genepi::read_chr_num(in[ColChr],chr_num))	exit_error("Failed to read "+in[ColChr]+" as a chromosome.");
				int	bp_beg;		if (!read_val_ge(in[ColBeg],bp_beg,1))			exit_error("Failed to read "+in[ColBeg]+" as a position in basepairs.");
				int	bp_end;		if (!read_val_ge(in[ColEnd],bp_end,1))			exit_error("Failed to read "+in[ColEnd]+" as a position in basepairs.");
				int	copy;		if (!read_val_ge(in[ColCpy],copy,0))			exit_error("Failed to read "+in[ColCpy]+" as a copy number.");
				if (copy>=MAX_CNV_CP) continue;
				if (del_only&&copy>1) continue;
				sv_data_woSymb[copy][in[ColSID]][chr_num].push_back(make_pair(bp_beg,bp_end));
			}
		}
		else
		{
			int ColChr=0;
			int ColBeg=1;
			int ColEnd=2;
			int ColSID=5;
			int ColCpy=6;
			for (Rows_in_File(in,filename,7))
			{
				int	chr_num;	if (!genepi::read_chr_num(in[ColChr],chr_num))	exit_error("Failed to read "+in[ColChr]+" as a chromosome.");
				int	bp_beg;		if (!read_val_ge(in[ColBeg],bp_beg,1))			exit_error("Failed to read "+in[ColBeg]+" as a position in basepairs.");
				int	bp_end;		if (!read_val_ge(in[ColEnd],bp_end,1))			exit_error("Failed to read "+in[ColEnd]+" as a position in basepairs.");
				int	copy;		if (!read_val_ge(in[ColCpy],copy,0))			exit_error("Failed to read "+in[ColCpy]+" as a copy number.");
				if (copy>=MAX_CNV_CP) continue;
				if (del_only&&copy>1) continue;
				sv_data_woSymb[copy][in[ColSID]][chr_num].push_back(make_pair(bp_beg,bp_end));
			}
		}
	}
	
	void sv_read_woSymb(const vector<string>& filenames)
	{
		for (auto &f:filenames) sv_read_woSymb(f);
	}
	
	int sv_copy(const int& chr_num, const string& symbol, const string& SeqID, const int& pos)
	{
		for (int i=0;i<MAX_CNV_CP;++i)
		{
			auto it1 = sv_data_wiSymb[i].find(chr_num);
			if (it1!=sv_data_wiSymb[i].end())
			{
				auto it2 = it1->second.find(symbol);
				if (it2!=it1->second.end())
				{
					auto it3 = it2->second.find(SeqID);
					if (it3!=it2->second.end())
					{
						for (auto &segment:it3->second)
						{
							int& bp_beg=segment.first;
							int& bp_end=segment.second;
							if (pos>=bp_beg && pos<=bp_end) return i;
						}
					}
				}
			}
		}
		return 2;
	}
	
	int sv_copy(const int& chr_num, const string& symbol, const string& SeqID)
	{
		int min_cp= std::numeric_limits<int>::max();
		int max_cp=-std::numeric_limits<int>::max();
		for (int i=0;i<MAX_CNV_CP;++i)
		{
			auto it1 = sv_data_wiSymb[i].find(chr_num);
			if (it1!=sv_data_wiSymb[i].end())
			{
				auto it2 = it1->second.find(symbol);
				if (it2!=it1->second.end())
				{
					auto it3 = it2->second.find(SeqID);
					if (it3!=it2->second.end())
					{
						if (i<min_cp) min_cp=i;
						if (i>max_cp) max_cp=i;
					}
				}
			}
		}
		if (min_cp<2) return min_cp;
		if (max_cp>2) return max_cp;
		return 2;
	}

	int sv_exist(const int& chr_num, const string& symbol, const vector<string>& SeqIDs)
	{
		int count = 0;
		for (int i=0;i<MAX_CNV_CP;++i)
		{
			auto it1 = sv_data_wiSymb[i].find(chr_num);
			if (it1!=sv_data_wiSymb[i].end())
			{
				auto it2 = it1->second.find(symbol);
				if (it2!=it1->second.end())
				{
					for (auto &SeqID:SeqIDs)
					{
						auto it3 = it2->second.find(SeqID);
						if (it3!=it2->second.end()) ++count;
					}
				}
			}
		}
		return count;
	}

};
