#include <tft/libfbj_file.hpp>
#include <tft/libfbj_base.hpp>
#include "vAnnGene_mes5.hpp"

using namespace std;

void convert(const string& seq, vector<int>& vv)
{
	vv.assign(seq.size(),0);
	for (size_t i=0;i<seq.size();++i)
	{
		if		(seq[i]=='C') vv[i]=1;
		else if (seq[i]=='G') vv[i]=2;
		else if (seq[i]=='T') vv[i]=3;
	}
}

int coordinate(vector<int>& vv, int begin, int end)
{
	int res=0;
	for (int i=begin, j=end-begin-1; i<end; ++i, --j)
	{
		if (vv[i]) res += pow(4,j) * vv[i];
	}
	return res;
}

void	MaxEntScan_score5ss::read(const std::string& filename)
{
	scores.assign(262144,0);
	for (Rows_in_File(in,filename,2))
	{
		double val;
		if (!read_val(in[1],val)) exit_error("Failed to read MaxEntScan score5ss value "+in[1]+" in "+filename);
		scores[in.RowNumber()]=val;
	}
}

double MaxEntScan_score5ss::score(const std::string& seq)
{
	vector<int> vv;
	convert(seq,vv);
	int c = coordinate(vv,0,vv.size());
	return scores[c];
}

double MaxEntScan_score5ss::max(const std::string& seq)
{
	static const int nn=9;
	vector<int> vv;
	convert(seq,vv);
	double maxscore = -std::numeric_limits<double>::max();
	for (int begin=0, end=nn; end<=vv.size(); ++begin, ++end)
	{
		int c = coordinate(vv,begin,end);
		double s = scores[c];
		if (s > maxscore) maxscore = s;
	}
	return maxscore;
}

void test ()
{
	MaxEntScan_score5ss mes5;
	mes5.read("/Users/bf/work/software/MaxEntScan/all9mer5ss.out");
	cout<<mes5.score("AAAAAAAGC")<<endl; // -9.06
	cout<<mes5.score("AAAAAAGCT")<<endl; // -16.44
	cout<<mes5.score("AAAAAGCTT")<<endl; // -19.12
	cout<<mes5.max("AAAAAAAGCTT")<<endl; // -9.06
}
