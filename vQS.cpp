#include <tft/libfbj_file.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);
	string method="kimber"; // kimber 3.09sd. Default is kimber for historical reasons, 3.09sd works closely or better.
	
	// other parameters
	vector<string> inputs;
	double lb=-6;
	double ub=0.5;
	
	// handle program options
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--lb"))		ReadArg(program.arg(),argi,lb);
		else if	(str_startsw(program.arg()[argi],"--ub"))		ReadArg(program.arg(),argi,ub);
		else if	(str_startsw(program.arg()[argi],"--method"))	ReadArg(program.arg(),argi,method);
		else inputs.push_back(program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_lb",ftos(lb));
	program.help_text_var("_Default_ub",ftos(ub));
	program.help_text_var("_Default_method",method);
	perch::check_arguments();

	// check errors
	boost::to_lower(method);
	if (method!="kimber"&&!str_endsw(method,"sd")) exit_error("--method must be kimber or #sd");
	lns<<showl<<"The method to calculate VQSLOD cutoff is "<<method<<flush_logger;
	
	// read file
	int NumLOD=0;								// for Ti/Tv.  total number of coding SNVs
	typedef double rank_t; 						// for Ti/Tv.  must be floating point because of ties
	map<double,rank_t> 			val_to_rank;	// for Ti/Tv.  VQSLOD => rank number
	map<double,pair<int,int> > 	titv;			// for Ti/Tv.  VQSLOD => number of Ti,Tv whose VQSLOD equals to the index value
	Values<double> 				snv, idl;		// for Kimber. VQSLOD
	for (Rows_in_File(in,inputs,4))
	{
		if (in[0]!="known") continue;
		double v;
		if (!read_val_noNaN(in[2],v)) exit_error("failed to read VQSLOD from "+in.FileName());
		if 		(in[1]=="SNV") snv.push_back(v);
		else if	(in[1]=="InDel") idl.push_back(v);
		else exit_error("failed to read variant type (SNV/InDel) from "+in.FileName());
		if (in[3]==".") continue;
		++NumLOD;
		++val_to_rank[v];
		if (in[3]=="coding_ti") { if (!exist_element(titv,v)) titv.insert(make_pair(v,pair<int,int>(1,0))); else titv[v].first++; }
		if (in[3]=="coding_tv") { if (!exist_element(titv,v)) titv.insert(make_pair(v,pair<int,int>(0,1))); else titv[v].second++;}
	}
	lns<<showl<<"There are  "<<snv.size()<<" known SNVs and "<<idl.size()<<" known InDels."<<flush_logger;
	double cutoff_snv = std::numeric_limits<double>::signaling_NaN();
	double cutoff_idl = std::numeric_limits<double>::signaling_NaN();

	// determine VQSLOD cutoff by Kimber_LF
	if (method=="kimber")
	{
		cutoff_snv = snv.get(STAT::Kimber_LF);
		cutoff_idl = idl.get(STAT::Kimber_LF);
	}
	
	// determine VQSLOD cutoff by SD
	if (str_endsw(method,"sd"))
	{
		double nosd;
		if (!read_val_noNaN(substr_before_find(method,"sd"),nosd)) exit_error("failed to read number from "+method);
		if (nosd<0) exit_error("For --method=#sd, # msut be >=0.");
		double median_snv=snv.get(STAT::MEDIAN);
		double median_idl=idl.get(STAT::MEDIAN);
		Values<double> snv2, idl2;
		for (size_t i=0;i<snv.size();++i) { if (snv[i]<median_snv) { snv2.push_back(snv[i]); snv2.push_back(median_snv+median_snv-snv[i]); } else if (snv[i]==median_snv) snv2.push_back(snv[i]); }
		for (size_t i=0;i<idl.size();++i) { if (idl[i]<median_idl) { idl2.push_back(idl[i]); idl2.push_back(median_idl+median_idl-idl[i]); } else if (idl[i]==median_idl) idl2.push_back(idl[i]); }
		double sd_snv = snv2.get(STAT::SPL_SD);
		double sd_idl = idl2.get(STAT::SPL_SD);
		cutoff_snv = median_snv-nosd*sd_snv; // 3.09 corresponds to p=0.001 one-tail
		cutoff_idl = median_idl-nosd*sd_idl;
		//cerr<<"# SNV: original N="<<snv.size()<<" N="<<snv2.size()<<" 1.96sd="<<median_snv-1.96*sd_snv<<" 3.09sd="<<median_snv-3.09*sd_snv<<" 4sd="<<median_snv-4*sd_snv<<" 5sd="<<median_snv-5*sd_snv<<endl;
		//cerr<<"# IND: original N="<<idl.size()<<" N="<<idl2.size()<<" 1.96sd="<<median_idl-1.96*sd_idl<<" 3.09sd="<<median_idl-3.09*sd_idl<<" 4sd="<<median_idl-4*sd_idl<<" 5sd="<<median_idl-5*sd_idl<<endl;
	}

	// output
	if (cutoff_snv<lb || cutoff_snv>ub) lns<<showw<<"# The calculated VQSLOD cutoff for SNV ("<<cutoff_snv<<") is out of bound ["<<lb<<","<<ub<<"]."<<flush_logger;
	else if (!std::isnan(cutoff_snv))	program.outf<<"--filt-vqs-snv="<<cutoff_snv<<endl;
	if (cutoff_idl<lb || cutoff_idl>ub) lns<<showw<<"# The calculated VQSLOD cutoff for InDel ("<<cutoff_idl<<") is out of bound ["<<lb<<","<<ub<<"]."<<flush_logger;
	else if (!std::isnan(cutoff_idl))	program.outf<<"--filt-vqs-indel="<<cutoff_idl<<endl;
	return 0;
	
	// determine VQSLOD cutoff by VQSLOD_RANK - Ti/Tv plot. Not working.
	// rank data
	rank_t total=0;
	for (each_element_rev(val_to_rank,it))
	{
		rank_t n = it->second;
		it->second = total + (1+n)/2. ;
		total += n;
	}
	// make Ti/Tv curve (X is VQSLOD ranks in reverse order, Y is Ti/Tv of varaints whose VQSLOD rank is ~x)
	map<double,double> rank_to_ratio;
	if (perch::_Debug) program.outf<<"Rank\tTi/Tv\n";
	double width=0.05;
	double steps=0.05;
	for (double mid=steps;mid<=(100-steps);mid+=steps)
	{
		double total_ti=0;
		double total_tv=0;
		for (auto &x:titv)
		{
			double pct = val_to_rank[x.first]*100.0/NumLOD;
			if (pct<(mid-width) || pct>(mid+width)) continue;
			total_ti+=x.second.first;
			total_tv+=x.second.second;
		}
		if (total_tv==0) continue;
		double ratio=total_ti/total_tv;
		rank_to_ratio[mid]=ratio;
		if (perch::_Debug) program.outf<<mid<<'\t'<<total_ti<<'\t'<<total_tv<<'\t'<<ratio<<endl;
	}
	// determin cutoff
	if (rank_to_ratio.size()<3) return 0;
	double x1=rank_to_ratio.begin()->first;
	double y1=rank_to_ratio.begin()->second;
	double xn=rank_to_ratio.rbegin()->first;
	double yn=rank_to_ratio.rbegin()->second;
	// solve y=ax+c
	double a =(yn-y1)/(xn-x1);
	double c =y1-a*x1;
	// cerr<<"y="<<a<<"x+"<<c<<endl;
	// solve Ax+By+C=0
	double A =a;
	double B =-1;
	double C =c;
	// cerr<<A<<"x+"<<B<<"y+"<<C<<"=0"<<endl;
	// find min d
	double max_dist = -std::numeric_limits<double>::max();
	double max_rank = std::numeric_limits<double>::signaling_NaN();
	for (auto &point:rank_to_ratio)
	{
		double m=point.first;
		double n=point.second;
		if (m==x1 || m==xn) continue;
		double d=std::abs(A*m+B*n+C)/sqrt(A*A+B*B);
		// cerr<<i<<" "<<m<<" "<<n<<" "<<d<<endl;
		if (d>max_dist) { max_dist=d; max_rank=m; }
	}
	if (!std::isnan(max_rank))
	{
		lns<<showl<<"# Using TiTv, the suggested VQSLOD cutoff for SNV is "<<max_rank<<flush_logger;
	}
	return 0;
}
