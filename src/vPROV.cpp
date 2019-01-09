#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_mtjobs.hpp>
#include "victor_par.hpp"
#include "vAnnGene_anc.hpp"

using namespace std;

std::mutex												out_mutex;
map<int, map<int, map<string, map<string,double> > > >	out_data;

// multithread jobs
struct JobData {
	string				tx;
	map<string,string>	vars; // vars[hgvs]=index
	JobData(const string& t) { tx=t; }
	void add(const string& hgvs, const string& index) { vars[hgvs]=index; }
};

void run_job(JobData& J)
{
	genepi::gene_info& g=genepi::gene_byTranscript[J.tx].begin()->second;
	if (g.cdsStart == g.cdsEnd) return;
	string RNA = genepi::RNA_seq(g,false);
	int offset=0;
	if (exist_element(ancestral_var_r,g.name))
	{
		for (auto it=ancestral_var_r[g.name].rbegin(); it!=ancestral_var_r[g.name].rend(); ++it)
		{
			int mRNA_loc = g.cdsLen+g.u5Len;	// location in RNA, 1-based
			const int& bgn = it->first.first;	// location in RNA, 1-based
			const int& end = it->first.second;	// location in RNA, 1-based
			const string& rfa = it->second.first;
			const string& ala = it->second.second;
			if (bgn<=g.u5Len) break;
			if		(rfa=="-")
			{
				RNA.insert(bgn,ala);
				if (end<=mRNA_loc) offset+=ala.size();
				//cerr<<g.name<<'\t'<<bgn<<'\t'<<end<<'\t'<<rfa<<'\t'<<ala<<endl;
			}
			else if (ala=="-")
			{
				if (!(bgn<=mRNA_loc && mRNA_loc<=end))
				{
					RNA.erase(bgn-1,rfa.size());
					if (end<mRNA_loc) offset-=rfa.size();
					//cerr<<g.name<<'\t'<<bgn<<'\t'<<end<<'\t'<<rfa<<'\t'<<ala<<endl;
				}
			}
			else // (rfa!="-" && ala!="-") => SNV => bgn=end
			{
				RNA[bgn-1]=ala[0];
				//cerr<<g.name<<'\t'<<bgn<<'\t'<<end<<'\t'<<rfa<<'\t'<<ala<<endl;
			}
		}
	}
	string CDS=RNA;
	CDS.erase(0,g.u5Len);
	CDS.resize(g.cdsLen+offset);
	string pr = genepi::translate_cds(CDS,true);
	if (str_endsw(pr,"*")) pr.pop_back();
	openOutFile_or_exit(fa_out,program.prefix()+"."+J.tx+".fa");
	fa_out<<'>'<<J.tx<<endl<<pr<<endl;
	closefile(fa_out);
	openOutFile_or_exit(va_out,program.prefix()+"."+J.tx+".var");
	for (auto& v:J.vars) va_out<<v.first<<endl;
	closefile(va_out);
	try { exec("provean.sh -q "+program.prefix()+"."+J.tx+".fa -v "+program.prefix()+"."+J.tx+".var > "+program.prefix()+"."+J.tx+".out", false); }
	catch (...) { exit_error("failed to execute provean.sh"); }
	map<string,double> scores; // scores[index]=provean_score
	tfile_format	format;
	format.set_delimiters("\t");
	format.forbid_nf_rpt();
	for (Rows_in_File(in,program.prefix()+"."+J.tx+".out",&format))
	{
		if (in[0]=="No variations entered") continue;
		if (str_startsw(in[0],"[")) continue;
		if (str_has(in[0],"reference AA does not match")) continue;
		if (str_has(in[0],"invalid position")) continue;
		if (in.NumFields()!=2) exit_error("Unknown file format of "+program.prefix()+"."+J.tx+".out");
		if (in[0].empty()) exit_error("file "+program.prefix()+"."+J.tx+".out column 1 is empty");
		if (!exist_element(J.vars,in[0])) exit_error("file "+program.prefix()+"."+J.tx+".out column 1 is unknown");
		double provean_score=std::numeric_limits<double>::signaling_NaN();
		if (!read_val_noNaN(in[1],provean_score)) exit_error("failed to read provean score "+in[1]+" in "+program.prefix()+"."+J.tx+".out");
		scores[J.vars[in[0]]]=provean_score;
	}
	exec("rm -f "+program.prefix()+"."+J.tx+".*", false);
	if (!scores.empty())
	{
		for (auto &i:scores)
		{
			vector<string> in;
			boost::split(in,i.first,boost::is_any_of(","));
			int chr_num = genepi::read_chr_num(in[0]); if (chr_num<1) exit_error("Failed to read "+in[0]+" as a chromosome.");
			int bp=-1; if (!read_val_ge(in[1],bp,1)) exit_error("Failed to read "+in[1]+" as a bp.");
			string& ref = in[2];
			string& alt = in[3];
			out_mutex.lock();
			if (exist_element(out_data[chr_num][bp][ref],alt))	{ if (out_data[chr_num][bp][ref][alt]<i.second) out_data[chr_num][bp][ref][alt]=i.second; }
			else												{ out_data[chr_num][bp][ref][alt]=i.second; }
			out_mutex.unlock();
		}
	}
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	vector<string>	var_in;
	string			anc_in = "provided";

	// handle program options
	program.set_prefix(program.name());
	program.enable_option("--prefix");
	program.enable_option("--nt");
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--anc"))			ReadArg(program.arg(),argi,anc_in);
		else var_in.push_back(program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_anc",anc_in);
	perch::check_arguments();

	// check errors
	if		(anc_in=="provided" && str_startsw(genepi::GDB_name,"refGene")) anc_in = perch::DBpath()+"refGene.anc";
	else if (anc_in=="provided" && str_startsw(genepi::GDB_name,"ensGene")) anc_in.clear();
	if (program.prefix().empty())	exit_error("--prefix cannot be empty");
	
	// read gene database
	lns<<showl<<"Read gene data from "<<genepi::GDB_name<<flush_logger;
	genepi::read_genes();

	// read ancestral variants
	if (!anc_in.empty())
	{
		lns<<showl<<"Read ancestral variants from "<<anc_in<<flush_logger;
		read_ancestral_var(anc_in);
	}
	
	// read variants
	map<string, map<string,string> > all_data; // all_data[tx][hgvs]=index
	for (Rows_in_File(in,var_in,3)) // index tx hgvs
	{
		if (!exist_element(genepi::gene_byTranscript,in[1])) continue;
		genepi::gene_info& g=genepi::gene_byTranscript[in[1]].begin()->second;
		if (g.cdsStart == g.cdsEnd) continue;
		all_data[in[1]][in[2]]=in[0];
	}
	
	// make jobs
	MtJobs<JobData> CompJobs;
	for (auto &x:all_data)
	{
		CompJobs.push_back(JobData(x.first));
		JobData& J = CompJobs.back();
		for (auto& v:x.second) J.add(v.first,v.second);
	}
	if (CompJobs.empty()) exit_error("There is no job created.");

	// run jobs
	lns<<showl<<"Runing "<<CompJobs.size()<<" jobs with "<<program.nt<<" threads."<<flush_logger;
	if (program.nt==1)	for (auto &J:CompJobs) run_job(J);
	else				CompJobs.run(run_job, program.nt);

	// write
	program.outf<<"#CHROM\tPOS\tREF\tALT\tPROVEAN"<<endl;
	for (auto &chr:out_data)
		for (auto &bp:chr.second)
			for (auto &ref:bp.second)
				for (auto &alt:ref.second)
					program.outf<<genepi::convert_chr_num(chr.first)<<'\t'<<bp.first<<'\t'<<ref.first<<'\t'<<alt.first<<'\t'<<alt.second<<endl;
	
	return 0;
}
