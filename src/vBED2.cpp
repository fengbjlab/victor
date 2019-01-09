#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_mtjobs.hpp>
#include "victor_par.hpp"

using namespace std;

// shared data
std::mutex 				lns_mutex;
std::mutex 				cr_vec_mutex;
vector< vector<int> > 	cr_vec; // cr_vec[chr][bp]=num_samples
vector<int>				max_bp;
bool					omit_unk_chr = true;
double					min_dp = 34; // http://www.ncbi.nlm.nih.gov/pubmed/25038816

// multithread jobs
struct JobData {
	char 	filetype;
	string 	filename;
	JobData(char ft, string fn) { filetype=ft; filename=fn; }
};

void run_job(JobData& J)
{
	lns_mutex.lock();
	lns<<showl<<"Reading "<<J.filename<<flush_logger;
	lns_mutex.unlock();
	if (J.filetype=='g')
	{
		for (Rows_in_File(in,J.filename,2))
		{
			if (in[0]=="Locus") continue;
			string chr_str, pos_str, dp_str;
			if (in.NumFields()<2) exit_error("Insufficient number of columns in "+J.filename);
			chr_str = substr_before_find(in[0],":");
			pos_str = substr_after_find(in[0],":");
			dp_str = in[1];
			int	chr_num;	if (!genepi::read_chr_num(chr_str,chr_num))	{ if (omit_unk_chr) continue; else exit_error("Failed to read "+chr_str+" as a chromosome."); }
			int	bp;			if (!read_val_ge(pos_str,bp,1))		exit_error("Failed to read "+pos_str+" as a 1-based basepair position.");
			double dp;		if (!read_val_ge(dp_str,dp,0.0))	exit_error("Failed to read "+dp_str+" as a read depth number.");
			if (dp<min_dp) continue;
			if (chr_num>=(int)max_bp.size()) exit_error("wrong chr "+itos(chr_num));
			if (bp>max_bp[chr_num]) exit_error("wrong bp "+pos_str+" for chr "+chr_str+" in "+J.filename+". should be <= "+itos(max_bp[chr_num]));
			cr_vec_mutex.lock();
			++cr_vec[chr_num][bp];
			cr_vec_mutex.unlock();
		}
	}
	else if (J.filetype=='s')
	{
		for (Rows_in_File(in,J.filename,3))
		{
			if (in[0]=="Locus"||in[0]=="#CHROM"||in[0]=="Chr") continue;
			string chr_str, pos_str, dp_str;
			if (in.NumFields()<3) exit_error("Insufficient number of columns in "+J.filename);
			chr_str = in[0];
			pos_str = in[1];
			dp_str = in[2];
			int	chr_num;	if (!genepi::read_chr_num(chr_str,chr_num))	{ if (omit_unk_chr) continue; else exit_error("Failed to read "+chr_str+" as a chromosome."); }
			int	bp;			if (!read_val_ge(pos_str,bp,1))		exit_error("Failed to read "+pos_str+" as a 1-based basepair position.");
			double dp;		if (!read_val_ge(dp_str,dp,0.0))	exit_error("Failed to read "+dp_str+" as a read depth number.");
			if (dp<min_dp) continue;
			if (chr_num>=(int)max_bp.size()) exit_error("wrong chr "+itos(chr_num));
			if (bp>max_bp[chr_num]) exit_error("wrong bp "+pos_str+" for chr "+chr_str+" in "+J.filename+". should be <= "+itos(max_bp[chr_num]));
			cr_vec_mutex.lock();
			++cr_vec[chr_num][bp];
			cr_vec_mutex.unlock();
		}
	}
	else if (J.filetype=='v')
	{
		for (Rows_in_File(in,J.filename,3))
		{
			if (in[0]=="Locus"||in[0]=="#CHROM"||in[0]=="Chr") continue;
			string chr_str, pos_str, spl_str;
			if (in.NumFields()<3) exit_error("Insufficient number of columns in "+J.filename);
			chr_str = in[0];
			pos_str = in[1];
			spl_str = in[2];
			int	chr_num;	if (!genepi::read_chr_num(chr_str,chr_num))	{ if (omit_unk_chr) continue; else exit_error("Failed to read "+chr_str+" as a chromosome."); }
			int	bp;			if (!read_val_ge(pos_str,bp,1))	exit_error("Failed to read "+pos_str+" as a 1-based basepair position.");
			int sp;			if (!read_val_ge(spl_str,sp,0))	exit_error("Failed to read "+spl_str+" as an integer number.");
			if (chr_num>=(int)max_bp.size()) exit_error("wrong chr "+itos(chr_num));
			if (bp>max_bp[chr_num]) exit_error("wrong bp "+pos_str+" for chr "+chr_str+" in "+J.filename+". should be <= "+itos(max_bp[chr_num]));
			cr_vec_mutex.lock();
			cr_vec[chr_num][bp]+=sp;
			cr_vec_mutex.unlock();
		}
	}
	else exit_error("wrong file type code");
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	vector<string>	inputs;
	int				MinSpl = 0;
	double			pc_cov = 0.8;
	string			f_type;
	string			l_file;
	string			desert; 	// eg 5k 1m 500bp. Usage: vBED --prefix=try --format=bed --split=500k maximum.bed
	int				desert_bp=0;
	string			inlist;
	int				TotSpl=0;
	
	// handle program options
	program.enable_option("--prefix");
	program.enable_option("--nt");
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--filelist"))	ReadArg(program.arg(),argi,inlist);
		else if	(str_startsw(program.arg()[argi],"--min-spl"))	ReadArg(program.arg(),argi,MinSpl);
		else if	(str_startsw(program.arg()[argi],"--format"))	ReadArg(program.arg(),argi,f_type);
		else if	(str_startsw(program.arg()[argi],"--pc"))		ReadArg(program.arg(),argi,pc_cov);
		else if	(str_startsw(program.arg()[argi],"--dp"))		ReadArg(program.arg(),argi,min_dp);
		else if (str_startsw(program.arg()[argi],"--log"))		ReadArg(program.arg(),argi,l_file);
		else if (str_startsw(program.arg()[argi],"--split"))	ReadArg(program.arg(),argi,desert);
		else if (str_startsw(program.arg()[argi],"--tot-spl"))	ReadArg(program.arg(),argi,TotSpl);
		else if (str_startsw(program.arg()[argi],"-"))			exit_error("unknown option "+program.arg()[argi]);
		else inputs.push_back(program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_format",f_type);
	program.help_text_var("_Default_PC",ftos(pc_cov));
	program.help_text_var("_Default_dp",ftos(min_dp));
	program.help_text_var("_Default_log",l_file);
	program.help_text_var("_Default_TotSpl",itos(TotSpl));
	program.help_text_var("_Default_split",desert);
	program.help_text_var("_Default_min_spl",itos(MinSpl));	
	perch::check_arguments();

	// check errors
	if (f_type!="vBED" && f_type!="gatk" && f_type!="samtools") exit_error("--format not set or wrong. Allowed values are gatk / samtools / vBED.");
	bool MedRes = (f_type=="vBED");
	if (MedRes && MinSpl==0) exit_error("please set --min-spl");
	if (!inlist.empty()) for (Rows_in_File(in,inlist,1)) inputs.push_back(in[0]);
	if (inputs.empty()) inputs.push_back(label_stdin());
	else for (auto &x:inputs) x=perch::find_file(x);
	if (!desert.empty())
	{
		double v;
		boost::to_lower(desert);
		if 		(str_endsw(desert,"bp")) { desert.pop_back(); desert.pop_back(); if (!read_val_noNaN(desert,v)) exit_error("failed to parse --split"); desert_bp=v; }
		else if (str_endsw(desert,"kb")) { desert.pop_back(); desert.pop_back(); if (!read_val_noNaN(desert,v)) exit_error("failed to parse --split"); desert_bp=v*1000; }
		else if (str_endsw(desert,"mb")) { desert.pop_back(); desert.pop_back(); if (!read_val_noNaN(desert,v)) exit_error("failed to parse --split"); desert_bp=v*1000000; }
		else if (str_endsw(desert,"k"))	 { desert.pop_back(); 					 if (!read_val_noNaN(desert,v)) exit_error("failed to parse --split"); desert_bp=v*1000; }
		else if (str_endsw(desert,"m"))	 { desert.pop_back(); 					 if (!read_val_noNaN(desert,v)) exit_error("failed to parse --split"); desert_bp=v*1000000; }
		else exit_error("the --split argument must ends with bp/kb/mb/k/m");
		if (program.prefix().empty()) exit_error("--split should be used together with --prefix");
	}

	// prepare data
	cr_vec.push_back(vector<int>()); 	for (int chr=1;chr<=genepi::MaxChrNum();++chr) cr_vec.push_back(vector<int>(genepi::chrlen_bp(chr)+1,0));
	max_bp.push_back(0);				for (int chr=1;chr<=genepi::MaxChrNum();++chr) max_bp.push_back(genepi::chrlen_bp(chr));
	char ft = f_type[0];
	
	// read inputs
	MtJobs<JobData> CompJobs;
	for (size_t i=0; i<inputs.size(); ++i)
		CompJobs.push_back(JobData(ft,inputs[i]));
	if (CompJobs.empty()) exit_error("There is no job created.");
	lns<<showl<<"Runing "<<CompJobs.size()<<" jobs with "<<program.nt<<" threads."<<flush_logger;
	if (program.nt==1)	for (auto &J:CompJobs) run_job(J);
	else				CompJobs.run(run_job, program.nt);
	
	// calculate
	lns<<showl<<"Merging data"<<flush_logger;
	int min_spl = ( MinSpl!=0 ? MinSpl : std::ceil(inputs.size()*pc_cov) );
	if (min_spl==0) min_spl=1;
	boost::iostreams::filtering_ostream logout;
	if (!l_file.empty() && !openOutFile(logout,l_file)) exit_error("cannot write to "+l_file);
	genepi::ChrRegions mrged_regions;
	for (int chr_num=1;chr_num<=genepi::MaxChrNum();++chr_num)
	{
		for (int bp=1;bp<=genepi::chrlen_bp(chr_num);++bp)
		{
			if (cr_vec[chr_num][bp]>=min_spl)
			{
				mrged_regions.add(chr_num,bp,bp);
				if (!l_file.empty())
				{
					if (TotSpl) logout << genepi::convert_chr_num(chr_num) << '\t' << bp << '\t' << cr_vec[chr_num][bp]/(double)TotSpl << endl;
					else		logout << genepi::convert_chr_num(chr_num) << '\t' << bp << '\t' << cr_vec[chr_num][bp] << endl;
				}
			}
		}
	}
	mrged_regions.finishing();
	if (!l_file.empty()) closefile(logout);
	
	// output
	lns<<showl<<"Writing regions with "<<min_dp<<"x coverage in "<<pc_cov*100<<"% of samples. ";
	lns<<"Total bp is "<<mrged_regions.total_bp()<<flush_logger;
	if (desert_bp)	mrged_regions.write(desert_bp,program.prefix());
	else			mrged_regions.write(program.outf,true);

	return 0;
}
