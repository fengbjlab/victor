#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include "victor_par.hpp"

using namespace std;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	vector<string>	inputs;
	set<double>		cutoff = {10,20,25,30,50,75,100};
	bool			omit_unk_chr=true;
	double			min_dp = 34; // http://www.ncbi.nlm.nih.gov/pubmed/25038816
	double			pc_cov = 0.8;
	string			f_type;
	string			t_file;
	string			l_file;
	string			mdFile;
	string			o_sumg;
	int				up_reg = 30;
	int				dn_reg = 0;
	int				in_reg = 12;
	string			desert; 	// eg 5k 1m 500bp. Usage: vBED --prefix=try --format=bed --split=500k maximum.bed
	int				desert_bp=0;
	int				pad_bp = 0;	// padding, for bed file only

	// handle program options
	genepi::canonic=true; // use canonical transcript only
	program.enable_option("--prefix");
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(str_startsw(program.arg()[argi],"--format"))	ReadArg(program.arg(),argi,f_type);
		else if (str_startsw(program.arg()[argi],"--pc"))		ReadArg(program.arg(),argi,pc_cov);
		else if	(str_startsw(program.arg()[argi],"--cutoff"))	ReadArg(program.arg(),argi,min_dp);
		else if (str_startsw(program.arg()[argi],"--dp"))		ReadSet(program.arg(),argi,cutoff);
		else if (str_startsw(program.arg()[argi],"--up"))		ReadArg(program.arg(),argi,up_reg);
		else if (str_startsw(program.arg()[argi],"--dn"))		ReadArg(program.arg(),argi,dn_reg);
		else if (str_startsw(program.arg()[argi],"--intron"))	ReadArg(program.arg(),argi,in_reg);
		else if (str_startsw(program.arg()[argi],"--on-target"))ReadArg(program.arg(),argi,t_file);
		else if (str_startsw(program.arg()[argi],"--per-gene"))	ReadArg(program.arg(),argi,o_sumg);
		else if (str_startsw(program.arg()[argi],"--log"))		ReadArg(program.arg(),argi,l_file);
		else if (str_startsw(program.arg()[argi],"--mean"))		ReadArg(program.arg(),argi,mdFile);
		else if (str_startsw(program.arg()[argi],"--split"))	ReadArg(program.arg(),argi,desert);
		else if (str_startsw(program.arg()[argi],"--padding"))	ReadArg(program.arg(),argi,pad_bp);
		else if (str_startsw(program.arg()[argi],"-"))			exit_error("unknown option "+program.arg()[argi]);
		else inputs.push_back(program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_DP",str_of_container(cutoff,',',false));
	program.help_text_var("_Default_PC",ftos(pc_cov));
	program.help_text_var("_Default_format",f_type);
	program.help_text_var("_Default_target",t_file);
	program.help_text_var("_Default_summary",o_sumg);
	program.help_text_var("_Default_log",l_file);
	program.help_text_var("_Default_mean",mdFile);
	program.help_text_var("_Default_cutoff",ftos(min_dp));
	program.help_text_var("_Default_up",itos(up_reg));
	program.help_text_var("_Default_dn",itos(dn_reg));
	program.help_text_var("_Default_intron",itos(in_reg));
	program.help_text_var("_Default_split",desert);
	program.help_text_var("_Default_padding",itos(pad_bp));
	perch::check_arguments();

	// check errors
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
	if (inputs.empty()) inputs.push_back(label_stdin());
	else for (auto &x:inputs) x=perch::find_file(x);
	
	// setup individual names
	vector<string> SeqID=inputs;
	if (inputs.size()>1)
	{
		vector< vector<string> > words;
		size_t size=0;
		bool size_diff=false;
		for (auto &infile:inputs)
		{
			if (infile.empty()) exit_error("empty string is not allow within a list of filenames");
			vector<string> line;
			boost::split(line,infile,boost::is_any_of("./_"));
			words.push_back(line);
			if (size==0) size=line.size();
			else if (size!=line.size()) { size_diff=true; break; }
		}
		if (!size_diff)
		{
			size_t first = 0;
			for (size_t i=0;i<size;++i)
			{
				bool str_diff=false;
				string str="dummy_string_that_should_never_happen";
				for (auto &person:words)
				{
					if (str=="dummy_string_that_should_never_happen") str=person[i];
					else if (str!=person[i]) { str_diff=true; break; }
				}
				if (str_diff) break;
				first += words[0][i].size();
				if (i) ++first;
			}
			size_t last = 0;
			for (int i=size-1;i>=0;--i)
			{
				bool str_diff=false;
				string str="dummy_string_that_should_never_happen";
				for (auto &person:words)
				{
					if (str=="dummy_string_that_should_never_happen") str=person[i];
					else if (str!=person[i]) { str_diff=true; break; }
				}
				if (str_diff) break;
				last += words[0][i].size();
				if (i!=(int)(size-1)) ++last;
			}
			for (auto &id:SeqID)
			{
				if (last) id.resize(id.size()-last-1);
				if (first) id=id.substr(first+1);
			}
		}
	}
	
	// data
	vector< vector<genepi::ChrRegions> > cr_vec;// cr_vec[level][person]=covered_by_level
	vector<double> levels;						// levels[level]
	size_t threshold_level=0;					// which level is the threshold for good coverage
	if (f_type=="gatk"||f_type=="samtools"||f_type=="multi")
	{
		cutoff.insert(min_dp);
		for (auto &c:cutoff) levels.push_back(c);
		for (threshold_level=0;levels[threshold_level]!=min_dp;++threshold_level) ;
		cr_vec.resize(levels.size());
	}
	else if (f_type=="bed")
	{
		levels.assign(1,-1);
		threshold_level=0;
		cr_vec.resize(1);
	}
	
	// setup regions
	vector<string>				g_id;			// g_id[gene]=gene symbols
	vector<genepi::ChrRegions>	g_cr;			// g_cr[gene]=gene chr regions
	vector<int>					g_sz;			// g_sz[gnee]=gene size
	vector< vector<int>	>		g_cv;			// g_cv[level][gene]=coveraged bp
	if (!t_file.empty())
	{
		g_id.push_back("TargetedRegion");
		g_cr.push_back(genepi::ChrRegions());
		genepi::ChrRegions& cr=g_cr.back();
		cr.setup(t_file,true,0,omit_unk_chr);
		g_sz.push_back(cr.total_bp());
	}
	if (!o_sumg.empty())
	{
		genepi::read_genes();
		for (auto &gene:genepi::gene_bySymbol)
		{
			g_id.push_back(gene.second.begin()->second.name2);
			g_cr.push_back(genepi::ChrRegions());
			genepi::ChrRegions &cr = g_cr.back();
			for (auto &tx:gene.second)
			{
				genepi::gene_info& g = tx.second;
				if (g.strand=='+')
					for (int i=0; i<g.exonCount; ++i)
					{
						if (i==0 && i==g.exonCount-1)	cr.add(g.chrNumPlink,std::max(g.exonStarts[i]-up_reg,0)+1,g.exonEnds[i]+dn_reg);
						else if	(i==0)					cr.add(g.chrNumPlink,std::max(g.exonStarts[i]-up_reg,0)+1,g.exonEnds[i]+in_reg);
						else if (i==g.exonCount-1)		cr.add(g.chrNumPlink,std::max(g.exonStarts[i]-in_reg,0)+1,g.exonEnds[i]+dn_reg);
						else							cr.add(g.chrNumPlink,std::max(g.exonStarts[i]-in_reg,0)+1,g.exonEnds[i]+in_reg);
					}
				else
					for (int i=g.exonCount-1; i>=0; --i)
					{
						if (i==0 && i==g.exonCount-1)	cr.add(g.chrNumPlink,std::max(g.exonStarts[i]-dn_reg,0)+1,g.exonEnds[i]+up_reg);
						else if	(i==0)					cr.add(g.chrNumPlink,std::max(g.exonStarts[i]-dn_reg,0)+1,g.exonEnds[i]+in_reg);
						else if (i==g.exonCount-1)		cr.add(g.chrNumPlink,std::max(g.exonStarts[i]-in_reg,0)+1,g.exonEnds[i]+up_reg);
						else							cr.add(g.chrNumPlink,std::max(g.exonStarts[i]-in_reg,0)+1,g.exonEnds[i]+in_reg);
					}
			}
			cr.finishing();
			g_sz.push_back(cr.total_bp());
		}
	}
	g_cv.assign(levels.size(),vector<int>(g_id.size(),0));
	
	// read inputs
	int num_inputs = inputs.size();
	map<int, map<int,double> > mean_depth;
	bool log_mean_depth = !mdFile.empty();
	if (f_type=="gatk"||f_type=="samtools")
	{
		bool is_gatk = (f_type=="gatk");
		for (size_t infile=0; infile<inputs.size(); ++infile)
		{
			lns<<showl<<"Read "<<inputs[infile]<<flush_logger;
			size_t ind = cr_vec[0].size();
			for (auto &i:cr_vec) i.push_back(genepi::ChrRegions());
			for (Rows_in_File(in,inputs[infile],2))
			{
				if (in[0]=="Locus") continue;
				string chr_str, pos_str, dp_str;
				if (is_gatk)
				{
					if (in.NumFields()<2) exit_error("Insufficient number of columns in "+inputs[infile]);
					chr_str = substr_before_find(in[0],":");
					pos_str = substr_after_find(in[0],":");
					dp_str = in[1];
				}
				else
				{
					if (in.NumFields()<3) exit_error("Insufficient number of columns in "+inputs[infile]);
					chr_str = in[0];
					pos_str = in[1];
					dp_str = in[2];
				}
				int	chr_num;	if (!genepi::read_chr_num(chr_str,chr_num))	{ if (omit_unk_chr) continue; else exit_error("Failed to read "+chr_str+" as a chromosome."); }
				int	bp;			if (!read_val_ge(pos_str,bp,1))		exit_error("Failed to read "+pos_str+" as a position in basepairs.");
				double dp;		if (!read_val_ge(dp_str,dp,0.0))	exit_error("Failed to read "+in[1]+" as read depth.");
				for (size_t lvl=0;lvl<levels.size();++lvl)
					if (dp>=levels[lvl]) cr_vec[lvl][ind].add(chr_num,bp,bp);
				if (log_mean_depth) mean_depth[chr_num][bp]+=dp;
			}
			for (auto &i:cr_vec) i.back().finishing();
		}
	}
	else if (f_type=="multi")
	{
		lns<<showl<<"Read inputs"<<flush_logger;
		for (Rows_in_File(in,inputs,3))
		{
			if (in.RowNumber()==0)
			{
				num_inputs = in.NumFields()-2;
				SeqID.clear(); for (int i=0;i<num_inputs;++i) SeqID.push_back("Spl"+itos(i+1));
				for (auto &i:cr_vec) i.assign(num_inputs,genepi::ChrRegions());
			}
			else
			{
				if (num_inputs != in.NumFields()-2) exit_error("unequal number of columns in inputs");
			}
			if (in[0]=="Locus" || in[0]=="Chr") continue;
			string& chr_str = in[0];
			string& pos_str = in[1];
			int	chr_num;	if (!genepi::read_chr_num(chr_str,chr_num))	{ if (omit_unk_chr) continue; else exit_error("Failed to read "+chr_str+" as a chromosome."); }
			int	bp;			if (!read_val_ge(pos_str,bp,1))		exit_error("Failed to read "+pos_str+" as a position in basepairs.");
			for (int ind=2; ind<in.NumFields(); ++ind)
			{
				string& dp_str = in[ind];
				double dp;	if (!read_val_ge(dp_str,dp,0.0))	exit_error("Failed to read "+in[ind]+" as read depth.");
				for (size_t lvl=0;lvl<levels.size();++lvl)
					if (dp>=levels[lvl]) cr_vec[lvl][ind-2].add(chr_num,bp,bp);
				if (log_mean_depth) mean_depth[chr_num][bp]+=dp;
			}
		}
		for (auto &lvl:cr_vec) for (auto &cr:lvl) cr.finishing();
	}
	else if (f_type=="bed")
	{
		for (auto &infile:inputs)
		{
			cr_vec[threshold_level].push_back(genepi::ChrRegions());
			genepi::ChrRegions &cr = cr_vec[threshold_level].back();
			cr.setup(infile,true,pad_bp,omit_unk_chr);
			lns<<showl<<cr.total_bp()<<" basepairs covered by "<<infile<<flush_logger;
		}
		if (num_inputs==1 && o_sumg.empty())
		{
			if (desert_bp)	cr_vec[threshold_level][0].write(desert_bp,program.prefix());
			else			cr_vec[threshold_level][0].write(program.outf,true);
			return 0;
		}
	}
	else exit_error("--format not set");
	
	// mean depth
	if (log_mean_depth)
	{
		openOutFile_or_exit(md_out,mdFile);
		for (auto &chr:mean_depth)
			for (auto &bp:chr.second)
				md_out<<chr.first<<'\t'<<bp.first<<'\t'<<bp.second/num_inputs<<endl;
		closefile(md_out);
	}
	
	// calculate
	boost::iostreams::filtering_ostream logout;
	if (!l_file.empty() && !openOutFile(logout,l_file)) exit_error("cannot write to "+l_file);
	vector<genepi::ChrRegions> mrged_regions(levels.size());
	for (size_t lvl=0; lvl<levels.size(); ++lvl)
	{
		bool to_log = (!l_file.empty() && lvl==threshold_level);
		lns<<showl<<"For "<<levels[lvl]<<"x"<<flush_logger;
		lns<<showl<<"   Calculate the union of all files, ";
		genepi::ChrRegions union_regions;
		for (auto &cr:cr_vec[lvl])
		{
			map< int, map<int,int> > data = cr.data();
			for (auto &chr:data)
				for (auto &seg:chr.second)
					union_regions.add(chr.first,seg.first,seg.second);
		}
		union_regions.finishing();
		lns<<"which has "<<union_regions.total_bp()<<" basepairs in total"<<flush_logger;
		
		int max_mis = std::ceil( (1-pc_cov) * num_inputs + std::numeric_limits<double>::min());
		if (pc_cov==0) // union_regions.total_bp()==0) // removed max_mis==num_inputs
		{
			mrged_regions[lvl]=union_regions;
			// lns<<showl<<"  "<<union_regions.total_bp()<<" basepairs covered after merging"<<flush_logger;
		}
		else
		{
			lns<<showl<<"Calculate the intersection of all files"<<flush_logger;
			progress_timer pt;
			if (!program.quiet) { pt.prefix("#   Runtime "); pt.start(union_regions.total_bp()); }
			map< int, map<int,int> > data = union_regions.data(); // consider to use reference
			for (auto &chr:data)
				for (auto &seg:chr.second)
					for (int bp=seg.first; bp<=seg.second; ++bp)
					{
						int num_mis=0;
						for (auto &cr:cr_vec[lvl])
							if (!cr.contain(chr.first,bp)) { ++num_mis; if (num_mis>=max_mis) break; }
						if (num_mis<max_mis)
						{
							mrged_regions[lvl].add(chr.first,bp,bp);
							if (to_log) logout << chr.first << '\t' << bp << '\t' << num_inputs-num_mis << endl;
							for (size_t g=0;g<g_cr.size();++g)
								if (g_cr[g].contain(chr.first,bp)) ++g_cv[lvl][g];
						}
						if (!program.quiet) ++pt;
					}
			if (!program.quiet) pt.finish();
			mrged_regions[lvl].finishing();
			lns<<showl<<"  Got "<<mrged_regions[lvl].total_bp()<<" basepairs in total"<<flush_logger;
		}
	}
	if (!l_file.empty()) closefile(logout);
	lns<<showl<<"Write the regions that "<<pc_cov*100<<"% of samples are above "<<min_dp<<"x, ";
	lns<<"which has "<<mrged_regions[threshold_level].total_bp()<<" basepairs"<<flush_logger;
	if (desert_bp)	mrged_regions[threshold_level].write(desert_bp,program.prefix());
	else			mrged_regions[threshold_level].write(program.outf,true);
	if (!o_sumg.empty())
	{
		openOutFile_or_exit(outf,o_sumg);
		outf<<"RegionOrGene"; for (size_t lvl=0; lvl<levels.size(); ++lvl) outf<<'\t'<<levels[lvl]<<'x'; outf<<endl;
		for (size_t g=0;g<g_id.size();++g)
		{
			bool covered = false;
			for (size_t lvl=0; lvl<levels.size(); ++lvl) if (g_cv[lvl][g]) { covered=true; break; }
			if (covered)
			{
				outf << g_id[g];
				for (size_t lvl=0; lvl<levels.size(); ++lvl) outf<<'\t'<<ftos_MaxWidth((double)g_cv[lvl][g]/g_sz[g]);
				outf << endl;
				// lns<<showl<<g_id[g]<<'\t'<<g_sz[g]<<'\t'<<g_cv[threshold_level][g]<<flush_logger;
			}
		}
		closefile(outf);
	}
	
	return 0;
}
