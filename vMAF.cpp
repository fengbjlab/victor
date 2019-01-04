/*
 Removed option descriptions:
 --wr S            Output header {_Default_wr}
 --rv D            For #RV counting: variants with MAF<=D are rare variants (0 = no #RV counting) {_Default_rv}
 --wc F            For #RV counting: well-covered regions {_Default_wc}
 --tot-spl I       For #RV counting: the total number of samples. If not set, will calculate from --pop AN. {_Default_tot_spl}
 --ext-af F        For #RV counting: external allele frequency field in INFO. {_Default_ext_af}
 */

#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_program.hpp>
#include "vDEL.hpp"
#include "victor_par.hpp"

using namespace std;
typedef genepi::genotype GTP;
typedef std::tuple<int,int,string,string> ID_t;

vector<double> het_probs; // for HWE test, no need to clear for each use.

// HWE test, returns p-value, modified from http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.c
double my_SNPHWE(const int obs_hets, const int obs_hom1, const int obs_hom2, vector<double>& het_probs)
{
	int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
	int rare_copies = 2 * obs_homr + obs_hets;
	int genotypes   = obs_hets + obs_homc + obs_homr;
	if (genotypes==0 || rare_copies==0) return 1;
	
	het_probs.assign (rare_copies+1, 0);
	
	/* start at midpoint */
	// int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
	// the above is not robust to large sample size (e.g., Het=25768, Hom=22444, AC=70656, AN=111176). So change to below
	int mid = rare_copies * ( 1.0 - rare_copies / (2.0 * genotypes) );

	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1)) mid++;
	
	int curr_hets = mid;
	int curr_homr = (rare_copies - mid) / 2;
	int curr_homc = genotypes - curr_hets - curr_homr;
	
	
	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
	{
		het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		sum += het_probs[curr_hets - 2];
		
		/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
		++curr_homr;
		++curr_homc;
	}
	
	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
	{
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));
		sum += het_probs[curr_hets + 2];
		
		/* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
		--curr_homr;
		--curr_homc;
	}
	
	for (int i = 0; i <= rare_copies; ++i) het_probs[i] /= sum;
	
	double p_hwe = 0.0;
	for (int i = 0; i <= rare_copies; ++i)
	{
		if (het_probs[i] > het_probs[obs_hets]) continue;
		p_hwe += het_probs[i];
	}
	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
	return p_hwe;
}

// right align a variant, return number of bp shifted. Only pure Del or pure Ins can be shifted. DelIns cannot. Assum no Ns in ref1 or alt1.
int right_normalize(const int chr_num, const int pos1, const string& ref1, const string& alt1, int& pos2, string& ref2, string& alt2)
{
	if (ref1==alt1) exit_error("REF and ALT should not be the same.");
	if (ref1.empty() || alt1.empty()) exit_error("REF or ALT cannot be empty.");
	if (ref1.find('N')!=std::string::npos) exit_error("'N' is not allowed in REF. Culprit: "+itos(chr_num)+" "+itos(pos1)+" "+ref1+" "+alt1);
	if (alt1.find('N')!=std::string::npos) exit_error("'N' is not allowed in ALT. Culprit: "+itos(chr_num)+" "+itos(pos1)+" "+ref1+" "+alt1);
	pos2=pos1;
	ref2=ref1;
	alt2=alt1;
	int shifted = 0;
	while (ref2.size()>1 && alt2.size()>1 && ref2.back()==alt2.back()) { ref2.pop_back(); alt2.pop_back(); }
	while (ref2.size()>1 && alt2.size()>1 && ref2.front()==alt2.front()) { ref2.erase(0,1); alt2.erase(0,1); ++pos2; }
	int len=ref2.size();
	if		(alt2.size()==1 && str_startsw(ref2,alt2)) // pure del
	{
		for (;;)
		{
			string next = genepi::DNA_seq(chr_num,pos2+len,1);
			if (next[0]==ref2[1])
			{
				ref2=ref2.substr(1)+next;
				alt2=next;
				++pos2;
				++shifted;
				continue;
			}
			break;
		}
		return shifted;
	}
	else if (ref2.size()==1 && str_startsw(alt2,ref2)) // pure ins
	{
		for (;;)
		{
			string next = genepi::DNA_seq(chr_num,pos2+len,1);
			if (next[0]==alt2[1])
			{
				alt2=alt2.substr(1)+next;
				ref2=next;
				++pos2;
				++shifted;
				continue;
			}
			break;
		}
		return shifted;
	}
	else return 0;
}

int main (int argc, char * const argv[])
{
	perch::h_afID = {	// header of allele frequency IDs
		perch::h_MxAF,		// It's important to include self (MaxAF).
		// "EspAA_A","EspEA_A","KgAmrAA","KgAsnAA","KgAfrAA","KgEurAA","G1K_AF","ESP_AF","ExAC_AF","UK10K_AF", // my annotations. removed to make it clean.
		"EAS_AF","AMR_AF","AFR_AF","EUR_AF","SAS_AF"};	// from KG. AF also in ExAC and KG, but was removed later because it considered all samples.
	// "AF_AFR","AF_AMR","AF_ASN","AF_EUR","AF_MAX" };	// from UK10k, but they are actually G1k frequencies, plus they cannot be correctly vSPLIT, so removed.;
	
	// basic parameters
	perch::prelude(argc,argv);

	// other parameters
	vector<string>	inputs;
	string			wcFile = "minimum.bed";								// region for summary statistics of RV
	string			ext_af;												// external af field
	string			log_fn;												// write variants filtered by MAC to this file
	int				doWhat=0;											// 0=annotate 1=merge_regions 2=grep_well_cov
	int				TotSpl=0;											// total number of samples for sum_RV and FltMAC
	double			sum_RV=0;											// def of rare var (use 0.05 to be robust to small cohorts becaues the prob of no obs from 200 chr is 0.000035)
	double			HW_pvl=0.000001;									// for QC: p-value threshold for HWE test, prv 0.000001
	bool			AddInf=true;										// add a field in INFO, otherwise add a column
	bool			r_norm=true;										// right normalize. Good for annotation, bad for finding a non-redundant variant set.
	bool			OutMis=false;										// write variants even it has a missing MaxAF
	string			wrWhat=perch::h_MxAF;										// header
	string			MisStr="NA";										// missing value
	set<string>		popIDs={"_AFR","_AMR","_EAS","_NFE","_SAS"};		// cal AF by ACpop/ANpop (ExAC). Removed "_ASJ" "_FIN" to avoid founder effect. Removed "_OTH" because it shows a founder effect too.
	int				popMin=200;											// minimum AN, previously 2
	int				FltMAC=0;											// filter variants by minor allele acount (0 = no filter), suggest 5
	
	// handle program options
	string			program_arguments;
	for (int argi=0;argi<argc;++argi) { if (argi) program_arguments.push_back(' '); program_arguments+=string(argv[argi]); }
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1;argi<program.arg().size();++argi)
	{
		if		(program.arg()[argi]=="--merge")	doWhat=1; // rm "argi==1" because it may not be the first argument
		else if	(program.arg()[argi]=="--grep")		doWhat=2; // if users set "alias vMAF='vMAF --ann=$BING/work/data/MaxAF/MaxAF_UGEnoTCGA.gz'"
		else if (str_startsw(program.arg()[argi],"--wc"))			ReadArg(program.arg(),argi,wcFile);
		else if (str_startsw(program.arg()[argi],"--wr"))			ReadArg(program.arg(),argi,wrWhat);
		else if (str_startsw(program.arg()[argi],"--ms"))			ReadArg(program.arg(),argi,MisStr);
		else if (str_startsw(program.arg()[argi],"--pop"))			ReadSet(program.arg(),argi,popIDs);
		else if (str_startsw(program.arg()[argi],"--filt-an"))		ReadArg(program.arg(),argi,popMin);
		else if (str_startsw(program.arg()[argi],"--filt-mac"))		ReadArg(program.arg(),argi,FltMAC);
		else if (str_startsw(program.arg()[argi],"--filt-hwe-pv"))	ReadArg(program.arg(),argi,HW_pvl);
		else if (str_startsw(program.arg()[argi],"--log"))			ReadArg(program.arg(),argi,log_fn);
		else if (str_startsw(program.arg()[argi],"--add-info"))		ReadArg(program.arg(),argi,AddInf);
		else if (str_startsw(program.arg()[argi],"--right-norm"))	ReadArg(program.arg(),argi,r_norm);
		else if (str_startsw(program.arg()[argi],"--tot-spl"))		ReadArg(program.arg(),argi,TotSpl);
		else if (str_startsw(program.arg()[argi],"--rv"))			ReadArg(program.arg(),argi,sum_RV);
		else if (str_startsw(program.arg()[argi],"--ext-af"))		ReadArg(program.arg(),argi,ext_af);
		else if (str_startsw(program.arg()[argi],"--out-ms"))		ReadArg(program.arg(),argi,OutMis);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else add_to_container(inputs,program.arg()[argi]);
	}
	
	// show help
	program.help_text_var("_Default_wc",wcFile);
	program.help_text_var("_Default_wr",wrWhat);
	program.help_text_var("_Default_ms",MisStr);
	program.help_text_var("_Default_pop",str_of_container(popIDs,',',false));
	program.help_text_var("_Default_filt_hwe_pv",ftos(HW_pvl));
	program.help_text_var("_Default_filt_an",itos(popMin));
	program.help_text_var("_Default_filt_mac",itos(FltMAC));
	program.help_text_var("_Default_out_log",log_fn);
	program.help_text_var("_Default_AddInfo",str_YesOrNo(AddInf));
	program.help_text_var("_Default_right_norm",str_YesOrNo(r_norm));
	program.help_text_var("_Default_out_ms",str_YesOrNo(OutMis));
	program.help_text_var("_Default_tot_spl",itos(TotSpl));
	program.help_text_var("_Default_rv",ftos(sum_RV));
	program.help_text_var("_Default_ext_af",ext_af);
	program.push_back_help(GTP::help_text());
	perch::check_arguments();
	
	// check error
	if (FltMAC)
	{
		if (!TotSpl)
			exit_error("--filt-mac must be used together with --tot-spl, because I need that to check whether your cutoff value is OK.");
		if (FltMAC>=TotSpl)
		{
			lns<<showw<<"sample size is too small; --filt-mac="<<FltMAC<<" is ignored. If you still want to use --filt-mac, change it to a smaller threshold."<<flush_logger;
			FltMAC=0;
		}
	}

	// data for #RV counting
	map<int,int>		Num_RV;		// Num_RV[chr_num]=number_of_RV
	genepi::ChrRegions	min_bed;
	if (wcFile.empty()) min_bed.set_whole_genome(true);
	else				min_bed.setup(perch::find_file(wcFile),true,0,true);

	if (doWhat==0)
	{
		boost::iostreams::filtering_ostream logout;
		if (!log_fn.empty()) { if (!openOutFile(logout, log_fn)) exit_error("cannot open log file."); }
		
		field_numbers	FldFrq(false,true);	// field numb for Allele Frequencies
		field_numbers	FldInf(false,true);	// field numb for INFO
		field_numbers	FldChr(false,true);	// field numb for #CHROM
		field_numbers	FldPos(false,true);	// field numb for POS
		field_numbers	FldRef(false,true);	// field numb for REF
		field_numbers	FldAlt(false,true);	// field numb for ALT
		field_numbers	FldRes(false,true);	// field numb for result
		tfile_format	format;
		format.set_delimiters("\t");
		format.set_option(SKIP_NOTES,false);
		bool header_not_read = true;
		map<string,int> pop_add,pop_del;
		for (Rows_in_File(in, inputs, &format))
		{
			if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
			{
				print_container(in.contents(),program.outf,' ',true);
				in.clear_nf();
				continue;
			}
			if (exist_any(perch::h_col1, in.contents()))
			{
				header_not_read=false;
				format.clear_field_nums();
				FldFrq.clear();
				FldInf.clear();
				FldChr.clear();
				FldPos.clear();
				FldRef.clear();
				FldAlt.clear();
				FldRes.clear();
				for (int i=0;i<in.NumFields();++i)
				{
					if (exist(perch::h_afID,in[i]))		FldFrq.push_back(i+1);
					if (in[i]=="INFO")					FldInf.push_back(i+1);
					if (in[i]=="Chr"   && FldChr.no_input())	FldChr.push_back(i+1); if (in[i]=="#CHROM")		{	FldChr.clear(); FldChr.push_back(i+1); }
					if (in[i]=="Start" && FldPos.no_input())	FldPos.push_back(i+1); if (in[i]=="POS")		{	FldPos.clear(); FldPos.push_back(i+1); }
					if (in[i]=="Ref"   && FldRef.no_input())	FldRef.push_back(i+1); if (in[i]=="REF")		{	FldRef.clear(); FldRef.push_back(i+1); }
					if (in[i]=="Alt"   && FldAlt.no_input())	FldAlt.push_back(i+1); if (in[i]=="ALT")		{	FldAlt.clear(); FldAlt.push_back(i+1); }
				}
				if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
				if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
				if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
				if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
				if (AddInf) { if (FldInf.no_input()) { FldRes.push_back(in.NumFields()+1); in.contents().push_back("INFO"); format.set_field_nums(FldRes,"",tfile_format::Expand); } else FldRes=FldInf; }
				else								 { FldRes.push_back(in.NumFields()+1); in.contents().push_back(wrWhat); format.set_field_nums(FldRes,"",tfile_format::Expand); }
				format.set_field_nums(FldFrq,"lines missing the columns for allele frequencies.",tfile_format::Continue);
				format.set_field_nums(FldChr,"lines missing the #CHROM column.",tfile_format::Continue);
				format.set_field_nums(FldPos,"lines missing the POS column.",tfile_format::Continue);
				print_container(in.contents(),program.outf,DLMTR,true);
				continue;
			}
			if (header_not_read) exit_error("Header lines missing.");
			
			// check file
			{
				if (in[FldAlt[0]].find(',')!=std::string::npos) exit_error("The input file has not been split by alternative alleles.");
				if (in[FldRef[0]]=="-" || in[FldAlt[0]]=="-") exit_error("Input is not left-normalized. Culprit: "+in[FldChr[0]]+" "+in[FldPos[0]]+" "+in[FldRef[0]]+" "+in[FldAlt[0]]);
				size_t RefLen = in[FldRef[0]].size();
				size_t AltLen = in[FldAlt[0]].size();
				if (RefLen==AltLen && RefLen>1)
				{
					int sub = 0;
					for (size_t i=0;i<RefLen;++i)
						if (in[FldRef[0]][i]!=in[FldAlt[0]][i]) ++sub;
					if (sub==1) exit_error("Input is not left-normalized. Culprit: "+in[FldChr[0]]+" "+in[FldPos[0]]+" "+in[FldRef[0]]+" "+in[FldAlt[0]]);
				}
			}
			
			// get allele frequencies from AF fields
			vector<double> afs;
			FldFrq.contents_to_doubles(in.contents(),true,afs,false);

			// basic info
			int	chr_num;	if (!genepi::read_chr_num(in[FldChr[0]],chr_num))	exit_error("Failed to read "+in[FldChr[0]]+" as a chromosome.");
			int	bp;			if (!read_val_ge(in[FldPos[0]],bp,1))				exit_error("Failed to read "+in[FldPos[0]]+" as a position in basepairs.");
			bool is_auto = genepi::is_autosomal(chr_num);
			
			double ext_af_val = 0;
			if (!FldInf.no_input())
			{
				// get INFO
				vector<string> INFO;
				if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
				bool info_modified=false;
				for (vector<string>::iterator it = INFO.begin(); it != INFO.end(); it++)
					if (str_startsw(*it,wrWhat+"=")) { INFO.erase(it); info_modified=true; break; }
				if (info_modified)
				{
					in[FldInf[0]] = str_of_container(INFO,';');
					if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
				}
				
				// get allele frequencies from AF in INFO
				for (auto& id:perch::h_afID)
					for (auto& f:INFO)
						if (str_startsw(f,id) && f[id.size()]=='=')
						{
							try { afs.push_back(boost::lexical_cast<double>(f.substr(id.size()+1))); break; }
							catch (...) {  }
						}
				
				// get allele frequencies from AC/AN in INFO
				int ThisSplSz=0;
				int Total_AC=0;
				int Total_AN=0;
				vector<double> temp_afs;
				for (auto &pop:popIDs)
				{
					int AC, AN;
					bool good = (read_val_ge(get_string(INFO,"AC"+pop),AC,0) &&
								 read_val_ge(get_string(INFO,"AN"+pop),AN,popMin));
					if (is_auto) ThisSplSz+=AN/2;
					if (HW_pvl && is_auto && good)
					{
						int Hom;
						if (read_val_ge(get_string(INFO,"Hom"+pop),Hom,0))
						{
							int Het = AC - Hom*2;
							if (Het<0) { good=false; ++pop_del[pop]; }
							else
							{
								double pValue = my_SNPHWE(Het,AN/2-Het-Hom,Hom,het_probs);
								if (pValue<=HW_pvl) { good=false; ++pop_del[pop]; }
							}
						}
					}
					if (good)
					{
						Total_AC+=AC;
						Total_AN+=AN;
					}
					if (good && AC && AN) { temp_afs.push_back((double)AC/AN); ++pop_add[pop]; }
				}
				int MAC = std::min(Total_AC,Total_AN-Total_AC);
				if ( FltMAC && MAC<FltMAC)
				{
					if (!log_fn.empty())
						for (auto &x:temp_afs)
							logout<<in[FldChr[0]]<<DLMTR<<in[FldPos[0]]<<DLMTR<<in[FldRef[0]]<<DLMTR<<in[FldAlt[0]]<<DLMTR<<x<<endl;
					afs.clear(); // This is for 1000 Genomes. It has EAS_AF,AMR_AF,AFR_AF,EUR_AF,SAS_AF read to afs. Only AC and AN to calculate MAC.
				}
				else
				{
					for (auto &x:temp_afs) afs.push_back(x);
				}
				if (ThisSplSz>TotSpl) TotSpl=ThisSplSz;
				
				// external AF
				if (!ext_af.empty())
				{
					for (auto& f:INFO)
						if (str_startsw(f,ext_af) && f[ext_af.size()]=='=')
							if (!read_val(f.substr(ext_af.size()+1),ext_af_val)) ext_af_val=0;
				}
			}
			
			// calculate MaxAF. Take the one closest to 0.5, and choose 1 over 0.
			// af=nan if 1) afs's empty 2) all afs is NaN, which is impossible. Otherwise [0,1].
			double af = -1;
			for (auto& x:afs)
			{
				if ( fabs(x-0.5) < fabs(af-0.5) ) af=x;
				if ( x==1 && af==0) af=x;
			}
			if (af==-1) af=std::numeric_limits<double>::signaling_NaN(); // no data

			// calculate MaxMAF.
			for (auto& x:afs) if (x>0.5) x=1-x;
			double maxmaf = std::numeric_limits<double>::signaling_NaN();
			if (!afs.empty()) maxmaf = *std::max_element(afs.begin(),afs.end());
			
			// calculate Num_RV
			bool is_exome_rv = (!std::isnan(maxmaf) && maxmaf && maxmaf<=sum_RV && ext_af_val<=sum_RV && min_bed.contain(chr_num,bp));
			if (TotSpl && sum_RV && is_exome_rv) Num_RV[chr_num]++;
			
			if (!OutMis && (std::isnan(af)||af==0)) continue;
			string result = ftos(af,MisStr);
			if (AddInf) in[FldRes[0]] += ";" + wrWhat + "=" + result;
			else		in[FldRes[0]]  = result;
			print_container(in.contents(),program.outf,DLMTR,true);
		}
		for (auto & pop:popIDs)
		{
			if (pop_add[pop] || pop_del[pop])
				lns<<showl<<pop_add[pop]<<" variants read and "<<pop_del[pop]<<" removed by HWD from "<<pop<<flush_logger;
		}
		closefile(logout);
		if (TotSpl && sum_RV)
		{
			lns << showl << "Total sample size is " << TotSpl << flush_logger;
			for (auto &x:Num_RV)
				lns << showl << "Chromosome " << x.first << " has " << x.second / (double)TotSpl << " rare variants per sample" << flush_logger;
		}
	}
	else if (doWhat==1)
	{
		map<ID_t,double> annotation; // index => MaxAF
		int warning1 = elog.get_token("lines omitted because allele frequency is 0.");
		for (auto &f:inputs)
		{
			lns<<showl<<"Read " << f << flush_logger;
			for (Rows_in_File(in,f,5))
			{
				if (exist_any(perch::h_col1, in.contents())) continue;
				int		ch; if (!genepi::read_chr_num(in[0],ch))	exit_error("Failed to read "+in[0]+" as a chromosome.");
				int		bp; if (!read_val_ge(in[1],bp,1))			exit_error("Failed to read "+in[1]+" as a position in basepairs.");
				double	af; if (!read_val_ge_le(in[4],af,0.0,1.0))	exit_error("Failed to read "+in[4]+" as an allele frequency.");
				if (af==0) { elog.add(warning1); continue; }
				string& ref=in[2];
				string& alt=in[3];
				if (ref.find('N')!=std::string::npos) continue;
				if (alt.find('N')!=std::string::npos) continue;
				{
					ID_t ID = std::make_tuple (ch, bp, ref, alt);
					map<ID_t,double>::iterator it = annotation.find(ID);
					if (it == annotation.end())
					{
						annotation[ID]=af;
					}
					else
					{
						double& pr = it->second; // previous AF (the closest to 0.5)
						if (std::abs(af-0.5)<=std::abs(pr-0.5)) pr=af;
					}
				}
				int l2;
				string rf2, al2;
				if (r_norm && (ref.size()!=1 || alt.size()!=1) && right_normalize(ch, bp, ref, alt, l2, rf2, al2))
				{
					ID_t ID = std::make_tuple (ch, l2, rf2, al2);
					map<ID_t,double>::iterator it = annotation.find(ID);
					if (it == annotation.end())
					{
						annotation[ID]=af;
					}
					else
					{
						double& pr = it->second; // previous AF (the closest to 0.5)
						if (std::abs(af-0.5)<=std::abs(pr-0.5)) pr=af;
					}
				}
			}
		}
		program.outf << "##created by "<<program_arguments<<endl;
		program.outf << "#CHROM\tPOS\tREF\tALT\t" << wrWhat << endl;
		for (auto &i:annotation)
		{
			program.outf<< genepi::convert_chr_num(std::get<0>(i.first)) << DLMTR
						<< std::get<1>(i.first) << DLMTR
						<< std::get<2>(i.first) << DLMTR
						<< std::get<3>(i.first) << DLMTR
						<< ftos(i.second,2) << endl;
		}
	}
	else if (doWhat==2)
	{
		tfile_format		format;
		format.set_delimiters("\t");
		format.set_option(SKIP_NOTES,false);
		format.forbid_nf_rpt();
		for (Rows_in_File(in, inputs, &format))
		{
			if (in[0].size()>0 && in[0][0]=='#')
			{
				print_container(in.contents(),program.outf,DLMTR,true);
				continue;
			}
			try {
				int	chr_num;	if (!genepi::read_chr_num(in[0],chr_num))	exit_error("Failed to read "+in[0]+" as a chromosome.");
				int	bp;			if (!read_val_ge(in[1],bp,1))				exit_error("Failed to read "+in[1]+" as a position in basepairs.");
				if (min_bed.contain(chr_num,bp)) print_container(in.contents(),program.outf,DLMTR,true);
			} catch (...) { print_container(in.contents(),program.outf,DLMTR,true); }
		}
	}
	else exit_error("Unknown run mode.");
	return 0;
}
