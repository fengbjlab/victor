/*
 --out-lof F        Output loss-of-function genes to F (col: Symbol SeqID RemainCopy) {_Default_del_log}
*/
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_regress.hpp>
#include "victor_par.hpp"

using namespace std;
typedef genepi::genotype GTP;

double truncated(double pen)
{
	if (pen>=1) return 0.9999999;
	if (pen<=0) return 0.0000001;
	return pen;
}

// These data must be set before any job and be constant afterward.
genepi::gtp_par 			gpar;
set<string>					keepSV = { "DEL", "DUP", "INS", "INV", };
double						filSAF = 0.01;			// skip genes if observed non-2 > filSAF, 0=noQC
double						MisPVl = 0.000001;		// skip genes if missing % between cs and ct has p-val < MisPVl, 0=noQC
vector<double>				CovORs;	// covariates OR/unit: CovORs[var_coor]
vector< vector<double> >	CovVar;	// covariates for REG: CovVar[spl_coor][var_coor]
vector<double>				DepVar;	// Dependent variable: DepVar[spl_coor]
vector<double>				strata;	// strata def by cov : strata[spl_coor]=RRc
vector<string>				record;	// sequenced recordID: record[spl_coor]
vector< vector<int> >		Tx_DEL; // No. of copies of a: Tx_DEL[MyTxID][spl_coor] = 0,1,2.. (diploidy) or 0,2,4.. (haploidy). -1=missing, -2=NA
vector< vector<int> >		Tx_DUP; // No. of copies of a: Tx_DUP[MyTxID][spl_coor] = 0,1,2.. (diploidy) or 0,2,4.. (haploidy). -1=missing, -2=NA
vector< genepi::gene_info >	TxInfo;	// gene info   of tx : TxInfo[MyTxID] = gene_info
boost::iostreams::filtering_ostream logfile;
string								loginfo;
bool								logopen=false;
bool			pdupid = false;		// partial dup is damaging

// genotype is 0/1/2 for number of copies of ALT. when ploidy=1, genotype is 0/2 (automatically 1=>2). delta is -1 for DEL and 1 for DUP.
void populate_genotype_long(vector< vector<int> >& db, int delta, const string& genotype, const string& GQ, int i, int sex, int chr_num, int bp, size_t g)
{
	int ploidy = genepi::expected_ploidy(chr_num,bp,sex);
	if (ploidy==0) { db[g][i]=-2; return; } // NA
	if (db[g][i]<0) return; // if already NA/missing, keep it
	double gq;
	if (!read_val(GQ,gq)) return;
	if (gq<GTP::GQ_cut) { db[g][i]=-1; return; }
	int gt;
	if (!read_val_ge_le(genotype,gt,0,2)) return;
	if (gt==0) return;
	if (delta==0) return;
	if (ploidy==1 && gt==1) gt=2;
	
	db[g][i] = 2 + gt * delta;
	if (logopen && db[g][i]!=2) logfile<<loginfo<<record[i]<<'\t'<<db[g][i]<<endl;
	if (db[g][i]<0) db[g][i]=0; // hope this won't happen
}

// genotype is 0- for number of copies of ALT. delta is -1 for DEL and 1 for DUP.
void populate_genotype_clamms(vector< vector<int> >& db, int delta, const string& genotype, const string& GQ, int i, int sex, int chr_num, int bp, size_t g)
{
	int ploidy = genepi::expected_ploidy(chr_num,bp,sex);
	if (ploidy==0) { db[g][i]=-2; return; } // NA
	if (db[g][i]<0) return; // if already NA/missing, keep it
	if (delta==0) return;
	double gq;
	if (!read_val(GQ,gq)) return;
	if (gq<GTP::GQ_cut) { db[g][i]=-1; return; }
	int gt;
	if (!read_val_ge(genotype,gt,0)) return;
	if (ploidy==1) gt*=2;
	if (gt==2) return;
	if (delta>0 && gt<2) return;
	if (delta<0 && gt>2) return;
	db[g][i] = gt;
	if (logopen) logfile<<loginfo<<record[i]<<'\t'<<db[g][i]<<endl;
}

void populate_genotype_G1K(vector< vector<int> >& db, int delta, tabular_file& in, field_numbers& FldSpl, vector<int>& SexSpl, int chr_num, int bp, size_t g)
{
	for (size_t i=0;i<FldSpl.size();++i)
	{
		int ploidy = genepi::expected_ploidy(chr_num,bp,SexSpl[i]);
		char c = GTP::read(in[FldSpl[i]],chr_num,bp,SexSpl[i],gpar);
		if (!GTP::valid(c)) { db[g][i]=-2; continue; } // NA
		if (!GTP::usable(c)) { db[g][i]=-1; continue; } // missing
		if (db[g][i]<0) continue; // if already NA/missing, keep it
		int gt = GTP::num_alt_recessive(c);
		if (gt==0) continue;
		if (delta==0) continue;
		if (ploidy==1 && gt==1) gt=2;
		
		db[g][i] = 2 + gt * delta;
		if (logopen && db[g][i]!=2) logfile<<loginfo<<record[i]<<'\t'<<db[g][i]<<endl;
		if (db[g][i]<0) db[g][i]=0; // hope this won't happen
	}
}

// delta has special meaning: 1=partial_cover 2=whole_gene
void populate_genotype_XHMM(vector< vector<int> >& db, int delta, tabular_file& in, field_numbers& FldSpl, vector<int>& SexSpl, int chr_num, int bp, size_t g)
{
	int svtype=delta;
	for (size_t i=0;i<FldSpl.size();++i)
	{
		int ploidy = genepi::expected_ploidy(chr_num,bp,SexSpl[i]);
		if (ploidy==0) { db[g][i]=-2; continue; } // NA
		vector<string> fd, eq, sq, nq; // all fields, EQ, NQ
		boost::split(fd,in[FldSpl[i]],boost::is_any_of(":"));
		int gt; if (!read_val_ge_le(fd[0],gt,0,2))	continue;
		boost::split(eq,fd[3],boost::is_any_of(","));
		boost::split(sq,fd[4],boost::is_any_of(","));
		boost::split(nq,fd[5],boost::is_any_of(","));
		int e1; if (!read_val_ge(eq[0],e1,0)) exit_error("cannot read EQ for A1");
		int e2; if (!read_val_ge(eq[1],e2,0)) exit_error("cannot read EQ for A2");
		int s1; if (!read_val_ge(sq[0],s1,0)) exit_error("cannot read SQ for A1");
		int s2; if (!read_val_ge(sq[1],s2,0)) exit_error("cannot read SQ for A2");
		int n1; if (!read_val_ge(nq[0],n1,0)) exit_error("cannot read NQ for A1");
		int n2; if (!read_val_ge(nq[1],n2,0)) exit_error("cannot read NQ for A2");
		int NAR = 0; // number of alt recessive model
		if (gt==0) { if (n1<GTP::GQ_cut||n2<GTP::GQ_cut) { Tx_DEL[g][i]=-1; Tx_DUP[g][i]=-1; continue; } else continue; }
		if (gt==1 && !exist_element(keepSV,"DEL")) continue;
		if (gt==2 && !exist_element(keepSV,"DUP")) continue;
		if (gt==2 && svtype==1 && !pdupid) continue;
		if (gt==1) 				{ vector< vector<int> >& db=Tx_DEL; delta=-1; if (s1<GTP::GQ_cut) { db[g][i]=-1; continue; } else if (ploidy==1) NAR=2; else NAR=1; }
		if (gt==2 && svtype==1) { vector< vector<int> >& db=Tx_DEL; delta=-1; if (s2<GTP::GQ_cut) { db[g][i]=-1; continue; } else if (ploidy==1) NAR=2; else NAR=1; }
		if (gt==2 && svtype==2) { vector< vector<int> >& db=Tx_DUP; delta= 1; if (s2<GTP::GQ_cut) { db[g][i]=-1; continue; } else if (ploidy==1) NAR=2; else NAR=1; }
		if (db[g][i]<0) continue; // if already NA/missing, keep it
		db[g][i] = 2 + NAR * delta;
		if (logopen && db[g][i]!=2) logfile<<loginfo<<record[i]<<'\t'<<db[g][i]<<endl;
		if (db[g][i]<0) db[g][i]=0; // hope this won't happen
	}
}

void populate_genotype_XHMM2(vector< vector<int> >& db, int delta, tabular_file& in, field_numbers& FldSpl, vector<int>& SexSpl, int chr_num, int bp, size_t g)
{
	for (size_t i=0;i<FldSpl.size();++i)
	{
		int ploidy = genepi::expected_ploidy(chr_num,bp,SexSpl[i]);
		if (ploidy==0) { db[g][i]=-2; continue; } // NA
		vector<string> fd, eq, sq, nq; // all fields, EQ, NQ
		boost::split(fd,in[FldSpl[i]],boost::is_any_of(":"));
		int gt; if (!read_val_ge_le(fd[0],gt,0,2))	continue;
		boost::split(eq,fd[3],boost::is_any_of(","));
		boost::split(sq,fd[4],boost::is_any_of(","));
		boost::split(nq,fd[5],boost::is_any_of(","));
		int e1; if (!read_val_ge(eq[0],e1,0)) exit_error("cannot read EQ for A1");
		int e2; if (!read_val_ge(eq[1],e2,0)) exit_error("cannot read EQ for A2");
		int s1; if (!read_val_ge(sq[0],s1,0)) exit_error("cannot read SQ for A1");
		int s2; if (!read_val_ge(sq[1],s2,0)) exit_error("cannot read SQ for A2");
		int n1; if (!read_val_ge(nq[0],n1,0)) exit_error("cannot read NQ for A1");
		int n2; if (!read_val_ge(nq[1],n2,0)) exit_error("cannot read NQ for A2");
		int NAR = 0; // number of alt recessive model
		int s=std::max(s1,s2);
		//cerr<<delta<<' '<<gt<<' '<<e1<<' '<<e2<<' '<<s1<<' '<<s2<<' '<<n1<<' '<<n2<<' '<<db[g][i]<<endl;
		if (gt==1)
		{
			if		(delta==-1) { if (s<GTP::GQ_cut) { db[g][i]=-1; continue; } else if (ploidy==1) NAR=2; else NAR=1; } // prv s1<GTP::GQ_cut
			else if (delta== 1) { if (s<GTP::GQ_cut) { db[g][i]=-1; continue; } else if (ploidy==1) NAR=2; else NAR=1; } // prv s2<GTP::GQ_cut
		}
		//cerr<<' '<<NAR<<endl;
		if (gt==0) continue; // { if (n1<GTP::GQ_cut||n2<GTP::GQ_cut) { db[g][i]=-1; continue; } else continue; }
		if (db[g][i]<0) continue; // if already NA/missing, keep it
		db[g][i] = 2 + NAR * delta;
		if (logopen && db[g][i]!=2) logfile<<loginfo<<record[i]<<'\t'<<db[g][i]<<endl;
		if (db[g][i]<0) db[g][i]=0; // hope this won't happen
		//cerr<<' '<<db[g][i]<<'\n';
	}
}

void (*populate_genotype)(vector< vector<int> >& db, int delta, tabular_file& in, field_numbers& FldSpl, vector<int>& SexSpl, int chr_num, int bp, size_t g)=NULL;

double pen(int copies)
{
	switch (copies)
	{
		case 0: return perch::penetr[2]; break;
		case 1: return perch::penetr[1]; break;
		case 2: return perch::penetr[0]; break;
		case 3: return perch::penetr[1]; break;
		case 4: return perch::penetr[2]; break;
		default:return perch::penetr[2]; break;
	}
	exit_error("impossible");
	return perch::penetr[2];
}

void test(vector< vector<int> >& db)
{
	for (size_t gene=0;gene<TxInfo.size();++gene)
	{
		// convert data, exclude records with missing genotypes
		struct CGinfo {				// collapsed genotype info
			int					cs;	// number of cases
			int					ct;	// number of controls
			CGinfo():cs(0),ct(0) {}
		};
		struct STinfo {				// strata info
			int					cs; // number of cases
			int					ct; // number of controls
			map<int,CGinfo>		cg;	// collapsed genotype
			STinfo():cs(0),ct(0) {}
		};
		map<double, STinfo>		sg;	// stratified genotypes sg[strata] (strata defined by relative risk)
		double t_cs=0, t_ct=0;		// total number of cases, controls w/o missing
		double m_cs=0, m_ct=0;		// total number of cases, controls w/  missing
		double a_cs=0, a_ct=0;		// total number of cases, controls w/  alternative allele (!=2)
		for (size_t i=0;i<db[gene].size();++i)
		{
			if (db[gene][i]==-2) continue;
			if (db[gene][i]==-1) {
				if (DepVar[i])	++m_cs;
				else			++m_ct;
				continue;
			}
			STinfo& s = sg[strata[i]];
			if (DepVar[i])	{ ++s.cg[db[gene][i]].cs; ++s.cs; ++t_cs; if (db[gene][i]!=2) ++a_cs; }
			else			{ ++s.cg[db[gene][i]].ct; ++s.ct; ++t_ct; if (db[gene][i]!=2) ++a_ct; }
		}
		
		if (!t_cs || !t_ct)
		{
			// std::numeric_limits<double>::signaling_NaN();
		}
		else
		{
			// QC
			double PrpCSs = m_cs/(m_cs+t_cs);
			double PrpCTs = m_ct/(m_ct+t_ct);
			double PrpAll = (m_cs+m_ct) / (m_cs+m_ct+t_cs+t_ct);
			double zScore = (PrpCSs-PrpCTs) / sqrt(PrpAll*(1-PrpAll)*(1.0/(m_cs+t_cs)+1.0/(m_ct+t_ct)));
			double pValue = cdf_norms_2sided_pv(zScore);
			if (perch::MisCut!=1)
			{
				if (perch::Mis_ea)
				{
					if (PrpCSs>perch::MisCut||PrpCTs>perch::MisCut) continue;
				}
				else
				{
					if (PrpAll>perch::MisCut) continue;
				}
			}
			if (MisPVl && pValue<=MisPVl) continue;
			if (filSAF && (a_cs+a_ct)/(t_cs+t_ct)>=filSAF ) continue;

			// calculate
			bool calculated=false;
			double score=0;	// log10LR
			for (auto &s:sg)
			{
				if (!s.second.cs || !s.second.ct) continue;
				if (s.second.cg.size()<2) continue;
				calculated=true;
				double newPre = truncated(perch::preval * s.first); // prevalence in this strata: P(s=1|c)
				// cerr<<perch::preval<<'\t'<<s.first<<endl;
				vector<double> fp,pn,fs,ft,ns,nt;
				for (auto &g:s.second.cg)
				{
					ns.push_back( g.second.cs);
					nt.push_back( g.second.ct);
					pn.push_back( pen(g.first) );
					fp.push_back( (double)g.second.cs/(double)s.second.cs * newPre + (double)g.second.ct/(double)s.second.ct * (1-newPre) );
					// cerr<<g.first<<'\t'<<ns.back()<<'\t'<<nt.back()<<'\t'<<fp.back()<<'\t'<<pn.back()<<endl;
				}
				MakeFraction(fp);
				genepi::up_dn_af(fp,pn,fs,ft);
				double llr=0; // log10 likelihood ratio (LOD) for this strata
				for (size_t i=0;i<fp.size();++i)
				{
					llr +=		  ns[i]  * log10(fs[i]) ;
					llr +=		  nt[i]  * log10(ft[i]) ;
					llr -= (ns[i]+nt[i]) * log10(fp[i]) ;
					// cerr<<ns[i]<<'/'<<s.second.cs<<'\t'<< nt[i]<<'/'<<s.second.ct<<'\t'<<fp[i]<<'\t'<<fs[i]<<'\t'<<ft[i]<<'\t'<<llr<<endl;
				}
				score += llr;
			}
			if (calculated) program.outf << TxInfo[gene].name2 << '\t' << score << endl;
		}
	}
}

int main (int argc, char * const argv[])
{
	// basic parameters
	const int SV_GQ_CUT=60;
	
	GTP::GQ_cut=SV_GQ_CUT;
	GTP::DP_cut=0;
	perch::prelude(argc,argv);
	
	// other parameters
	string			spl_in;				// sample file with header: SeqID Aff/QT Cov's
	string			vcf_in;				// input genotype file in VCF format
	string			del_lg;				// output log file for damaged genes.   col:                    Symbol SeqID Copy
	string			det_lg;				// output log file for variant details. col: Chr Start End Type Symbol SeqID Copy
	bool			clamms = false;		// input is CLAMMs format
	bool			inLong = false;		// input is long format
	bool			penncn = false;		// input is PennCNV format
	double			filVQS = 1.2168;	// filter variants by VQSLOD>filVQS. Always filter, even when filVQS=0.
	bool			ToWait = false;		// wait for vSIM to print the first comment line of the genotype file
	string			testWh = "no";		// test what
	int				up_reg = 35;		// upstream flanking region
	int				dn_reg = 0;			// downstream flanking region
	int				in_reg = 12;		// intronic flanking region
	
	// handle program arguments
	genepi::canonic=true; // use canonical transcript only
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"--spl"))				ReadArg(program.arg(),argi,spl_in);
		else if (str_startsw(program.arg()[argi],"--filt-vqs"))			ReadArg(program.arg(),argi,filVQS);
		else if (str_startsw(program.arg()[argi],"--test"))				ReadArg(program.arg(),argi,testWh);
		else if (str_startsw(program.arg()[argi],"--filt-miss-pv"))		ReadArg(program.arg(),argi,MisPVl);
		else if (str_startsw(program.arg()[argi],"--filt-SplAF"))	{	ReadArg(program.arg(),argi,filSAF); if (filSAF>=0.5) filSAF=0; }
		else if (str_startsw(program.arg()[argi],"--up"))				ReadArg(program.arg(),argi,up_reg);
		else if (str_startsw(program.arg()[argi],"--dn"))				ReadArg(program.arg(),argi,dn_reg);
		else if (str_startsw(program.arg()[argi],"--intron"))			ReadArg(program.arg(),argi,in_reg);
		else if (str_startsw(program.arg()[argi],"--long"))				ReadArg(program.arg(),argi,inLong);
		else if (str_startsw(program.arg()[argi],"--clamms"))		{	ReadArg(program.arg(),argi,clamms); if (GTP::GQ_cut==SV_GQ_CUT) GTP::GQ_cut=2; }
		else if (str_startsw(program.arg()[argi],"--penn"))			{	ReadArg(program.arg(),argi,penncn); if (GTP::GQ_cut==SV_GQ_CUT) GTP::GQ_cut=10; }
		else if (str_startsw(program.arg()[argi],"--sv-type"))			ReadSet(program.arg(),argi,keepSV);
		else if (str_startsw(program.arg()[argi],"--out-lof"))			ReadArg(program.arg(),argi,del_lg);
		else if (str_startsw(program.arg()[argi],"--out-LoF"))			ReadArg(program.arg(),argi,del_lg);
		else if (str_startsw(program.arg()[argi],"--detail"))			ReadArg(program.arg(),argi,det_lg);
		else if (str_startsw(program.arg()[argi],"--pDup-damage"))		ReadArg(program.arg(),argi,pdupid);
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		//else if (spl_in.empty()) spl_in=program.arg()[argi];
		else if (vcf_in.empty()) vcf_in=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help text
	program.help_text_var("_Default_sample_file",spl_in);
	program.help_text_var("_Default_filt_VQS",ftos(filVQS));
	program.help_text_var("_Default_test",testWh);
	program.help_text_var("_Default_filt_miss_d",ftos(MisPVl));
	program.help_text_var("_Default_filt_SplAF",ftos(filSAF));
	program.help_text_var("_Default_up",itos(up_reg));
	program.help_text_var("_Default_dn",itos(dn_reg));
	program.help_text_var("_Default_intron",itos(in_reg));
	program.help_text_var("_Default_clamms",str_YesOrNo(clamms));
	program.help_text_var("_Default_long",str_YesOrNo(inLong));
	program.help_text_var("_Default_penn",str_YesOrNo(penncn));
	program.help_text_var("_Default_del_log",del_lg);
	program.help_text_var("_Default_det_log",det_lg);
	program.help_text_var("_Default_sv_type",str_of_container(keepSV,',',false));
	program.help_text_var("_Default_pDupDamage",str_YesOrNo(pdupid));
	program.push_back_help(GTP::help_text());
	perch::check_arguments();
	
	// check parameters
	if (perch::penetr.empty())	exit_error("penetrance is not set.");
	if (perch::preval==0)		exit_error("prevalence is not set");
	if (spl_in.empty()) exit_error("Sample File not set.");
	int	DoWhat = -1;
	boost::to_lower(testWh);
	if		(testWh=="no")	DoWhat=0;
	else if	(testWh=="damaging") DoWhat=1;
	else if	(testWh=="enhancing") DoWhat=2;
	//else if (testWh=="all") DoWhat=3; // removed because the test doesn't make sense. It's better to do linear regression.
	else exit_error("--test consequence should be del/dup/all/no.");
	bool log_only = (DoWhat==0 && !del_lg.empty());

	// read gene db
	genepi::read_genes();

	// read sample file
	map< string, double >			dep_db;		// dep_db[SeqID]           = dependent variable
	map< string, vector<double> >	cov_db;		// cov_db[SeqID][variable] = covariate variable
	map< string, int >				SexMap;		// SexMap[SeqID]           = gender (1 for male, 2 for female)
	set< string >					h_csID;		// SeqID of cases
	set< string >					h_ctID;		// SeqID of controls
	double							factor=1;	// 1/denominator for RRc calculation
	if (!spl_in.empty())
	{
		int			cols=0;		// number of columns
		const int	covBgn=3;	// covariates start from this column (0-based)
		const int	minCol=3;	// minimum number of columns in a sample file
		const int	Spl_ID=0;
		const int	SplSex=1;
		const int	SplAff=2;
		
		if (ToWait)
		{
			string first_line;
			safeGetline(cin,first_line);
			if (first_line!="## vSIM created genotype file")
				exit_error("Don't use --wait if genotype file is not in standard input or not created by vSIM.");
		}
		for (Rows_in_File(in,spl_in,0)) // required columns: SeqID Outcome
		{
			// read and clean data
			if (in.RowNumber()==0)
			{
				cols=in.NumFields();
				if (cols<minCol)															exit_error("Insufficient number of columns in the Sample File "+spl_in);
				boost::to_lower(in[Spl_ID]); if (in[Spl_ID]!="seqid"&&in[Spl_ID]!="sample")	exit_error("The first  column of a Sample File should be SeqID/sample.");
				boost::to_lower(in[SplSex]); if (in[SplSex]!="sex"&&in[SplSex]!="gender")	exit_error("The second column of a Sample File should be sex/gender.");
				for (int c=covBgn;c<cols;++c) { double OR; if (ReadStr(substr_after_find(in[c],"_OR="),OR,0)) CovORs.push_back(OR); else CovORs.push_back(0); }
				continue;
			}
			if (cols!=in.NumFields()) exit_error("inconsistent number of columns in "+spl_in);
			if (exist_element(SexMap,in[Spl_ID])) exit_error("duplicated sample "+in[Spl_ID]+" in "+spl_in);
			if (exist_element(perch::rm_ind,in[Spl_ID])) continue; // skip samples to be removed

			int		sex = perch::read_sex(in[SplSex]);
			double	dep = perch::read_aff(in[SplAff]); if (log_only) dep=1;
			vector<double>	cov;	// covariate variable
			SexMap[in[Spl_ID]]=sex;
			if (std::isnan(dep)) continue; // skip samples with missing outcome
			
			// read covariates
			if (cols>covBgn)
			{
				bool mss=false;	// this individual is missing for at least 1 covariate
				for (int c=covBgn;c<cols;++c)
				{
					double val;
					if (!ReadStr(in[c],val,0)) { mss=true; break; }
					cov.push_back(val);
				}
				if (mss) continue; // skip samples with 1+ missing covariate
			}
			
			// update
			if		(dep==2) h_csID.insert(in[Spl_ID]);
			else if (dep==1) h_ctID.insert(in[Spl_ID]);
			dep_db[in[Spl_ID]]=dep;
			cov_db[in[Spl_ID]]=cov;
		}
		if (!perch::is_qtl())
			for (auto &i:dep_db) i.second-=1;

		// prepare covariates for HLR, which only works for case-control data
		if (!h_csID.empty() && !h_ctID.empty())
		{
			// do logistic regression to obtain ORs
			{
				const int ni = h_csID.size() + h_ctID.size();	// number of individuals (observations)
				const int np = cov_db.begin()->second.size()+1;	// number of parameters
				Eigen::VectorXd c(np);			// beta
				Eigen::MatrixXd cov(np,np);		// covariance matrix
				Eigen::VectorXd y(ni);			// y
				Eigen::MatrixXd X(ni,np);		// X
				size_t j=0;
				for (auto &id:h_csID)
				{
					y(j)=dep_db[id];
					X(j,0)=1;
					for (size_t c=0;c<cov_db[id].size();++c) X(j,c+1)=cov_db[id][c];
					++j;
				}
				for (auto &id:h_ctID)
				{
					y(j)=dep_db[id];
					X(j,0)=1;
					for (size_t c=0;c<cov_db[id].size();++c) X(j,c+1)=cov_db[id][c];
					++j;
				}
				if (multifit_logistic_validate(X,y,1))
				{
					double chisq;
					multifit_logistic_workspace work(ni,np);
					multifit_logistic(X,y,c,cov,chisq,work);
					for(int i=1;i<np;++i) { CovORs[i-1]=exp(c(i)); }
				}
				else
					exit_error("Failed to do logistic regression on covariates only.");
			}
			
			// prepare covariates for HLR, assumes independence between covariates
			factor = 1;
			for (size_t c=0;c<CovORs.size();++c)
			{
				map<double,pair<double,double> > ValCnt; // ValCnt[value] = <#cs,#ct>
				for (auto &i:cov_db)
				{
					if (dep_db[i.first]) { ++ValCnt[i.second[c]].first;  }
					else				 { ++ValCnt[i.second[c]].second; }
				}
				double denominator = 0;
				double sum_of_prob = 0;
				for (auto &v:ValCnt)
				{
					double prob = v.second.first / h_csID.size() * perch::preval + v.second.second / h_ctID.size() * (1-perch::preval); // P(value)
					sum_of_prob += prob;
					denominator += prob * pow(CovORs[c],v.first);
				}
				denominator /= sum_of_prob;
				factor /= denominator;
			}
		}
	}
	
	if (!det_lg.empty())
	{
		openOutFile(logfile,det_lg);
		logopen=true;
		logfile<<"#CHROM"<<DLMTR<<"StartBP"<<DLMTR<<"EndBP"<<DLMTR<<"Type"<<DLMTR<<perch::h_symb<<DLMTR<<"SeqID"<<DLMTR<<"Copy"<<endl;
	}
	int ErrMsg1 = elog.get_token(" variant skipped because END is not an integer.");
	int ErrMsg2 = elog.get_token(" variant skipped because of an unknown SVTYPE.");
	if (inLong) // long format
	{
		map<string,int> coordinate; // coordinate[sampleID] = coordinate in DepVar etc.
		for (auto &x:dep_db) if (x.second!=0) { CovVar.push_back(cov_db[x.first]); DepVar.push_back(dep_db[x.first]); record.push_back(x.first); } // cs
		for (auto &x:dep_db) if (x.second==0) { CovVar.push_back(cov_db[x.first]); DepVar.push_back(dep_db[x.first]); record.push_back(x.first); } // ct
		for (size_t i=0;i<DepVar.size();++i)
		{
			double RRc = factor; // Relative Risk as compared to the general population
			for (size_t c=0;c<CovORs.size();++c) RRc *= pow(CovORs[c],CovVar[i][c]);
			strata.push_back(boost::lexical_cast<double>(ftos(RRc,1))); // use 1 significant digits to merge as many samples as possible
		}
		for (auto &g:genepi::gene_byTranscript) for (auto &t:g.second)
		{
			Tx_DEL.push_back(vector<int>(DepVar.size(),2));
			Tx_DUP.push_back(vector<int>(DepVar.size(),2));
			TxInfo.push_back(t.second);
		}
		for (size_t i=0;i<record.size();++i) coordinate[record[i]]=i;
		
		int ColInd=-1;	// field numb for ID
		int ColChr=-1;	// field numb for Chr
		int ColPos=-1;	// field numb for Start
		int ColEnd=-1;	// field numb for End
		int ColTyp=-1;	// field numb for Type
		int ColQua=-1;	// field numb for Qual
		tfile_format	format;
		format.set_option(SKIP_NOTES,false);
		bool header_not_read = true;
		for (Rows_in_File(in,vcf_in,&format))
		{
			// read header
			if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
			{
				in.clear_nf();
				continue;
			}
			if (exist_any(perch::h_col1, in.contents()))
			{
				header_not_read=false;
				format.clear_field_nums();
				ColInd=-1;
				ColChr=-1;
				ColPos=-1;
				ColEnd=-1;
				ColTyp=-1;
				ColQua=-1;
				for (int i=0;i<in.NumFields();++i)
				{
					string s=boost::to_lower_copy(in[i]);
					if (s=="sampleid" || s=="sample_id" || s=="sample" || s=="seqid" || s=="seq_id" || s=="indid" || s=="ind_id" || s=="id" || s=="iid" || s=="ind") ColInd=i;
					if (s=="type")	ColTyp=i;
					if (s=="qual")	ColQua=i;
					if (s=="end")	ColEnd=i;
					if (s=="chr"   && ColChr==-1) ColChr=i; if (s=="#chrom")	ColChr=i;
					if (s=="start" && ColPos==-1) ColPos=i; if (s=="pos")		ColPos=i;
				}
				if (ColInd==-1) exit_error("The ID column is missing.");		format.set_field_nums(ColInd,"lack the ID    field",tfile_format::Break);
				if (ColChr==-1) exit_error("The #CHROM/Chr column is missing.");format.set_field_nums(ColChr,"lack the Chr   field",tfile_format::Break);
				if (ColPos==-1) exit_error("The POS/Start column is missing.");	format.set_field_nums(ColPos,"lack the Start field",tfile_format::Break);
				if (ColEnd==-1) exit_error("The END/End column is missing.");	format.set_field_nums(ColEnd,"lack the End   field",tfile_format::Break);
				if (ColTyp==-1) exit_error("The Type column is missing.");		format.set_field_nums(ColTyp,"lack the Type  field",tfile_format::Break);
				if (ColQua==-1) exit_error("The Qual column is missing.");		format.set_field_nums(ColQua,"lack the Qual  field",tfile_format::Break);
				continue;
			}
			if (header_not_read) exit_error("Header lines missing.");
			
			if (!exist_element(SexMap,in[ColInd])) continue; // assume dep_db and cov_db have the same number of samples as SexMap
			string SVTYPE = boost::to_upper_copy(in[ColTyp]);
			string END = in[ColEnd];
			int	chr_num;	if (!genepi::read_chr_num(in[ColChr],chr_num))	exit_error("Failed to read "+in[ColChr]+" as a chromosome.");
			int	bp;			if (!read_val_ge(in[ColPos],bp,1))				exit_error("Failed to read "+in[ColPos]+" as a position in basepairs.");
			if (SVTYPE=="INV")
			{
				if (!exist_element(keepSV,"INV")) continue;
				int b2; if (!ReadStr(END,b2,ErrMsg1)) continue;
				for (size_t g=0;g<TxInfo.size();++g)
					if ( genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg) &&
						!genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
					{
						loginfo = in[ColChr]+DLMTR+in[ColPos]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
						populate_genotype_long(Tx_DEL,-1,"1",in[ColQua],coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
					}
			}
			else if (SVTYPE=="INS"||SVTYPE=="SVA"||SVTYPE=="ALU"||SVTYPE=="LINE1")
			{
				if (!exist_element(keepSV,"INS")) continue;
				for (size_t g=0;g<TxInfo.size();++g)
					if (genepi::overlap_exon(TxInfo[g],chr_num,bp,bp,up_reg,dn_reg,in_reg))
					{
						loginfo = in[ColChr]+DLMTR+in[ColPos]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
						populate_genotype_long(Tx_DEL,-1,"1",in[ColQua],coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
					}
			}
			else if (SVTYPE=="DEL")
			{
				if (!exist_element(keepSV,"DEL")) continue;
				int b2; if (!ReadStr(END,b2,ErrMsg1)) continue;
				for (size_t g=0;g<TxInfo.size();++g)
					if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg))
					{
						loginfo = in[ColChr]+DLMTR+in[ColPos]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
						populate_genotype_long(Tx_DEL,-1,"1",in[ColQua],coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
					}
			}
			else if (SVTYPE=="DUP") // equivalent of <CN2> with genotype 0:1 in G1K format or <DUP> in XHMM format (3 copies in total)
			{
				if (!exist_element(keepSV,"DUP")) continue;
				int b2; if (!ReadStr(END,b2,ErrMsg1)) continue;
				if (pdupid)
				{
					for (size_t g=0;g<TxInfo.size();++g)
						if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg) && !genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
						{
							loginfo = in[ColChr]+DLMTR+in[ColPos]+DLMTR+END+DLMTR+"p"+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
							populate_genotype_long(Tx_DEL,-1,"1",in[ColQua],coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
						}
				}
				for (size_t g=0;g<TxInfo.size();++g)
					if (genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
					{
						loginfo = in[ColChr]+DLMTR+in[ColPos]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
						populate_genotype_long(Tx_DUP,1,"1",in[ColQua],coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
					}
			}
			else { elog.add(ErrMsg2); continue; }
		}
	}
	else if (clamms) // CLAMMs format
	{
		map<string,int> coordinate; // coordinate[sampleID] = coordinate in DepVar etc.
		for (auto &x:dep_db) if (x.second!=0) { CovVar.push_back(cov_db[x.first]); DepVar.push_back(dep_db[x.first]); record.push_back(x.first); } // cs
		for (auto &x:dep_db) if (x.second==0) { CovVar.push_back(cov_db[x.first]); DepVar.push_back(dep_db[x.first]); record.push_back(x.first); } // ct
		for (size_t i=0;i<DepVar.size();++i)
		{
			double RRc = factor; // Relative Risk as compared to the general population
			for (size_t c=0;c<CovORs.size();++c) RRc *= pow(CovORs[c],CovVar[i][c]);
			strata.push_back(boost::lexical_cast<double>(ftos(RRc,1))); // use 1 significant digits to merge as many samples as possible
		}
		for (auto &g:genepi::gene_byTranscript) for (auto &t:g.second)
		{
			Tx_DEL.push_back(vector<int>(DepVar.size(),2));
			Tx_DUP.push_back(vector<int>(DepVar.size(),2));
			TxInfo.push_back(t.second);
		}
		for (size_t i=0;i<record.size();++i) coordinate[record[i]]=i;
		
		/*
		 Columns (tab separated):
		 1: Chromosome
		 2: Start (5' exon breakpoint)
		 3: End (3' exon breakpoint)
		 4: SampleID
		 5: DEL/DUP
		 6: Copy Number
		 7: # of exons
		 8: Quality (confidence) of non-diploid call (Qnon_dip, phred scaled)
		 9: Quality (confidence) of exact copy number call (Qexact, not phred scaled)
		 10: # hom SNPs
		 11: # SNPs
		 12: Allele balance of heterozygous SNPs
		 13: QC score (0-3)
		 14: Locus definition tag (maps to 1st 3 columns of loci.annotated.bed file)
		 15: ENSEMBL75 overlapping genes (comma seperated)
		 */
		int ColChr=0;	// field numb for Chr
		int ColPos=1;	// field numb for Start
		int ColEnd=2;	// field numb for End
		int ColInd=3;	// field numb for ID
		int ColTyp=4;	// field numb for Type
		int ColCpy=5;	// field numb for Copy. In RGC: [0,1] for DEL; [2,6] for DUP. When Typ=DUP and Cpy=2, Chr is either X or Y.
		int ColQua=12;	// QC score. 0-1 is low confidence; 2-3 is high confidence. I know this because the RGC "*.high-confidence.cnv.annotated.bed" only have 2 and 3.
		for (Rows_in_File(in,vcf_in,6))
		{
			if (!exist_element(SexMap,in[ColInd])) continue; // assume dep_db and cov_db have the same number of samples as SexMap
			string SVTYPE = boost::to_upper_copy(in[ColTyp]);
			int	chr_num;	if (!genepi::read_chr_num(in[ColChr],chr_num))	exit_error("Failed to read "+in[ColChr]+" as a chromosome.");
			int	bp;			if (!read_val_ge(in[ColPos],bp,1))				exit_error("Failed to read "+in[ColPos]+" as a position in basepairs.");
			int	b2;			if (!read_val_ge(in[ColEnd],b2,1))				exit_error("Failed to read "+in[ColEnd]+" as a position in basepairs.");
			if (SVTYPE=="DEL")
			{
				if (!exist_element(keepSV,"DEL")) continue;
				for (size_t g=0;g<TxInfo.size();++g)
					if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg))
					{
						loginfo = in[ColChr]+DLMTR+in[ColPos]+DLMTR+in[ColEnd]+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
						populate_genotype_clamms(Tx_DEL,-1,in[ColCpy],in[ColQua],coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
					}
			}
			else if (SVTYPE=="DUP")
			{
				if (!exist_element(keepSV,"DUP")) continue;
				if (pdupid)
				{
					for (size_t g=0;g<TxInfo.size();++g)
						if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg) && !genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
						{
							loginfo = in[ColChr]+DLMTR+in[ColPos]+DLMTR+in[ColEnd]+DLMTR+"p"+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
							populate_genotype_clamms(Tx_DEL,-1,in[ColCpy],in[ColQua],coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
						}
				}
				for (size_t g=0;g<TxInfo.size();++g)
					if (genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
					{
						loginfo = in[ColChr]+DLMTR+in[ColPos]+DLMTR+in[ColEnd]+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
						populate_genotype_clamms(Tx_DUP,1,in[ColCpy],in[ColQua],coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
					}
			}
			else { elog.add(ErrMsg2); continue; }
		}
	}
	else if (penncn)
	{
		map<string,int> coordinate; // coordinate[sampleID] = coordinate in DepVar etc.
		for (auto &x:dep_db) if (x.second!=0) { CovVar.push_back(cov_db[x.first]); DepVar.push_back(dep_db[x.first]); record.push_back(x.first); } // cs
		for (auto &x:dep_db) if (x.second==0) { CovVar.push_back(cov_db[x.first]); DepVar.push_back(dep_db[x.first]); record.push_back(x.first); } // ct
		for (size_t i=0;i<DepVar.size();++i)
		{
			double RRc = factor; // Relative Risk as compared to the general population
			for (size_t c=0;c<CovORs.size();++c) RRc *= pow(CovORs[c],CovVar[i][c]);
			strata.push_back(boost::lexical_cast<double>(ftos(RRc,1))); // use 1 significant digits to merge as many samples as possible
		}
		for (auto &g:genepi::gene_byTranscript) for (auto &t:g.second)
		{
			Tx_DEL.push_back(vector<int>(DepVar.size(),2));
			Tx_DUP.push_back(vector<int>(DepVar.size(),2));
			TxInfo.push_back(t.second);
		}
		for (size_t i=0;i<record.size();++i) coordinate[record[i]]=i;
		
		int ColInd=4;	// field numb for ID
		tfile_format	format;
		format.set_option(SKIP_NOTES,false);
		format.set_option(SUCCESSIVE_DELIMITERS_AS_ONE,true);
		for (Rows_in_File(in,vcf_in,&format))
		{
			if (!exist_element(SexMap,in[ColInd])) continue; // assume dep_db and cov_db have the same number of samples as SexMap
			int copy; if (!read_val(in[3].substr(10),copy)) continue;
			string qual_str = in[7].substr(5);
			string chr_str = substr_before_find(in[0],":");
			string bp_str = substr_before_find(substr_after_find(in[0],":"),"-");
			string b2_str = substr_after_find(substr_after_find(in[0],":"),"-");
			int bp; if (!read_val(bp_str,bp)) continue;
			int b2; if (!read_val(b2_str,b2)) continue;
			int chr_num = genepi::read_chr_num(chr_str); if (chr_num<=0) continue;
			int delta;
			string geno;
			if		(copy==0) { delta=-1; geno="2"; if (!exist_element(keepSV,"DEL")) continue; }
			else if (copy==1) { delta=-1; geno="1"; if (!exist_element(keepSV,"DEL")) continue; }
			else if (copy==3) { delta= 1; geno="1"; if (!exist_element(keepSV,"DUP")) continue; }
			else if (copy==4) { delta= 1; geno="2"; if (!exist_element(keepSV,"DUP")) continue; }
			else continue;
			if (delta==-1)
			{
				for (size_t g=0;g<TxInfo.size();++g)
					if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg))
					{
						loginfo = chr_str+DLMTR+bp_str+DLMTR+b2_str+DLMTR+"DEL"+DLMTR+TxInfo[g].name2+DLMTR;
						populate_genotype_long(Tx_DEL,delta,geno,qual_str,coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
					}
			}
			if (delta==1)
			{
				if (pdupid)
				{
					for (size_t g=0;g<TxInfo.size();++g)
						if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg) && !genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
						{
							loginfo = chr_str+DLMTR+bp_str+DLMTR+b2_str+DLMTR+"pDUP"+DLMTR+TxInfo[g].name2+DLMTR;
							populate_genotype_long(Tx_DEL,-1,geno,qual_str,coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
						}
				}
				for (size_t g=0;g<TxInfo.size();++g)
					if (genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
					{
						loginfo = chr_str+DLMTR+bp_str+DLMTR+b2_str+DLMTR+"DUP"+DLMTR+TxInfo[g].name2+DLMTR;
						populate_genotype_long(Tx_DUP,delta,geno,qual_str,coordinate[in[ColInd]],SexMap[in[ColInd]],chr_num,bp,g);
					}
			}
		}
	}
	else // VCF or XHMM format
	{
		vector<int>		SexCSs;				// same order as FldCSs
		vector<int>		SexCTs;				// same order as FldCTs
		vector<int>		SexSpl;				// same order as FldCSs + FldCTs
		field_numbers	FldCSs(false,true);	// field numb for cases
		field_numbers	FldCTs(false,true);	// field numb for controls
		field_numbers	FldSpl(false,true);	// field numb for cases + controls
		field_numbers	FldIdx(false,true);	// field numb for Chr,Start,Ref,Alt / #CHROM,POS,REF,ALT
		field_numbers	FldChr(false,true);	// field numb for #CHROM
		field_numbers	FldPos(false,true);	// field numb for POS
		field_numbers	FldRef(false,true);	// field numb for REF
		field_numbers	FldAlt(false,true);	// field numb for ALT
		field_numbers	FldFlt(false,true);	// field numb for FILTER
		field_numbers	FldInf(false,true);	// field numb for INFO
		field_numbers	FldFmt(false,true);	// field numb for FORMAT
		tfile_format	format;
		format.set_delimiters("\t");
		format.set_option(SKIP_NOTES,false);
		bool header_not_read = true;
		
		for (Rows_in_File(in, vcf_in, &format))
		{
			// read header
			if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
			{
				in.clear_nf();
				continue;
			}
			if (exist_any(perch::h_col1, in.contents()))
			{
				if (!header_not_read) { continue; }
				header_not_read=false;
				format.clear_field_nums();
				FldCSs.clear();
				FldCTs.clear();
				FldSpl.clear();
				FldIdx.clear();
				FldChr.clear();
				FldPos.clear();
				FldRef.clear();
				FldAlt.clear();
				FldFlt.clear();
				FldInf.clear();
				FldFmt.clear();
				for (int i=0;i<in.NumFields();++i)
				{
					if (exist_element(h_csID,in[i])) {	FldCSs.push_back(i+1); SexCSs.push_back(SexMap[in[i]]); }
					if (exist_element(h_ctID,in[i])) {	FldCTs.push_back(i+1); SexCTs.push_back(SexMap[in[i]]); }
					if (in[i]=="FILTER")				FldFlt.push_back(i+1);
					if (in[i]=="INFO")					FldInf.push_back(i+1);
					if (in[i]=="FORMAT")				FldFmt.push_back(i+1);
					if (in[i]=="Chr"   && FldIdx.no_input()) {	FldIdx.push_back(i+1); FldIdx.push_back(i+2); FldIdx.push_back(i+4); FldIdx.push_back(i+5); }
					if (in[i]=="#CHROM")	{	FldIdx.clear(); FldIdx.push_back(i+1); FldIdx.push_back(i+2); FldIdx.push_back(i+4); FldIdx.push_back(i+5); }
					if (in[i]=="Chr"   && FldChr.no_input())	FldChr.push_back(i+1); if (in[i]=="#CHROM")		{	FldChr.clear(); FldChr.push_back(i+1); }
					if (in[i]=="Start" && FldPos.no_input())	FldPos.push_back(i+1); if (in[i]=="POS")		{	FldPos.clear(); FldPos.push_back(i+1); }
					if (in[i]=="Ref"   && FldRef.no_input())	FldRef.push_back(i+1); if (in[i]=="REF")		{	FldRef.clear(); FldRef.push_back(i+1); }
					if (in[i]=="Alt"   && FldAlt.no_input())	FldAlt.push_back(i+1); if (in[i]=="ALT")		{	FldAlt.clear(); FldAlt.push_back(i+1); }
				}
				if (FldIdx.no_input()) exit_error("The #CHROM-POS-ID-REF-ALT/Chr-Start-End-Ref-Alt column is missing.");
				if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
				if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
				if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
				if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
				if (FldFlt.no_input()) exit_error("The FILTER column is missing.");
				if (FldInf.no_input()) exit_error("The INFO column is missing.");
				if (FldCSs.no_input() && FldCTs.no_input()) exit_error("Found no sample columns.");
				if (TxInfo.empty())
				{
					SexSpl.reserve( SexCSs.size() + SexCTs.size() ); // preallocate memory
					SexSpl.insert( SexSpl.end(), SexCSs.begin(), SexCSs.end() );
					SexSpl.insert( SexSpl.end(), SexCTs.begin(), SexCTs.end() );
					FldSpl = FldCSs + FldCTs;
					for (size_t i=0;i<FldCSs.size();++i) { CovVar.push_back(cov_db[in[FldCSs[i]]]); DepVar.push_back(dep_db[in[FldCSs[i]]]); record.push_back(in[FldCSs[i]]); }
					for (size_t i=0;i<FldCTs.size();++i) { CovVar.push_back(cov_db[in[FldCTs[i]]]); DepVar.push_back(dep_db[in[FldCTs[i]]]); record.push_back(in[FldCTs[i]]); }
					for (size_t i=0;i<FldSpl.size();++i)
					{
						double RRc = factor; // Relative Risk as compared to the general population
						for (size_t c=0;c<CovORs.size();++c) RRc *= pow(CovORs[c],CovVar[i][c]);
						strata.push_back(boost::lexical_cast<double>(ftos(RRc,1))); // use 1 significant digits to merge as many samples as possible
					}
					for (auto &g:genepi::gene_byTranscript)
						for (auto &t:g.second)
						{
							Tx_DEL.push_back(vector<int>(FldSpl.size(),2));
							Tx_DUP.push_back(vector<int>(FldSpl.size(),2));
							TxInfo.push_back(t.second);
						}
					// cerr<<"# read "<<TxInfo.size()<<" transcripts.\n";
				}
				continue;
			}
			if (header_not_read) exit_error("Header lines missing.");
			
			// basic info
			string& ref=in[FldRef[0]];
			string& alt=in[FldAlt[0]];
			
			// check file
			int FileType=0; // 0=KG/WHAM 1=XHMM 2=XHMM2
			if (ref=="<DIP>" && (alt=="<DEL>" || alt=="<DUP>") ) // XHMM2
			{
				if (!str_startsw(in[FldFmt[0]],"GT:NDQ:DQ:EQ:SQ:NQ"))
					exit_error("I read XHMM files that FORMAT starts with GT:NDQ:DQ:EQ:SQ:NQ. Please contact me so that I can change my program to read your file.");
				FileType=2;
				populate_genotype=populate_genotype_XHMM2;
			}
			else if (ref=="<DIP>" && alt=="<DEL>,<DUP>") // XHMM
			{
				if (!str_startsw(in[FldFmt[0]],"GT:NDQ:DQ:EQ:SQ:NQ"))
					exit_error("I read XHMM files that FORMAT starts with GT:NDQ:DQ:EQ:SQ:NQ. Please contact me so that I can change my program to read your file.");
				FileType=1;
				populate_genotype=populate_genotype_XHMM;
			}
			else // KG or WHAM
			{
				FileType=0;
				populate_genotype=populate_genotype_G1K;
				if (alt.find(',')!=std::string::npos) exit_error("The input file has not been split by alternative alleles.");
				if (!FldFmt.no_input()) gpar.read(in[FldFmt[0]]);
			}
			
			// read INFO
			vector<string> INFO;
			boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));

			// skip by QC. Consider not to do it because MaxAF is not reliable.
			if (!FldInf.no_input())
			{
				string vQC = get_string(INFO,"VICTOR_QC");
				if (!vQC.empty() && vQC!="PASS") { continue; }
			} //*/
			
			// skip if not an SV, skip by SVTYPE
			string SVTYPE=get_string(INFO,"SVTYPE");
			if (SVTYPE.empty() || !str_has(alt,"<")) { continue; }
			string END = get_string(INFO,"END");
			
			// skip by FILTER
			if (!FldFlt.no_input())
				if (!perch::filflt.empty() && !exist_element(perch::filflt,in[FldFlt[0]])) { continue; }
			double VQSLOD = std::numeric_limits<double>::signaling_NaN();
			try { VQSLOD=get_value(INFO,"VQSLOD"); } catch (...) {}
			if (VQSLOD<filVQS) {  continue; }
			
			// read genotypes
			int bp=-1; try { bp = boost::lexical_cast<int>(in[FldPos[0]]); } catch (...) { exit_error("Failed to read "+in[FldPos[0]]+" as a position in basepairs."); }
			int chr_num = genepi::read_chr_num(in[FldChr[0]]);
			if (SVTYPE=="INV")
			{
				if (!exist_element(keepSV,"INV")) continue;
				int b2; if (!ReadStr(END,b2,ErrMsg1)) continue;
				for (size_t g=0;g<TxInfo.size();++g)
					if ( genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg) &&
						!genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
					{
						loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
						(*populate_genotype)(Tx_DEL,-1,in,FldSpl,SexSpl,chr_num,bp,g);
					}
			}
			else if (SVTYPE=="INS"||SVTYPE=="SVA"||SVTYPE=="ALU"||SVTYPE=="LINE1")
			{
				if (!exist_element(keepSV,"INS")) continue;
				for (size_t g=0;g<TxInfo.size();++g)
					if (genepi::overlap_exon(TxInfo[g],chr_num,bp,bp,up_reg,dn_reg,in_reg))
					{
						loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
						(*populate_genotype)(Tx_DEL,-1,in,FldSpl,SexSpl,chr_num,bp,g);
					}
			}
			else if (SVTYPE=="DEL"||SVTYPE=="DUP"||SVTYPE=="CNV")
			{
				int b2; if (!ReadStr(END,b2,ErrMsg1)) continue;
				if (FileType==2)
				{
					if (alt=="<DEL>")
					{
						if (!exist_element(keepSV,"DEL")) continue;
						for (size_t g=0;g<TxInfo.size();++g)
							if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg))
							{
								loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+"DEL"+DLMTR+TxInfo[g].name2+DLMTR;
								(*populate_genotype)(Tx_DEL,-1,in,FldSpl,SexSpl,chr_num,bp,g);
							}
					}
					else if (alt=="<DUP>")
					{
						if (!exist_element(keepSV,"DUP")) continue;
						if (pdupid)
						{
							for (size_t g=0;g<TxInfo.size();++g)
								if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg) && !genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
								{
									loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+"pDUP"+DLMTR+TxInfo[g].name2+DLMTR;
									(*populate_genotype)(Tx_DEL,-1,in,FldSpl,SexSpl,chr_num,bp,g);
								}
						}
						for (size_t g=0;g<TxInfo.size();++g)
							if (genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
							{
								loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+"DUP"+DLMTR+TxInfo[g].name2+DLMTR;
								(*populate_genotype)(Tx_DUP,1,in,FldSpl,SexSpl,chr_num,bp,g);
							}
					}
				}
				else if (FileType==1)
				{
					for (size_t g=0;g<TxInfo.size();++g)
					{
						if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg))
						{
							if (genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
							{
								loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
								(*populate_genotype)(Tx_DEL,2,in,FldSpl,SexSpl,chr_num,bp,g);
							}
							else
							{
								loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
								(*populate_genotype)(Tx_DEL,1,in,FldSpl,SexSpl,chr_num,bp,g);
							}
						}
					}
				}
				else if (alt=="<CN0>" || alt=="<DEL>") // corresponds to KG || WHAM, respectively
				{
					if (!exist_element(keepSV,"DEL")) continue;
					for (size_t g=0;g<TxInfo.size();++g)
						if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg))
						{
							loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
							(*populate_genotype)(Tx_DEL,-1,in,FldSpl,SexSpl,chr_num,bp,g);
						}
				}
				else // <CN#> (#>=2)
				{
					if (!exist_element(keepSV,"DUP")) continue;
					int dt=0;
					if (alt=="<DUP>") dt=2; // <DUP> in WHAM is equal to <CN2> in KG
					else if (!ReadStr(alt.substr(3,alt.size()-4),dt,0)) exit_error("ALT is not <CN#>.");
					--dt;
					if (pdupid)
					{
						for (size_t g=0;g<TxInfo.size();++g)
							if (genepi::overlap_exon(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg,in_reg) && !genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
							{
								loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+"p"+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
								(*populate_genotype)(Tx_DEL,-1,in,FldSpl,SexSpl,chr_num,bp,g);
							}
					}
					for (size_t g=0;g<TxInfo.size();++g)
						if (genepi::overlap_StoE(TxInfo[g],chr_num,bp,b2,up_reg,dn_reg))
						{
							loginfo = in[FldChr[0]]+DLMTR+in[FldPos[0]]+DLMTR+END+DLMTR+SVTYPE+DLMTR+TxInfo[g].name2+DLMTR;
							(*populate_genotype)(Tx_DUP,dt,in,FldSpl,SexSpl,chr_num,bp,g);
						}
					// if (genepi::overlap_exon() && !genepi::overlap_StoE()) (*populate_genotype)(Tx_DEL,in,FldSpl,SexSpl,chr_num,bp,g);
					// Removed because: 1) it's unlikely and not necesary; 2) needs to consider the boundary of a gene; 3) it creates DEL & DUP in the same person.
				}
			}
			else { elog.add(ErrMsg2); continue; }
		}
	}
	if (logopen) closefile(logfile);
	
	if (!del_lg.empty())
	{
		openOutFile_or_exit(file,del_lg);
		for (size_t tx=0; tx<Tx_DEL.size(); ++tx)
			for (size_t id=0; id<Tx_DEL[tx].size(); ++id)
				if (Tx_DEL[tx][id]==0||Tx_DEL[tx][id]==1)
					file << TxInfo[tx].name2 << '\t' << record[id] << '\t' << Tx_DEL[tx][id] << '\n';
		closefile(file);
	}
	
	// analysis
	if		(DoWhat==0) ;
	else if	(DoWhat==1) { test(Tx_DEL); }
	else if (DoWhat==2) { test(Tx_DUP); }
	else if (DoWhat==3)
	{
		vector< vector<int> > Tx_CMB(TxInfo.size(),vector<int>(record.size(),2));
		for (size_t gene=0;gene<TxInfo.size();++gene)
			for (size_t i=0;i<record.size();++i)
			{
				if		(Tx_DUP[gene][i]<0  || Tx_DEL[gene][i]<0)	Tx_CMB[gene][i] = -1;
				else												Tx_CMB[gene][i] = Tx_DUP[gene][i] + Tx_DEL[gene][i] - 2;
				// if (Tx_DEL[gene][i]!=2 || Tx_DUP[gene][i]!=2 || Tx_CMB[gene][i]!=2) cerr<<gene<<' '<<i<<' '<<Tx_DEL[gene][i]<<' '<<Tx_DUP[gene][i]<<' '<<Tx_CMB[gene][i]<<endl;
			}
		test(Tx_CMB);
	}
	else exit_error("wrong test method.");
	return 0;
}
