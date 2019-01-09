/*
 Changes from ppi-gba:
 1) add weight for seed genes
 */
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_mtjobs.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/logic/tribool.hpp>
#include <mutex>		// requires c++11. g++47 needs it! clang++4.0/g++44 doesn't.
#include <Eigen/Core>
#include <Eigen/Dense>
#include "victor_par.hpp"
#define EIGEN_NO_DEBUG	// disable Eigen range checking for faster access
#define NDEBUG 			// this is a c++ standard macro. https://stackoverflow.com/questions/32605838/error-eigen-library-on-linux

#include "gnGBA.hpp"
using namespace std;

typedef double			rslt_t;				// for results (prediction, rank, auc)
typedef float			cell_t;				// for matrix, use float (32bit) to reduce memory usage
typedef Eigen::MatrixXf EMtx_t;				// ..correspondingly
typedef Eigen::VectorXf EVec_t;				// ..correspondingly
EMtx_t							W, A;		// original matrix W (weights), D, I (identity). type match to cell_t (float)
std::map<std::string,size_t>	gene2coor;	// coordinate starts from 0
std::vector<std::string>		coor2gene;	// coordinate starts from 0
std::vector<std::set<size_t> >	geneGroup;	// coordinate grouped by user (normally by range of strength)
size_t num_genes;							// total number of genes in matrix
size_t MinSeeds=1;							// mininum known disease genes for cross-validation
size_t NumFolds=10;							// cross validation folds
bool NB_r1=false;							// InterConnectedness revision 1
bool NB_r3=false;							// InterConnectedness revision 3

// basic network measures for each gene
vector<rslt_t>			nm_k;			// degree
vector<rslt_t>			nm_s;			// strength = sum(Wij)
multimap<rslt_t,size_t>	s2co;			// strength to coordinate
rslt_t					minS =  std::numeric_limits<cell_t>::max();
rslt_t					maxS = -std::numeric_limits<cell_t>::max();
rslt_t					BinWidth=0.05;	// bin width to find resemblance for calibration

void BasicNetworkMeasures()
{
	nm_k.assign(num_genes,0); // degree
	nm_s.assign(num_genes,0); // strength = sum(Wij)
	
	for (size_t i=0;i<num_genes;++i)
		for (size_t j=0;j<num_genes;++j)
			if (i!=j && W(j,i)) { nm_s[i]+=W(j,i); ++nm_k[i]; }
	for (size_t i=0;i<num_genes;++i)
	{
		if (nm_s[i]<minS) minS=nm_s[i];
		if (nm_s[i]>maxS) maxS=nm_s[i];
		s2co.insert(pair<rslt_t,size_t>(nm_s[i],i));
	}
}

/* A job is the minimum unit of cross-validation analysis.
Suppose we have 5 genes (gene1..gene5), among which 2 are associated (1,4).
First, rank gene1 (target) among 1,2,3,5 by their closeness to gene4 (DzGene).
Then,  rank gene4 (target) among 2,3,4,5 by their closeness to gene1 (DzGene).
From each analysis, calculate AUC and store the result inside the JobData.
Repeat it for all associated genes, then calculate mean AUC for a trait.
Repeat it for all diseases, then calcualte the overall mean and SD of AUCs.

      1 :   +-------
        :   |
        :   |
Sens.   :   | AUCa = (1-a)*1
        :   |      = 1-a
        :   |
      0 +---+.......
        0   a      1
           1-spec.
Fig1: AUC for an associated gene of a disease.
Note: a is the proportion of unassociated genes whose rank is <= the target.

      1 :      +----
        :      |||||  <- (1-c)/3
     2/3:   +--+
Sens.   :   |\\\\\\\  <- (1-b)/3  AUC = (1-a)/3+(1-b)/3+(1-c)/3
	 1/3: +-+                         = (1-a + 1-b + 1-c)/3
        : |/////////  <- (1-a)/3      = mean(AUCa,AUCb,AUCc)
      0 +-+.........
        0 a b  c   1
           1-spec.
Fig2: AUC for all associated genes of a disease, eg, when there're 3 of them.
Note: a,b,c are the proportion of unassociated genes whose rank is <= the targets. 
However, this assumes that prediction scale doesn't change, which is not true.
The change is small if the #DiseaseGenes is large, and the score is standardized.
But still, it makes more sense to use gene-wise AUC instead of disease-wise AUC.*/

struct JobData {
	bool			nankg;			// input:  make s[i]=nan for known genes
	string			trait;			// input:  trait name
	vector<size_t>	coor_target;	// input:  coordinates of target genes. If empty, auc=nan.
	vector<size_t>	coor_KnGene;	// input:  coordinates of known  genes
	vector<double>	wght_KnGene;	// input:  func weight of known  genes
	rslt_t			auc;			// output: area under ROC curve
	rslt_t			tss;			// output: sum of squares
	JobData():nankg(true),auc(-1),tss(-1){}
};

struct AllJobsResultType {
	// data used by add()
	bool									wr_bf;	// write Bayes factor
	bool									wr_As;	// write prediction scores for all genes
	bool									mode2;	// whether it's MODE_2
	Interpolate								intpo;	// linear interpolation to get BF
	std::mutex								m;		// mutex for AllJobsResultType
	// data calcualted by summarize()
	std::map<std::string, Values<rslt_t> >	AUCs;	// AUCs[trait_name] = a Values of AUCs,
	std::map<std::string, boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::lazy_variance> > >	VARs;	// VARs[trait_name] = variances,
	std::map<std::string, boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::sum> > >			TSSs;	// TSSs[trait_name] = sum of squares,
	boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::lazy_variance> >	VARa;	// VARa = variances,
	boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::sum> >			TSSa;	// TSSa = sum of squares,
	
	AllJobsResultType():wr_bf(false),wr_As(false),mode2(false) {}
	
	template <typename PREDICTION_CONTAINER> // vector<rslt_t> S / Eigen::VectorXf FF
	void add(PREDICTION_CONTAINER& S, vector<int> labels, JobData& p)
	{
		m.lock();
		AUCs[p.trait].push_back(p.auc);
		
		if (mode2)
		{
			for (size_t i=0; i<num_genes; ++i) if (!std::isnan(S[i])) { VARa(S[i]); }
			TSSa(p.tss);
		}
		else
		{
			boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::lazy_variance> >& VARb = VARs[p.trait];
			boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::sum> >&			TSSb = TSSs[p.trait];
			for (size_t i=0; i<num_genes; ++i) if (!std::isnan(S[i])) { VARb(S[i]); }
			TSSb(p.tss);
		}

		if (wr_As)
		{
			for (size_t i=0; i<num_genes; ++i)
			{
				cout << S[i] << '\t' << labels[i] << '\t' << p.trait << '\t';
				print_container(p.coor_target,cout,',',true);
			}
		}
		if (wr_bf)
		{
			if (AUCs.size()==1 && AUCs.begin()->second.size()==1)
			{
				cout<<"Trait\\Gene"; for (size_t i=0; i<num_genes; ++i) cout << '\t' << coor2gene[i]; cout << endl;
			}
			cout << p.trait; for (size_t i=0; i<num_genes; ++i) cout << '\t' << log10(intpo.solve(S[i])); cout << endl;
		}
		m.unlock();
	}

	void summarize() // useful only for cross-validation
	{
		// AUCs by trait
		Values<rslt_t> summary;
		for (each_element(AUCs,it)) summary.push_back(it->second.get(STAT::MEAN));
		lns << showl << "AUC per disease (roughly): N="<<summary.get(STAT::N)<<" mean="<<summary.get(STAT::MEAN);
		if (summary.size()==1)	lns<<" SD not computable because there is only 1 trait." << flush_logger;
		else					lns<<" SD="<<summary.get(STAT::SPL_SD) << flush_logger;
		
		// AUCs by gene
		summary.clear();
		for (each_element(AUCs,it)) summary.push_back(it->second.begin(),it->second.end());
		lns << showl << "AUC per round of analysis: N="<<summary.get(STAT::N)<<" mean="<<summary.get(STAT::MEAN);
		if (summary.size()==1)	lns<<" SD not computable because there is only 1 gene." << flush_logger;
		else					lns<<" SD="<<summary.get(STAT::SPL_SD) << flush_logger;

		// var analysis
		if (mode2) // % explained by traits (valid only for MODE_2, ie each trait is run only once)
		{
			int		 N = boost::accumulators::count(VARa);
			double var = boost::accumulators::variance(VARa);
			double tss = var * N; // var is population variance
			double l1o = boost::accumulators::sum(TSSa);
			double explained = (1-l1o/tss);
			lns << showl << "Variance explained by diseases: " << explained << flush_logger;
		}
		else // % explained by cross-validation
		{
			summary.clear();
			for (auto &t:VARs)
			{
				int		 N = boost::accumulators::count(t.second);
				double var = boost::accumulators::variance(t.second);
				double tss = var * N; // var is population variance
				double l1o = boost::accumulators::sum(TSSs[t.first]);
				summary.push_back(1-l1o/tss);
			}
			lns << showl << "Variance explained by cross-validation: " << summary.get(STAT::MEAN) << flush_logger;
			for (auto &v:summary) lns << writel << "#VarExpl\t" << v << flush_logger;
		}
		
		// write individual AUCs
		lns << writel << "Trait\tmeanAUC\tAUCs" << flush_logger;
		for (each_element(AUCs,it))
			lns << writel << it->first << '\t' << it->second.get(STAT::MEAN) << '\t' << str_of_container(it->second,',') << flush_logger;
	}

} all_results;

void GBA_NaiveBayes(JobData& p)
{
	// prepare
	vector<rslt_t> weights(num_genes,0);
	if (NB_r3)
	{
		rslt_t half = (maxS-minS) * BinWidth / 2;
		for (auto &j:p.coor_KnGene)
		{
			multimap<rslt_t,size_t>::iterator its = s2co.lower_bound(nm_s[j] - half);
			multimap<rslt_t,size_t>::iterator ite = s2co.upper_bound(nm_s[j] + half);
//			for (multimap<rslt_t,size_t>::iterator it=its; it!=ite; ++it) weights[it->second] += 1; // uniform distribution
			for (multimap<rslt_t,size_t>::iterator it=its; it!=ite; ++it) weights[it->second] += 1 - std::abs(nm_s[j]-it->first)/half; // triangle distribution
		}
		MakeFraction(weights);
	}
	
	// calculate
	// method 1: w/o weighting
	// rslt_t factor1 = 1.0/p.coor_KnGene.size() ;
	// method 2: w/  weighting
	rslt_t factor1=0; for (auto &x:p.wght_KnGene) factor1+=x; factor1=1/factor1;
	
	vector<rslt_t> S(num_genes,0);
	for (size_t i=0;i<num_genes;++i)
	{
		// method 1: w/o weighting
		// for (auto &j:p.coor_KnGene) S[i] += W(i,j); // must be i,j for _RC matrix
		// method 2: w/  weighting
		for (size_t j=0; j<p.coor_KnGene.size(); ++j) S[i] += W(i,p.coor_KnGene[j])*p.wght_KnGene[j];
		
		if (NB_r1)
		{
			S[i] /= nm_s[i];
		}
		if (NB_r3)
		{
			rslt_t denominator = 0;
			for (size_t j=0;j<num_genes;++j) denominator += W(i,j)*weights[j];
			S[i] /= denominator;
		}
		S[i] *= factor1;
	}
	standardize(S);
	
	// summary
	vector<int>	labels(num_genes,0);
	if (p.nankg) for (each_element(p.coor_KnGene, j)) S[*j]=std::numeric_limits<rslt_t>::quiet_NaN();
	for (each_element(p.coor_target, j)) labels[*j]=1;
	p.auc = AUC_lite(labels,S);
	p.tss = tss(S);
	all_results.add(S,labels,p);
}

void GBA_GaussianSmoothing(JobData& p)
{
	cell_t pr_d = (cell_t)p.coor_KnGene.size() / num_genes;
	EVec_t ff ( EVec_t::Zero(num_genes) );			// f-final
	EVec_t f0 ( EVec_t::Constant(num_genes,pr_d) );	// f-original
	for (each_element(p.coor_KnGene,it)) f0(*it)=1;
	ff = A * f0;
	array_standardize(ff);
	
	vector<int>	labels(num_genes,0);
	if (p.nankg) for (each_element(p.coor_KnGene, j)) ff(*j)=std::numeric_limits<rslt_t>::quiet_NaN();
	for (each_element(p.coor_target, j)) labels[*j]=1;
	p.auc = AUC_lite(labels,ff);
	p.tss = array_tss(ff);
	all_results.add(ff,labels,p);
}

void OtherNetworkMeasures(JobData& p)
{
	vector<rslt_t> nm_c(num_genes,0); // weighted clustering coefficient = sum(Wij+Wih)/Si/(Ki-1)/2
	vector<rslt_t> nm_n(num_genes,0); // weighted average nearest-neighbors degree = sum(Wij*Kj)/Si

	// calculate (slow)
	for (size_t i=0;i<num_genes;++i)
		for (size_t j=0;j<num_genes;++j)
			for (size_t h=j+1;h<num_genes;++h)
				if (i!=j && i!=h && W(j,i) && W(h,i) && W(j,h)) nm_c[i] += ( W(j,i)+W(h,i) );
	for (size_t i=0;i<num_genes;++i)
		if (nm_c[i]) { rslt_t d = nm_s[i] * (nm_k[i]-1); nm_c[i]/=d; }
	
	for (size_t i=0;i<num_genes;++i)
		for (size_t j=0;j<num_genes;++j)
			if (i!=j && W(j,i)) { nm_n[i] += W(j,i) * nm_k[j]; }
	for (size_t i=0;i<num_genes;++i)
		if (nm_s[i]) nm_n[i] /= nm_s[i];

	// output all calculated network measures
	for (size_t i=0;i<num_genes;++i)
		cout << coor2gene[i] << '\t' << nm_k[i] << '\t' << nm_s[i] << '\t' << nm_c[i] << '\t' << nm_n[i] << '\n';
}

void calculate_all(string DzName, vector<size_t>& DzGenes, vector<size_t>& SdGenes, vector<double>& WtGenes, MtJobs<JobData>& CompJobs)
{
	if (DzName.empty()) return;
	if (DzGenes.size()<MinSeeds) return;
	CompJobs.push_back(JobData());
	JobData& prj = CompJobs.back();
	prj.trait = DzName;
	prj.nankg = false;
	prj.wght_KnGene = WtGenes;
	prj.coor_KnGene = SdGenes;
	prj.coor_target = DzGenes;
}

void leave_one_out(string DzName, vector<size_t>& DzGenes, vector<size_t>& SdGenes, vector<double>& WtGenes, MtJobs<JobData>& CompJobs)
{
	if (DzName.empty()) return;
	if (DzGenes.size()<MinSeeds) return;
	for (size_t i=0; i<DzGenes.size(); ++i)
	{
		CompJobs.push_back(JobData());
		JobData& prj = CompJobs.back();
		prj.trait = DzName;
		prj.nankg = true;
		for (size_t j=0; j<DzGenes.size(); ++j)
		{
			if (i==j)	prj.coor_target.push_back(DzGenes[j]);
			else	{	prj.coor_KnGene.push_back(SdGenes[j]); prj.wght_KnGene.push_back(WtGenes[j]); }
		}
	}
}

void cross_validat(string DzName, vector<size_t>& DzGenes, vector<size_t>& SdGenes, vector<double>& WtGenes, MtJobs<JobData>& CompJobs)
{
	if (DzName.empty()) return;
	if (DzGenes.size()<MinSeeds) return;
	vector<size_t>	IDs(DzGenes.size());
	for (size_t i=0; i<IDs.size(); ++i) IDs[i]=i%NumFolds;
	std::random_shuffle ( IDs.begin(), IDs.end(), MyRandom);
	for (size_t i=0; i<NumFolds; ++i)
	{
		CompJobs.push_back(JobData());
		JobData& prj = CompJobs.back();
		prj.trait = DzName;
		prj.nankg = true;
		for (size_t j=0; j<DzGenes.size(); ++j)
			if (IDs[j]==i)	prj.coor_target.push_back(DzGenes[j]);
			else		{	prj.coor_KnGene.push_back(SdGenes[j]); prj.wght_KnGene.push_back(WtGenes[j]); }
	}
}

void (*GBA)(JobData&) = NULL;
void (*JOB)(string, vector<size_t>&, vector<size_t>&, vector<double>& WtGenes, MtJobs<JobData>& ) = &calculate_all;

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);

	// parameters
	std::string genefile="g_symbol"; // gene list for all human genes
	std::string mtrxfile="g_mtx_lo"; // matrix file
	std::string dzgbfile; // disease gene bundled file - 1 file for all diseases
	std::string GtoBfile; // GBA score -> Bayes factor translation file
	std::string RbyGfile; // Random seeds by Group,   one group of GeneIDs per line, no header row/column
	std::string PbyGfile; // Permute GeneIDs by Group,one group of GeneIDs per line, no header row/column
	bool					l1oForKg=false;
	boost::logic::tribool	MatrixIC(boost::indeterminate);
	boost::logic::tribool	MatrixRC(boost::indeterminate);
	
	// handle program options
#ifdef PERCH_DISTRIBUTION
#else
	program.enable_option("--nt");
#endif
	program.forbid_option("-o");
	program.enable_option("--prefix");
	program.read_arguments(argc,argv,true,true);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"-g"))		ReadArg(program.arg(),argi,genefile);
		else if	(str_startsw(program.arg()[argi],"-m"))		ReadArg(program.arg(),argi,mtrxfile);
		else if	(str_startsw(program.arg()[argi],"-s"))		ReadArg(program.arg(),argi,dzgbfile);
		else if	(str_startsw(program.arg()[argi],"--seeds"))ReadArg(program.arg(),argi,dzgbfile);
		else if (str_startsw(program.arg()[argi],"-b"))		ReadArg(program.arg(),argi,GtoBfile);
		else if (str_startsw(program.arg()[argi],"-r"))		ReadArg(program.arg(),argi,RbyGfile);
		else if (str_startsw(program.arg()[argi],"-p"))		ReadArg(program.arg(),argi,PbyGfile);
		else if (str_startsw(program.arg()[argi],"-sg"))	ReadArg(program.arg(),argi,MinSeeds);
		else if (str_startsw(program.arg()[argi],"-cv")) {	ReadArg(program.arg(),argi,NumFolds); JOB = &cross_validat; }
		else if (program.arg()[argi]=="-l1o")										JOB = &leave_one_out;
		else if (program.arg()[argi]=="--ic"){	GBA = &GBA_NaiveBayes;			MatrixIC=true;  MatrixRC=false; NB_r1=true;  NB_r3=false; all_results.wr_bf=true; }
		else if (program.arg()[argi]=="--gs"){	GBA = &GBA_GaussianSmoothing;	MatrixIC=false; MatrixRC=false; NB_r1=false; NB_r3=false; all_results.wr_bf=true; }
		else if (program.arg()[argi]=="-NM") {	GBA = &OtherNetworkMeasures;					MatrixRC=false; NB_r1=false; NB_r3=false; }
		else if (program.arg()[argi]=="-GS") {	GBA = &GBA_GaussianSmoothing;	MatrixIC=false; MatrixRC=false; NB_r1=false; NB_r3=false; }
		else if (program.arg()[argi]=="-NB") {	GBA = &GBA_NaiveBayes;			MatrixIC=false; MatrixRC=false; NB_r1=false; NB_r3=false; }
		else if (program.arg()[argi]=="-N1") {	GBA = &GBA_NaiveBayes;			MatrixIC=false; MatrixRC=false; NB_r1=true;  NB_r3=false; }
		else if (program.arg()[argi]=="-N2") {	GBA = &GBA_NaiveBayes;			MatrixIC=false; MatrixRC=true;  NB_r1=false; NB_r3=false; }
		else if (program.arg()[argi]=="-N3") {	GBA = &GBA_NaiveBayes;			MatrixIC=false; MatrixRC=false; NB_r1=false; NB_r3=true;  }
		else if (program.arg()[argi]=="-IC") {	GBA = &GBA_NaiveBayes;			MatrixIC=true;  MatrixRC=false; NB_r1=false; NB_r3=false; }
		else if (program.arg()[argi]=="-I1") {	GBA = &GBA_NaiveBayes;			MatrixIC=true;  MatrixRC=false; NB_r1=true;  NB_r3=false; }
		else if (program.arg()[argi]=="-I2") {	GBA = &GBA_NaiveBayes;			MatrixIC=true;  MatrixRC=true;  NB_r1=false; NB_r3=false; }
		else if (program.arg()[argi]=="-I3") {	GBA = &GBA_NaiveBayes;			MatrixIC=true;  MatrixRC=false; NB_r1=false; NB_r3=true;  }
		else if (program.arg()[argi]=="--wr-a") all_results.wr_As=true;
		else if (program.arg()[argi]=="--wr-b") all_results.wr_bf=true;
		else if (program.arg()[argi]=="-l1o-sg")l1oForKg=true;
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	perch::check_arguments();

	if (genefile.empty()) exit_error("option -g is required.");
	if (mtrxfile.empty()) exit_error("Option -m is required.");
	if (dzgbfile.empty()) exit_error("Option -s is required.");
#ifndef PERCH_DISTRIBUTION
	if (!MatrixIC &&  str_has(mtrxfile,"_IC")) exit_error("You are not supposed to use the _IC matrix.");
	if (!MatrixIC &&  str_has(mtrxfile,"_RC")) exit_error("You are not supposed to use the _RC matrix.");
	if (!MatrixRC && !str_has(mtrxfile,"_lo")) mtrxfile+="_lo";
	if ( MatrixIC && !str_has(mtrxfile,"_IC")) mtrxfile+="_IC";
	if ( MatrixRC && !str_has(mtrxfile,"_RC")) mtrxfile+="_RC";
#endif
	if (JOB==&leave_one_out||JOB==&cross_validat) // MODE_1
	{
		if (GBA==&OtherNetworkMeasures) exit_error("Option -NM is for MODE_2 only.");
		if (GBA==&GBA_NaiveBayes) l1oForKg=true; // -l1o-sg is for MODE_2 only; in MODE_1 it's always true
	}
	else // MODE_2
	{
		all_results.mode2=true;
	}
	if (!GtoBfile.empty()) all_results.wr_bf=true;
	if (all_results.wr_bf && all_results.wr_As) exit_error("--wr-b and --wr-a cannot be both true.");
	if (GBA!=&GBA_NaiveBayes && l1oForKg) exit_error("-l1o-sg is for -NB/-I1/-I2 only.");

	// for distribution (no GBA options / -cv / -l1o / --wr-x / -l1o-sg etc), default -I1
	if (GBA==NULL)
	{
		// GBA=&GBA_NaiveBayes; MatrixIC=true;  MatrixRC=false; NB_r1=true;  NB_r3=false; all_results.wr_bf=true; // old 
		GBA=&GBA_GaussianSmoothing;	MatrixIC=false; MatrixRC=false; NB_r1=false; NB_r3=false; all_results.wr_bf=true;
	}
	
	// log program options
#ifndef PERCH_DISTRIBUTION
	lns.open(program.prefix()+".log", ios::out);
#endif
	if (JOB==&cross_validat && MinSeeds<NumFolds)
	{	MinSeeds=NumFolds;
		lns<<showl<<"-sg was automatically increased to "<<NumFolds<<" for -cv"<<NumFolds<<" to work."<<flush_logger; }
	if		(GBA==&GBA_GaussianSmoothing)		lns<<showl<<"Calculate GBA by Gaussian Smoothing ";
	else if (GBA==&GBA_NaiveBayes &&  MatrixIC)	lns<<showl<<"Calculate GBA by InterConnectedness ";
	else if (GBA==&GBA_NaiveBayes && !MatrixIC)	lns<<showl<<"Calculate GBA by Naive Bayes label propagation ";
	else if (GBA==&OtherNetworkMeasures)		lns<<showl<<"Calculate network measures ";
	if (NB_r1)		lns<<"{calibrated}"; // for both NB and IC
	if (MatrixRC)	lns<<"[calibrated]"; // for both NB and IC
	if (NB_r3)		lns<<"(calibrated)"; // for both NB and IC
	lns << flush_logger;
	lns<<showl<<"Number of threads = "<<program.nt << flush_logger;
	if (all_results.wr_As) lns<<showl<<"output 'GBA_Score Label TraitName TargetGenes' for all genes" << flush_logger;
	if (all_results.wr_bf) lns<<showl<<"output log10(LR) table (row=traits; col=genes; both have a header)" << flush_logger;
	if (GBA==&OtherNetworkMeasures) lns<<showl<<"Output 'GeneID Degree Strength Clustering NearestNeigborsDegree' for all genes" << flush_logger;
	
	// read Bayes factor curve
	if (!GtoBfile.empty())
	{
		lns<<showl<<"Use "<<GtoBfile<<" to calculate LR from GBA." << flush_logger;
		all_results.intpo.setup(GtoBfile);
	}
	else if (all_results.wr_bf)
	{
		if		(GBA==&GBA_NaiveBayes && MatrixIC==true && MatrixRC==false && NB_r1==true && NB_r3==false)	{ stringstream ss(default_I2_to_BF); all_results.intpo.setup(ss); }
		else if	(GBA==&GBA_GaussianSmoothing)																{ stringstream ss(default_GS_to_BF); all_results.intpo.setup(ss); }
		else exit_error("No reference data for LR calculation from GBA. Please use -b or -I1.");
		lns<<showl<<"Use the default reference to calculate LR from GBA." << flush_logger;
	}
	
	// read gene names
	lns<<showl<<"Read "<<genefile<<" ... "<<flush_logger;
	for (Rows_in_File(in,perch::find_file(genefile),1))
	{
		if (!exist_element(gene2coor, in[0])) { size_t c=gene2coor.size(); gene2coor[in[0]]=c; coor2gene.push_back(in[0]); }
		else exit_error("Duplicated gene "+in[0]+" in "+genefile);
	}
	num_genes = gene2coor.size();
	W = EMtx_t::Zero(num_genes, num_genes);
	lns<<showl<<"  "<<coor2gene.size()<<" genes read." << flush_logger;
	
	// read gene group for gene permutation
	if (!PbyGfile.empty())
	{
		lns<<showl<<"Read "<<PbyGfile<<" for stratified matrix permutation ... " << flush_logger;
		tfile_format fmt;
		fmt.forbid_nf_rpt();
		for (Rows_in_File(in,PbyGfile,&fmt))
		{
			geneGroup.push_back(set<size_t>());
			for (auto &i : in.contents()) geneGroup.back().insert(gene2coor[i]);
		}
		lns<<showl<<"  "<<geneGroup.size()<<" groups read." << flush_logger;
		
		std::vector<std::string> c2gPermut(num_genes);
		for (auto &gr:geneGroup)
		{
			vector<size_t> original; for (auto &c:gr) original.push_back(c);
			vector<size_t> shuffled = original;
			std::random_shuffle ( shuffled.begin(), shuffled.end(), MyRandom );
			for (size_t i=0;i<original.size();++i) c2gPermut[original[i]]=coor2gene[shuffled[i]];
		}
		coor2gene=c2gPermut;
		gene2coor.clear(); for (size_t i=0;i<num_genes;++i) gene2coor[coor2gene[i]]=i;
		geneGroup.clear();
	}
	
	// read gene group for random seed picking
	if (RbyGfile=="all")
	{
		geneGroup.push_back(set<size_t>());
		for (size_t i=0;i<gene2coor.size();++i) geneGroup.back().insert(i);
	}
	else if (!RbyGfile.empty())
	{
		lns<<showl<<"Read "<<RbyGfile<<" for random pseudo-null seed-picking ... " << flush_logger;
		tfile_format fmt;
		fmt.forbid_nf_rpt();
		for (Rows_in_File(in,RbyGfile,&fmt))
		{
			geneGroup.push_back(set<size_t>());
			for (auto &i : in.contents()) geneGroup.back().insert(gene2coor[i]);
		}
		lns<<showl<<"  "<<geneGroup.size()<<" groups read.";
	}

	// read matrix and prepare matrix
	lns<<showl<<"Read "<<mtrxfile<<" ... "<<flush_logger;
	cell_t maxcell = 0;
	{
		if		(str_has(mtrxfile,"_up"))
		{
			exit_error("Upper triangle is not supported due to the confusing formats.");
		}
		else if (str_has(mtrxfile,"_lo"))
		{
			tfile_format fmt;
			fmt.forbid_nf_rpt();
			for (Rows_in_File(in, perch::find_file(mtrxfile), &fmt))
			{
				int row = in.RowNumber()+1;
				for (int col=0;col<row;++col)
				{
					cell_t r = 0;
					try { r=boost::lexical_cast<cell_t>(in[col]); }
					catch (...) { exit_error("Error reading "+mtrxfile+": failed to read "+in[col]+" as a number."); }
					if (r>maxcell) maxcell=r;
					W(row,col)=r;
					W(col,row)=r;
				}
			}
		}
		else
		{
			tfile_format fmt;
			fmt.set_field_nums(num_genes, "lines do not have enough fields.",tfile_format::Continue);
			for (Rows_in_File(in, perch::find_file(mtrxfile), &fmt))
			{
				size_t row = in.RowNumber();
				for (size_t col=0;col<num_genes;++col)
				{
					cell_t r = 0;
					try { r=boost::lexical_cast<cell_t>(in[col]); }
					catch (...) { exit_error("Error reading "+mtrxfile+": failed to read "+in[col]+" as a number."); }
					if (row==col && r!=0) exit_error("The matrix' diagonal should be 0.");
					if (r>maxcell) maxcell=r;
					W(row,col)=r;
				}
			}
		}
	}
	lns<<showl<<"  done"<<flush_logger;

	if (GBA==&GBA_GaussianSmoothing)
	{
		lns<<showl<<"Prepare matrix ... "<<flush_logger;
		EMtx_t D (EMtx_t::Zero(num_genes, num_genes));
		EMtx_t I (EMtx_t::Identity(num_genes, num_genes));
		for (size_t i=0;i<num_genes;++i) D(i,i)=W.row(i).sum();
		// method 1 (00.629479)
		A = EMtx_t::Zero(num_genes, num_genes);
		A = (I+(D-W)).inverse(); //*/
		/*/ method 2 (even slower than method 1 -- 01.303636)
		// requires #include <Eigen/LU>	// for Eigen::FullPivLU
		EMtx_t B (EMtx_t::Zero(num_genes, num_genes));
		B = I+(D-W);
		Eigen::FullPivLU<EMtx_t> lu(B);
		A = lu.inverse(); //*/
		lns<<showl<<"  done"<<flush_logger;
	}

	lns<<showl<<"Calculate basic network measures ... "<<flush_logger;
	BasicNetworkMeasures();
	lns<<showl<<"  done"<<flush_logger;
	
	if (GBA==&GBA_NaiveBayes && !l1oForKg) for (size_t i=0;i<num_genes;++i) W(i,i)=maxcell;
	// diagonal=max => known genes have high BF but others are unaffected;
	// without it, BFs for known genes are equal to leave-one-out analysis.
	// Think what should it be for calibrated matrix.
	
	// job data
	MtJobs<JobData> CompJobs;
	struct DzGeneDataType
	{
		string				last_dz;
		map<size_t,double>	DzGuniq;
		void add_jobs_to(MtJobs<JobData>& CJ)
		{
			vector<size_t> DzGenes;
			vector<double> WtGenes;
			vector<size_t> SdGenes;
			for (auto &x:DzGuniq) { DzGenes.push_back(x.first); WtGenes.push_back(x.second); }
			bool failed=false;
			if		(geneGroup.size()==1)
			{
				set<size_t> RndSeed;
				while (RndSeed.size() < DzGuniq.size())
				{
					size_t j = MyRandom(num_genes);
					if (!exist_element(DzGuniq,j)&&!exist_element(RndSeed,j)) RndSeed.insert(j);
				}
				for (auto &x:RndSeed) SdGenes.push_back(x);
			}
			else if (geneGroup.size()>1)
			{
				vector< set<size_t> > candidates = geneGroup;
				for (size_t i=0; i<candidates.size(); ++i)
				{
					size_t need = 0;
					for (set<size_t>::iterator iter=candidates[i].begin(); iter!=candidates[i].end(); )
						if (exist_element(DzGuniq,*iter)) { candidates[i].erase(iter++); ++need; }
						else ++iter;
					if (need > candidates[i].size()) { failed=true; break; }
					for (;need;--need) {
						set<size_t>::iterator iter = candidates[i].begin();
						std::advance( iter, MyRandom(candidates[i].size()) );
						SdGenes.push_back(*iter);
						candidates[i].erase(iter);
					}
				}
			}
			else
			{
				SdGenes=DzGenes;
			}
			if (!failed) (*JOB)(last_dz, DzGenes, SdGenes, WtGenes, CJ);
			DzGuniq.clear();
		}
		void add_data(const string& trait, size_t gene_coor, double weight, MtJobs<JobData>& CJ)
		{
			if (trait!=last_dz) add_jobs_to(CJ);
			last_dz=trait;
			DzGuniq[gene_coor]=weight;
		}
	} DzGeneData;
	
	// read trait file and create jobs
	lns<<showl<<"Read "<<dzgbfile<<" ... "<<flush_logger;
	int num_fields=0;
	for (Rows_in_File(in,dzgbfile,2)) // 2 col: dz,gene. Lines are sorted by dz.
	{
		// decide file format
		if (num_fields)
		{
			if (num_fields!=in.NumFields()) exit_error("seed gene file has different number of columns");
		}
		else
		{
			num_fields=in.NumFields();
			if (num_fields!=2 && num_fields!=3) exit_error("seed gene file should have either 2 or 3 columns");
		}
		
		// read data
		if		(num_fields==2)
		{
			if (exist_element(gene2coor,in[1]))
				DzGeneData.add_data(in[0],gene2coor[in[1]],1,CompJobs);
		}
		else if (num_fields==3)
		{
			double w;
			if (!read_val_gt_le(in[2],w,0.0,1.0)) exit_error("seed gene weight should be in (0,1]");
			if (exist_element(gene2coor,in[1]))
				DzGeneData.add_data(in[0],gene2coor[in[1]],w,CompJobs);
		}
		else exit_error("seed gene file should have either 2 or 3 columns");
	}
	DzGeneData.add_jobs_to(CompJobs);
	if (CompJobs.empty()) exit_error("There is no job.");
	lns<<showl<<"  "<<CompJobs.size()<<" jobs created." << flush_logger;
	
	// do jobs
	lns<<showl<<"Calculating ... "<<flush_logger;
	if (GBA==&OtherNetworkMeasures) cout << "#Degree\tStrength\tClusteringCoefficient\tAverageNearestNeighborsDegree" << endl;
	if (program.nt==1)	for (each_element(CompJobs, prj)) (*GBA)(*prj);
	else				CompJobs.run(*GBA, program.nt);
	lns<<showl<<"  done" << flush_logger;
	
	// summarize
	if (JOB != &calculate_all) all_results.summarize();
	
	return 0;
}
