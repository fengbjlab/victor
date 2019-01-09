// body file: vAnnGene_knClinSig.cpp

#ifndef PERCH_VAG_VKS
#define PERCH_VAG_VKS

#include <set>
#include <map>

int GranthamScore(char AA1, char AA2);
int Blosum62Score(char AA1, char AA2);

class known_ClinSig_var {
private:
	std::set<std::string>				rmFH;				// remove full HGVS from HGVS_id and genomic_id
	std::set<std::string>				rmPH;				// remove partial HGVS from pHGVS_BayesDel, pHGVS_Grantham, pHGVS_Blosum62
	std::map<std::string,std::string>	genomic_id[2];		// left-normalized chr_pos_ref_alt
	std::set<std::string>				HGVS_id[2];			// full HGVS nomenclature
	std::map<std::string,double>		pHGVS_BayesDel[2];	// nm123:X456 -> minDel for path, maxDel for benign
	std::map<std::string,double>		pHGVS_Grantham[2];	// nm123:X456 -> minDel for path, maxDel for benign
	std::map<std::string,double>		pHGVS_Blosum62[2];	// nm123:X456 -> maxDel for path, minDel for benign
public:
	void		read(const std::string& filename);	// 9 columns: #CHROM,POS,REF,ALT,Symbol,Type,HGVS,DelScore,ClinSig
	std::string	test(const std::string& id, const std::string& FuncType, const std::string& hgvs, const double& del);
	// return [01].[abcd]. a:genomic_match b:protein_match c:same_aa_stronger_ev d:same_aa_weaker_ev
};

#endif
