// body file: vAnnGene_anc.cpp

#ifndef PERCH_VAG_ANC
#define PERCH_VAG_ANC

#include <map>
#include <set>
#include <vector>

extern std::map<std::string, std::map< std::pair<int,int>,std::pair<std::string,std::string> > >	ancestral_var_r; // ancestral var tx (HGVS r.)
extern std::set<std::string>																		ancestral_var_g; // ancestral var g  (HGVS g.)
void read_ancestral_var(const std::string& input);

#endif
