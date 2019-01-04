// body file: vAnnGene_mes5.cpp

#ifndef PERCH_VAG_MES
#define PERCH_VAG_MES

#include <map>
#include <vector>

class MaxEntScan_score5ss {
private:
	std::vector<double>	scores;
public:
	void	read(const std::string& filename);	// 2 columns: seq score
	double	score(const std::string& seq);
	double	max(const std::string& seq);
};

#endif
