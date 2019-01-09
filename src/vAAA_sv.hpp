// body file: vAAA_sv.cpp

#ifndef PERCH_AAA_SV
#define PERCH_AAA_SV
#include <vector>

namespace AAA_SV
{	
	void sv_read_wiSymb(const std::string& filename); 				// read output of "vSVA --detail filename"
	void sv_read_wiSymb(const std::vector<std::string>& filenames); // read output of "vSVA --detail filename"
	void sv_read_woSymb(const std::string& filename); 				// read output of "vSVA --detail filename"
	void sv_read_woSymb(const std::vector<std::string>& filenames); // read output of "vSVA --detail filename"
	int sv_copy(const int& chr_num, const std::string& symbol, const std::string& SeqID, const int& bp);
	int sv_copy(const int& chr_num, const std::string& symbol, const std::string& SeqID);
	int sv_exist(const int& chr_num, const std::string& symbol, const std::vector<std::string>& SeqIDs);
};

#endif
