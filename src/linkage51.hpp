// body file: linkage51_makeped.cpp
// body file: linkage51_unknown.cpp
// body file: linkage51_mlink.cpp

#ifndef LINKAGE51_LIB
#define LINKAGE51_LIB

#include <iostream>

namespace linkage
{
	enum MlinkResultType { LOD, HLOD, PLOD, LOG10L };
	extern MlinkResultType defaultResultType;
}

void makeped_program(std::istream& pedfile, std::ostream& pedout,int no_ques);
void unknown_program(std::istream& datafile, std::istream& pedfile, std::ostream& ipedfile, std::ostream& speedfile, int& errors);
void mlink_calculate(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& theta, double& alpha, double& result);
int makeped_main(int argc, const char * argv[]);
int unknown_main();
int mlink_main();

#ifdef LINKAGE_OPTIMIZE
void mlink_opt_theta(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& theta, double& result);
void mlink_opt_scale(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& scale, double& result);
void mlink_opt_domin(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& phenocopy, double& penetrance, double& result);
#endif

#endif
