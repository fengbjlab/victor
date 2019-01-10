// body file: libfbj_fet.cpp

#ifndef RC_FISHER_EXACT_TEST
#define RC_FISHER_EXACT_TEST

#define FET_MAX_C 100
#define FET_MAX_R 100

#include <iostream>
#include <cmath>
#include <cstring> // to compile memset memcpy memcmp by "gcc/4.7.2/bin/g++ -std=c++11" but not "gcc/4.1.2/bin/g++"
#include "libfbj_factorial.hpp"
#include "libfbj_base.hpp"

/* 
STEPS:
1) FisherExactTest fet;
2) fet.set_empty_table(r,c);
   fet.input(r,c,count);
2) fet.input_RxC(table,r,c);
3) p=fet.result();
NOTES: 
1) Data won't change after result(), so it's safe to input() again without set_empty_table().
2) If data are non-negative, it still runs and result will be wrong. So be careful.
 
 To compare
   this:  contingency-table 8 3 102 0 --fet
   STATA: csi 8 102 3 0 ,  exact
*/

class FisherExactTest {
private:
	static int lc;
	int FET_r;							//row
	int FET_c;							//column
	int FET_a[FET_MAX_R][FET_MAX_C];	//array
	int FET_rm[FET_MAX_R];				//row marginal total
	int FET_cm[FET_MAX_C];				//column marginal total
	int FET_gt;							//grand total
	long double FET_pv_this_tail;		//p value
	long double FET_pv_othr_tail;		//p value
	long double FET_pv_both_tail;		//p value
	bool modified;						//whether table modified after calculation
	bool mid_p;							//whether to perform mid p-value test
	void minimize(int fix_r,int fix_c,int next);
	int increment();
	FACTORIAL pvalue(FACTORIAL common_factor);
	void test();

public:
	
	FisherExactTest(bool mp):FET_r(0),FET_c(0),FET_pv_this_tail(-9999),FET_pv_othr_tail(-9999),FET_pv_both_tail(-9999),modified(false) { mid_p=mp; }
	
	inline void set_empty_table(int r,int c) {
		if (r>FET_MAX_R) { exit_error("FET:Max number of rows is "+itos(FET_MAX_R)); }
		if (c>FET_MAX_C) { exit_error("FET:Max number of columns is "+itos(FET_MAX_C)); }
		if (r<=0 || c<=0){ exit_error("FET:Number of rows/columns must be >=1."); }
		memset(FET_a, 0, sizeof(int)*FET_MAX_R*FET_MAX_C);
		FET_r=r;
		FET_c=c;
		modified=true;
	}
	
	inline void input(int r, int c,int count)
	{
		FET_a[r][c]=count;
		modified=true;
	}
	
	template <typename T>
	inline void input_RxC(T table[], int nr, int nc)
	{
		for (int ij=0,r=0;r<nr;++r)
			for (int  c=0;c<nc;++c,++ij)
				FET_a[r][c]=table[ij];
		FET_r=nr;
		FET_c=nc;
		modified=true;
	}
	
	inline long double result(const char& side='2') // the p value is two-tailed by default
	{
		if (FET_r<=0 || FET_c<=0) { exit_error("FET:Number of rows/columns <=0. Please call set_empty_table()."); }
		if (modified) test();
		if		(side=='2')	return FET_pv_both_tail;
		else if (side=='L')	return FET_pv_this_tail;
		else if (side=='R')	return FET_pv_othr_tail;
		else if (side=='l')	return FET_pv_this_tail;
		else if (side=='r')	return FET_pv_othr_tail;
		else if (side=='1')	return std::min(FET_pv_this_tail,FET_pv_othr_tail);
		else { exit_error("The side for FET:result(side) must be 1/2/L/R."); }
		return std::numeric_limits<double>::signaling_NaN(); // won't happen. Just to shut up the compiler for warnings.
	}
};

template <typename T>
inline double Fishers_exact_test_2x2(T table[], bool mp)
{
	FisherExactTest fet(mp);
	fet.input_RxC(table,2,2);
	return fet.result();
}

template <typename T>
inline double Fishers_exact_test_RxC(T table[], int nr, int nc, bool mp)
{
	FisherExactTest fet(mp);
	fet.input_RxC(table,nr,nc);
	return fet.result();
}

/*/ debug FET
double p;

int tb1[2][2]={ 3,7,15,18 };
p=Fishers_exact_test_kx2(tb1,2); // or Fishers_exact_test_2x2(tbl,p);
cout<<"p="<<p<<std::endl; // p=0.479932

int tb2[2][2]={ 3,7,18,15 };
p=Fishers_exact_test_kx2(tb2,2);
cout<<"p="<<p<<std::endl; // p=0.280579

int tb3[2][2]={ 0,0,7556,3480 };
p=Fishers_exact_test_kx2(tb3,2);
cout<<"p="<<p<<std::endl; // p=1

int tb4[2][2]={ 0,12,10,0 };
p=Fishers_exact_test_kx2(tb4,2);
cout<<"p="<<p<<std::endl; // p=1.54644e-06, some website output 0 but that's not correct
// it depends on whether test() should be p<=FET_th or p<FET_th. Correct one is <=.

int tb5[2][2]={ 0,7556,3480,0 };
p=Fishers_exact_test_kx2(tb5,2);
cout<<"p="<<p<<std::endl; // p=inf and is very slow. When FACTORIAL use direct calculation (no expansion) p=0, which is more reasonable.
//*/

#endif

