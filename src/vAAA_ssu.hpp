// body file: vAAA_ssu.cpp

#ifndef LIBFBJ_SumSqU_TESTS
#define LIBFBJ_SumSqU_TESTS

#include <Eigen/Core>
#include <Eigen/Dense>
#define EIGEN_NO_DEBUG	// disable Eigen range checking for faster access

int ssu_validate(const Eigen::MatrixXd & X,
				 const Eigen::VectorXd & y);

// ----------------- SumSqUs test -----------------

struct SumSqUs_workspace {
	int ni;
	int np;
	double			Ybar;
	Eigen::VectorXd Xbar;
	Eigen::VectorXd Y_Ybar;
	Eigen::MatrixXd X_Xbar;
	Eigen::VectorXd U;
	Eigen::MatrixXd CovS;
	Eigen::VectorXd cr;
	SumSqUs_workspace(int NumObs, int NumPar) {
		ni=NumObs;
		np=NumPar;
		Xbar	. setZero(np);
		Y_Ybar	. setZero(ni);
		X_Xbar	. setZero(ni,np);
		U		. setZero(np);
		CovS	. setZero(np,np);
		cr		. setZero(np);
	}
};

double SumSqUs_test(const Eigen::MatrixXd &X, const Eigen::VectorXd & Y, SumSqUs_workspace & ws);

// ----------------- SumSqUw test -----------------

// Only the variables independent of Y is stored here for thread-safety.
// Since X is stable, so they remain unchanged during permutation of Y.
// Sometimes you permute one of the Xs in permutation, when there're covariates.
// It's not efficient to make a Ystable version because many things depend on X.
// In that case you cannot use the last function.

struct SumSqUw_Xstable {
	int ni;
	int np;
	double			Ybar;
	Eigen::VectorXd Xbar;
	Eigen::MatrixXd Xtrn;
	Eigen::MatrixXd X_Xbar;
	Eigen::MatrixXd CovS;
	Eigen::MatrixXd diagCovS;
	SumSqUw_Xstable(int NumObs, int NumPar) {
		ni=NumObs;
		np=NumPar;
		Xbar	. setZero(np);
		Xtrn	. setZero(np,ni);
		X_Xbar	. setZero(ni,np);
		CovS    . setZero(np,np);
		diagCovS. setZero(np,np);
	}
};

// SSUw test, weight is the inverse of the sample SD of the variant
double SumSqUw_test(const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, SumSqUw_Xstable& ws);

// SSUw test with customized weights
double SumSqUw_test(const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, const std::vector<double> &W, SumSqUw_Xstable& ws);

// For permutation of Y only (X unchanged). Use the same WorkSpace.
double SumSqUw_test(const Eigen::VectorXd &Y, SumSqUw_Xstable& ws);

extern std::vector<double> pv_SSUw_default_empty_weights;
double pv_SSUw(const Eigen::MatrixXd & X, const Eigen::VectorXd	& y, const std::vector<double> &W=pv_SSUw_default_empty_weights);
double pv_SSUc(const Eigen::MatrixXd & X, const Eigen::VectorXd	& y, const std::vector<double> &W);

#endif
