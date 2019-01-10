// body file: libfbj_regress.cpp

#ifndef TFT_REGRESSION
#define TFT_REGRESSION

#include <Eigen/Core>
#include <Eigen/Dense>
#define EIGEN_NO_DEBUG	// disable Eigen range checking for faster access

struct multifit_logistic_workspace {
	int ni;
	int np;
	Eigen::VectorXd p ;		// prob of response
	Eigen::VectorXd V ;		// logit(p)
	Eigen::VectorXd t ;		// temp vector
	Eigen::MatrixXd X_;		// X~
	multifit_logistic_workspace(int NumObs, int NumPar) {
		ni=NumObs;
		np=NumPar;
		p .setZero(ni);
		V .setZero(ni);
		t .setZero(ni);
		X_.setZero(ni,np);
	}
};

// return 1 if good.
int multifit_logistic_validate (const Eigen::MatrixXd & X,
								const Eigen::VectorXd & y,
								const int p_start=1);	// if X[][0] is constant term, skip it by p_start=1

void multifit_logistic (const Eigen::MatrixXd & X,			// predictor matrix
						const Eigen::VectorXd & y,			// response vector
						Eigen::VectorXd & c,				// coefficients
						Eigen::MatrixXd & cov,				// covariance
						double & chisq,						// likelihood ratio test
						multifit_logistic_workspace& work);	// workspace

struct multifit_linear_workspace {
	int ni;
	int np;
	Eigen::MatrixXd t;		// temporary matrix (X'X)^-1X'
	Eigen::MatrixXd I;		// Indentity matrix
	Eigen::MatrixXd H;		// Hat matrix
	Eigen::VectorXd e;		// residuals
	multifit_linear_workspace(int NumObs, int NumPar) {
		ni=NumObs;
		np=NumPar;
		t . setZero(np,np);
		I = Eigen::MatrixXd::Identity(ni,ni);
		H . setZero(ni,ni);
		e . setZero(ni);
	}
};

// return 1 if good.
int multifit_linear_validate(const Eigen::MatrixXd & X,
							 const Eigen::VectorXd & y,
							 const int p_start); // if X[][0] is constant term, skip it by p_start=1, otherwise 0

void multifit_linear(const Eigen::MatrixXd		& X,
					 const Eigen::VectorXd		& y,
					 Eigen::VectorXd			& c,
					 Eigen::MatrixXd			& cov,
					 double						& chisq,
					 multifit_linear_workspace	& work );

// Do regression with the constance term. Return p-value for the first beta.
double pv_1st_linear	(const Eigen::MatrixXd	& X, const Eigen::VectorXd	& y, bool one_sided);
double pv_1st_logistic	(const Eigen::MatrixXd	& X, const Eigen::VectorXd	& y, bool one_sided);
double or_1st_logistic	(const Eigen::MatrixXd	& X, const Eigen::VectorXd	& y, bool one_sided);

#endif
