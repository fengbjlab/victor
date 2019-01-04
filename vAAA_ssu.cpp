#include <Eigen/Eigenvalues>
#include <tft/libfbj_math.hpp>
#include "vAAA_ssu.hpp"

/* R code from http://www.biostat.umn.edu/~weip/prog/BasuPanGE11/aSumTest.r
 # Y: disease lables; =0 for controls, =1 for cases;
 # X: genotype data; row for subject, column for SNP (#copies for an allele);
 SumSqUs_Test <- function (Y, X, B=500, alpha0=0.1) {
 Xg <- X 
 Xbar<-apply(Xg, 2, mean) 
 Xgb<-Xg
 for(i in 1:nrow(Xg)) Xgb[i,]<-Xg[i,]-Xbar
 U<-t(Xg) %*% (Y-mean(Y))
 
 #SumSqU:
 Tg1<- t(U) %*% U
 ##cov of the score stats:
 CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)
 ##distr of Tg1 is sum of cr Chisq_1:
 cr<-eigen(CovS, only.values=TRUE)$values
 ##approximate the distri by alpha Chisq_d + beta:
 alpha1<-sum(cr*cr*cr)/sum(cr*cr)
 beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
 d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
 alpha1<-as.double(alpha1)
 beta1<-as.double(beta1)
 d1<-as.double(d1)
 pTg1<-as.numeric(1-pchisq((Tg1-beta1)/alpha1, d1))
 
 #SumSqUw:
 diagCovS<-diag(CovS)
 diagCovS<-ifelse(diagCovS>1e-10, diagCovS, 1e-10)
 Tg2<- t(U) %*%  diag(1/diagCovS) %*% U
 ##distr of Tg1 is sum of cr Chisq_1:
 cr<-eigen(CovS %*% diag(1/diagCovS), only.values=TRUE)$values
 ##approximate the distri by alpha Chisq_d + beta:
 alpha2<-sum(cr*cr*cr)/sum(cr*cr)
 beta2<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
 d2<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
 alpha2<-as.double(alpha2)
 beta2<-as.double(beta2)
 d2<-as.double(d2)
 pTg2<-as.numeric(1-pchisq((Tg2-beta2)/alpha2, d2))
 
 return(c(pTg1,pTg2))
 }
 to use: 
 Y<-as.matrix(read.table("y"))
 X<-as.matrix(read.table("x"))
 SumSqUs_Test(Y,X)
 */

// approximate the distribution by   alpha * Chisq(d) + beta
// cr is vector of eigen values?
double scale_shift_chisq_q(double x2, Eigen::VectorXd& cr)
{
	Eigen::VectorXd bk = cr;		double sum_cr1 = cr.sum();
	cr = cr.array() * cr.array();	double sum_cr2 = cr.sum();
	cr = cr.array() * bk.array();	double sum_cr3 = cr.sum();
	double alpha = sum_cr3 / sum_cr2;
	double beta  = sum_cr1 - (sum_cr2 * sum_cr2)/sum_cr3;
	double d = (sum_cr2*sum_cr2*sum_cr2)/(sum_cr3*sum_cr3);
	return cdf_chisq_q((x2-beta)/alpha, d); // p-value
}

int ssu_validate(const Eigen::MatrixXd & X,
				 const Eigen::VectorXd & y)
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	
	// check validity of matrix X
	bool vary = false;
	for (int j=0;j<np;j++)
	{
		double minV= DBL_MAX;
		double maxV=-DBL_MAX;
		for (int i=0;i<ni;i++)
		{
			double value = X(i,j);
			if (value>maxV) maxV=value;
			if (value<minV) minV=value;
		}
		if (maxV!=minV) { vary=true; break; }
	}
	return vary;
}

// ----------------- SumSqUs test -----------------

double SumSqUs_test(const Eigen::MatrixXd &X, const Eigen::VectorXd & Y, SumSqUs_workspace & ws)
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	ws.Ybar = Y.mean();
	ws.Xbar = X.colwise().mean();
	ws.Y_Ybar = Y.array() - ws.Ybar;
	for (int col=0; col<np; ++col)
	{
		double Xb = ws.Xbar(col);
		for (int row=0; row<ni; ++row) ws.X_Xbar(row,col) = X(row,col)-Xb;
	}
	ws.U = X.transpose() * ws.Y_Ybar;
	double Tg1 = ws.U.transpose() * ws.U;
	ws.CovS = ws.Ybar * (1-ws.Ybar) * ws.X_Xbar.transpose() * ws.X_Xbar;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ws.CovS,false);
	if (es.info() != Eigen::Success) return std::numeric_limits<double>::signaling_NaN(); // exit_error("Eigen solver failed.");;
	ws.cr = es.eigenvalues();
	double p = scale_shift_chisq_q(Tg1,ws.cr); // pTg1
	return p;
}

// ----------------- SumSqUw test -----------------

// SSUw test, weight is the inverse of the sample SD of the variant
double SumSqUw_test(const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, SumSqUw_Xstable& ws)
{
	const int ni = X.rows();	// number of individuals (observations)
	const int np = X.cols();	// number of parameters
	Eigen::VectorXd Y_Ybar;		// (ni)
	Eigen::VectorXd U;			// (np)
	Eigen::VectorXd cr;			// (np)
	
	ws.Ybar = Y.mean();
	ws.Xbar = X.colwise().mean();
	Y_Ybar = Y.array() - ws.Ybar;
	for (int col=0; col<np; col++)
	{
		double Xb = ws.Xbar(col);
		for (int row=0; row<ni; row++) ws.X_Xbar(row,col) = X(row,col)-Xb;
	}
	ws.Xtrn = X.transpose();
	U = ws.Xtrn * Y_Ybar;
	ws.CovS = ws.Ybar * (1-ws.Ybar) * ws.X_Xbar.transpose() * ws.X_Xbar;
	ws.diagCovS.setZero(np,np);
	for (int i=0; i<np; i++)
	{
		double d = ws.CovS(i,i);
		if (d<=1e-10) d=1e-10;
		ws.diagCovS(i,i)=1/d;
	}
	double Tg2 = U.transpose() * ws.diagCovS * U;
	Eigen::EigenSolver<Eigen::MatrixXd> es(ws.CovS * ws.diagCovS,false);
	if (es.info() != Eigen::Success) return std::numeric_limits<double>::signaling_NaN(); // exit_error("Eigen solver failed.");;
	cr = es.eigenvalues().real();
//	std::sort(cr.data(),cr.data()+cr.size(),std::greater<double>()); // sorting seems to has no effect
	return scale_shift_chisq_q(Tg2,cr); // p-valaue
}

// SSUw test with customized weights
double SumSqUw_test(const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, const std::vector<double> &W, SumSqUw_Xstable& ws)
{
	const int ni = X.rows();	// number of individuals (observations)
	const int np = X.cols();	// number of parameters

	if ((int)W.size()!=(int)X.cols()) exit_error("Number of weights doesn't match number of variables.");
	ws.diagCovS.setZero(np,np);
	bool allZero = true;
	for (int i=0; i<np; i++) { ws.diagCovS(i,i)=W[i]; if (W[i]) allZero=false; }
	if (allZero) exit_error("Weights for SSU cannot be all zeros.");
	
	Eigen::VectorXd Y_Ybar;		// (ni)
	Eigen::VectorXd U;			// (np)
	Eigen::VectorXd cr;			// (np)
	
	ws.Ybar = Y.mean();
	ws.Xbar = X.colwise().mean();
	Y_Ybar = Y.array() - ws.Ybar;
	for (int col=0; col<np; col++)
	{
		double Xb = ws.Xbar(col);
		for (int row=0; row<ni; row++) ws.X_Xbar(row,col) = X(row,col)-Xb;
	}
	ws.Xtrn = X.transpose();
	U = ws.Xtrn * Y_Ybar;
	ws.CovS = ws.Ybar * (1-ws.Ybar) * ws.X_Xbar.transpose() * ws.X_Xbar;
	double Tg2 = U.transpose() * ws.diagCovS * U;
	Eigen::EigenSolver<Eigen::MatrixXd> es(ws.CovS * ws.diagCovS,false);
	if (es.info() != Eigen::Success) return std::numeric_limits<double>::signaling_NaN(); // exit_error("Eigen solver failed.");;
	cr = es.eigenvalues().real();
	//	std::sort(cr.data(),cr.data()+cr.size(),std::greater<double>()); // sorting seems to has no effect
	return scale_shift_chisq_q(Tg2,cr); // p-valaue
}

// for simulations, use the same WorkSpace, X unchanged while Y changed.
double SumSqUw_test(const Eigen::VectorXd &Y, SumSqUw_Xstable& ws)
{
	Eigen::VectorXd Y_Ybar;		// (ni)
	Eigen::VectorXd U;			// (np)
	Eigen::VectorXd cr;			// (np)

	Y_Ybar = Y.array() - ws.Ybar;
	U = ws.Xtrn * Y_Ybar;
	double Tg2 = U.transpose() * ws.diagCovS * U;
	Eigen::EigenSolver<Eigen::MatrixXd> es(ws.CovS * ws.diagCovS,false);
	if (es.info() != Eigen::Success) return std::numeric_limits<double>::signaling_NaN(); // exit_error("Eigen solver failed.");;
	cr = es.eigenvalues().real();
	//	std::sort(cr.data(),cr.data()+cr.size(),std::greater<double>()); // sorting seems to has no effect
	return scale_shift_chisq_q(Tg2,cr); // p-valaue
}

std::vector<double> pv_SSUw_default_empty_weights;

double pv_SSUw(const Eigen::MatrixXd & X, const Eigen::VectorXd	& y, const std::vector<double> &W)
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	if (ssu_validate(X,y))
	{
		SumSqUw_Xstable work(ni,np);
		double pv = SumSqUw_test(X,y,work);
		return pv;
	}
	return std::numeric_limits<double>::signaling_NaN();
}

double pv_SSUc(const Eigen::MatrixXd & X, const Eigen::VectorXd	& y, const std::vector<double> &W)
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	if (ssu_validate(X,y))
	{
		SumSqUw_Xstable work(ni,np);
		double pv = SumSqUw_test(X,y,W,work);
		return pv;
	}
	return std::numeric_limits<double>::signaling_NaN();
}

