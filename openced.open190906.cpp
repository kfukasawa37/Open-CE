#include <Rcpp.h>
using Rcpp::as;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, integer

double capt(double c, double E,double area)
{
	return 1-exp(-c*E/area);
}

// [[Rcpp::export]]
double sum_lp(VectorXd lp)
{
	double max_lp=lp.maxCoeff();
	int size_lp=lp.size();
	VectorXd explp(size_lp);
	for(int i=0;i<size_lp;i++){
		explp(i)=exp(lp(i)-max_lp);
        //explp(i)=exp(lp(i));
	}
	return log(explp.sum())+max_lp;
    //return log(explp.sum());

}

// [[Rcpp::export]]
double dbetabin_p_log(int k,int n, double p, double beta){
    if((k==0)&(n==0)){
        return 0;
    }else if(k>n){
        return -99999;
    }else{
        double alpha=p*beta/(1-p);
	    double res=std::lgamma(n+1)+std::lgamma(k+alpha)+std::lgamma(n-k+beta)+std::lgamma(alpha+beta)-std::lgamma(k+1)-std::lgamma(n-k+1)-std::lgamma(n+alpha+beta)-std::lgamma(alpha)-std::lgamma(beta);
        return(res);
    }
}

// [[Rcpp::export]]
double dgammapois_mu_log(int k, double mu, double phi){
	if(mu==0&k==0){
		return 0;
	}
	if(mu==0&k!=0){
		return -99999;
	}
    if(k<0){
        return -99999;
     }else{
	    double res=std::lgamma(k+phi)-std::lgamma(phi)-std::lgamma(k+1)+k*log(mu)-k*log(mu+phi)+phi*log(phi)-phi*log(mu+phi);
        return(res);
    }
}


// [[Rcpp::export]]
VectorXd std_lp(VectorXd lp)
{
	int nlp=lp.size();
	VectorXd p(nlp);
	for(int i=0;i<nlp;i++){
		p(i) = exp(lp(i)-lp.maxCoeff());
	}
	double p_sum=p.sum();
	VectorXd std_p(nlp);
	for(int i=0;i<nlp;i++){
		std_p(i)=log(p(i)/p_sum);
	}
	return std_p;
}

//use R::dbinom R::dpois (not Rcpp::) when scalar is applied.



// [[Rcpp::export]]
double opence_lf(	NumericVector par0,
				NumericMatrix C0,
				NumericMatrix CE0,
				NumericVector Coth0,
				NumericVector iscap,
				int Smax,
				double area)
{
	VectorXd par(as<VectorXd> (par0));
	MatrixXd C(as<MatrixXd> (C0));
	MatrixXd CE(as<MatrixXd> (CE0));
	VectorXd Coth(as<VectorXd> (Coth0));
	int npar=par.size();
	int ny=C.rows();
	int nc=C.cols();
	double logD=par(0);
	double logitem=par(1);
	VectorXd logc=par.segment(2,nc);
    VectorXd logbeta=par.tail(nc);
		
	double D=exp(logD);
	double em=1/(1+exp((-1)*logitem));
	double f=0;
	
	//Transition Probability 1:emigration
	MatrixXd TP1(Smax+1,Smax+1);
	double p=1-em;
	for(int i=0;i<Smax+1;i++){
		for(int j=0;j<Smax+1;j++){
			TP1(j,i)=R::dbinom(i,j,p,true);
		}
	}

	//Transition Probability 2:immigration
	MatrixXd TP2(Smax+1,Smax+1);
	//VectorXd stdTP(Smax+1);
	TP2.fill(-99999);
	for(int i=0;i<Smax+1;i++){
		int x=0;
		for(int j=i;j<Smax+1;j++){
			TP2(i,j)=R::dpois(x,D*area*em,true);
			x=x+1;
		}
		//stdTP=std_lp(TP2.block(0,i))
	}

	//t=0
	//set prior
	MatrixXd filt(Smax+1,ny+1);
    
	for(int i=0;i<Smax+1;i++){
		filt(i,0)=R::dpois(i,D*area,true);
	}


	//Recursive filtering
	MatrixXd temp(Smax+1,Smax+1);
	MatrixXd pred(Smax+1,ny);
    MatrixXd monit(Smax+1,ny);
	int Ctotal;
	for(int t=0;t<ny;t++){
		Ctotal=Coth(t);
		for(int j=0;j<nc;j++){
			Ctotal=Ctotal+C(t,j)*iscap(j);
		}

		for(int i=0;i<Smax+1;i++){
			for(int j=0;j<Smax+1;j++){
				temp(i,j)=TP1(j,i)+filt(j,t);
			}
 		}
        for(int i=0;i<Smax+1;i++){
    		pred(i,t)=sum_lp(temp.row(i));
            if(std::isnan(pred(i,t))){
                pred(i,t)=-99999;
            }
            //monit(i,t)=pred(i,t);
        }

		pred.col(t)=std_lp(pred.col(t));
		for(int i=0;i<Smax+1;i++){
			for(int j=0;j<Smax+1;j++){
				temp(i,j)=TP2(j,i)+pred(j,t);
			}
			pred(i,t)=sum_lp(temp.row(i));
            if(std::isnan(pred(i,t))){
                pred(i,t)=-99999;
            }
		}
		pred.col(t)=std_lp(pred.col(t));
		//Rcpp::Rcout << pred(100) <<"\n";
		
		for(int i=0;i<Smax+1;i++){
			filt(i,t+1)=pred(i,t);
		}
		
		for(int j=0;j<nc;j++){
			if(iscap(j)==1){
				if(CE(t,j)>0){
					for(int i=0;i<Smax+1;i++){
						filt(i,t+1)=filt(i,t+1)+dbetabin_p_log(C(t,j),i,capt(exp(logc(j)),CE(t,j),area),exp(logbeta(j)));
					}
				}
			}else{
				if(CE(t,j)>0){
					for(int i=0;i<Smax+1;i++){
						filt(i,t+1)=filt(i,t+1)+dgammapois_mu_log(C(t,j),exp(logc(j))*i*CE(t,j)/area,exp(logbeta(j)));
					}
				}
			}
		}
		f += sum_lp(filt.col(t+1));
		if(Ctotal>0){
			for(int i=0;i<Smax-Ctotal+1;i++){
				filt(i,t+1) = filt(i+Ctotal,t+1);
			}
			for(int i=0;i<Ctotal;i++){
				filt(i+Smax-Ctotal+1,t+1) = -99999;
			}
		}	
		filt.col(t+1)=std_lp(filt.col(t+1));

	}
	return f;
	//return monit;
}
