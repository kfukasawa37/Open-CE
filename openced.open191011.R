##########�J���̌QCatch-effort model
#����
#����R�R�[�h��C++�\�[�X�R�[�hopenced.open190906.cpp�𓯈�t�H���_�Ɋi�[���A
#R��ʂ�source("{�p�X��}/openced.open191011.R",chdir=T)�����s����Ɗ֐������p�\�ɂȂ�

#���C�u�����Ǎ�
if(!require(Rcpp)){
	install.packages("Rcpp")
	require(Rcpp)
}

if(!require(RcppEigen)){
	install.packages("RcppEigen")
	require(RcppEigen)
}

#Rcpp�R���p�C��
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("openced.open190906.cpp")


#�ΐ��ޓx�v�Z���b�p�[�֐�
opence_lf2<-function(par,Cmat,Emat,A,Coth,iscap,Smax=1000,capod=TRUE,monod=TRUE){
	par2<-par
	if(!capod){
		par2<-c(par2,rep(log(10000),sum(iscap)))
	}
	if(!monod){
		par2<-c(par2,rep(log(10000),sum(1-iscap)))
	}
	res<-opence_lf(par2,Cmat,Emat,Coth,iscap,Smax,A)

	return(res)
}


#####OpenCE����֐�
opence<-function(C,Eh,A,Y=NULL,Em=NULL,Coth=NULL,capod=TRUE,monod=TRUE,initpar=NULL,Smax=1000){
	Cmat<-as.matrix(C)
	Cname<-colnames(Cmat)
	colnames(Cmat)<-NULL
	nC<-ncol(Cmat)
	if(is.null(Cname)){
		Cname<-1:nC
	}
	ntC<-nrow(Cmat)

	Ehmat<-as.matrix(Eh)
	Ehname<-colnames(Ehmat)
	colnames(Ehmat)<-NULL
	nEh<-ncol(Ehmat)
	ntEh<-nrow(Ehmat)

	if((nC!=nEh)|(ntC!=ntEh)){
		stop("size of C and Eh are not the same.")
	}

    nt<-ntC
	CYmat<-Cmat


	if(is.null(Y)){
		mobs<-NULL
        nY<-0
		Yname<-NULL
		CYmat<-Cmat
		Ehmmat<-Ehmat
	}else{
		Ymat<-as.matrix(Y)
		Yname<-colnames(Ymat)
		colnames(Ymat)<-NULL
		nY<-ncol(Ymat)
		ntY<-nrow(Ymat)
		
		Emmat<-as.matrix(Em)
		Emname<-colnames(Emmat)
		colnames(Emmat)<-NULL
		nEm<-ncol(Emmat)
		ntEm<-nrow(Emmat)

		if((nY!=nEm)|(ntY!=ntEm)){
			stop("size of Y and Em are not the same.")
		}
    	if(ntC!=ntY){
	    	stop("row number of C and Y are not the same.")
	    }
		if(is.null(Yname)){
			Yname<-1:nY
		}

		CYmat<-cbind(Cmat,Ymat)
		Ehmmat<-cbind(Ehmat,Emmat)
	}
	
	iscap<-c(rep(1,nC),rep(0,nY))
	
	if(is.null(Coth)){
		Coth<-rep(0,ntC)
	}
	

	if(ntC!=length(Coth)){
		stop("row number of C and Y are not the same.")
	}

	if(is.null(initpar)){
		initpar<-c(3,-6,rep(-5,nC+nY),rep(0,nC*capod+nY*monod))
	}
	
	res<-nlm(function(...){opence_lf2(...)*(-1)},initpar,hessian=T,iterlim=1000,print.level=2,Cmat=CYmat,Emat=Ehmmat,A=A,Coth=Coth,iscap=iscap,Smax=Smax,capod=capod,monod=monod)
	parnames<-c("logD","logitem",
								paste("logc_",Cname,sep=""),
								ifelse(!is.null(Yname),paste("logs_",Yname,sep=""),NULL),
								ifelse(capod,paste("logbeta_",Cname,sep=""),NULL),
								ifelse(capod&(!is.null(Yname)),paste("logphi_",Yname,sep=""),NULL)
							)
	names(res$estimate)<-parnames
	dimnames(res$hessian)<-list(parnames,parnames)
	return(res)
}


