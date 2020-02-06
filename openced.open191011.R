##########開放個体群Catch-effort model
#準備
#このRコードとC++ソースコードopenced.open190906.cppを同一フォルダに格納し、
#R画面でsource("{パス名}/openced.open191011.R",chdir=T)を実行すると関数が利用可能になる

#ライブラリ読込
if(!require(Rcpp)){
	install.packages("Rcpp")
	require(Rcpp)
}

if(!require(RcppEigen)){
	install.packages("RcppEigen")
	require(RcppEigen)
}

#Rcppコンパイル
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("openced.open190906.cpp")


#対数尤度計算ラッパー関数
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


#####OpenCE推定関数
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

	if(!is.null(Yname)){
		name.logs<-paste("logs_",Yname,sep="")
	}else{
		name.logs<-character(0)
	}
	
	if(capod){
		name.logbeta<-paste("logbeta_",Cname,sep="")
	}else{
		name.logbeta<-character(0)
	}

	if(capod&(!is.null(Yname))){
		name.logphi<-paste("logphi_",Yname,sep="")
	}else{
		name.logphi<-character(0)
	}

	parnames<-c("logD","logitem",
								paste("logc_",Cname,sep=""),
								name.logs,
								name.logbeta,
								name.logphi
							)
	names(res$estimate)<-parnames
	dimnames(res$hessian)<-list(parnames,parnames)
	return(res)
}



