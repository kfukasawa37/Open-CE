####Control-removal design catch-effort model
setwd("L:\\catcheffort")	#任意のパスを設定



#ダミーデータ生成
#設定：半径2kmの円形の捕獲事業区の外側に、半径5kmのドーナツ型の対照区があり、
#対照区ではモニタリングのみを実施すると仮定する。
#(形状については任意だが、対象区とその外側の移出入は等しく無視できると仮定)
nt<-50	#50日捕獲実施
Eh<-cbind(rep(20,nt),rep(0,nt))	#捕獲努力量。2列T行の行列。1列目が捕獲事業区、2列目が周辺のバッファゾーン。1行1日。
Em<-cbind(rep(20,nt),rep(50,nt))	#モニタリング努力量(カメラ日など)。
A<-c(2^2*pi,5^2*pi-2^2*pi)	#捕獲事業区とバッファの面積
P<-2*2*pi	#捕獲事業区の周囲長

#真のパラメータ
logc<-(-4.5)	#対数捕獲効率
c<-exp(logc)	#捕獲効率
logm<-(-5)	#対数検出効率
m<-exp(logm)	#検出効率
logd<-(-6)	#ロジット移出入率
d<-exp(logd)	#移出入率(/km)
N0<-round(A*20)	#初期個体数(20頭/km2)

#データ生成

C<-Y<-Ne<-M<-matrix(NA,nt,2)
N<-matrix(NA,nt+1,2)
N[1,]<-N0
for(t in 1:nt){
	Ne[t,]<-rpois(2,N[t,]/A*P*d)
	M[t,]<-N[t,]-Ne[t,]+rev(Ne[t,])
	Y[t,]<-rpois(2,M[t,]*Em[t,]*m)
	ph<-1-exp(-c*Eh[t,]/A)
	C[t,]<-rbinom(2,M[t,],ph)
	N[t+1,]<-M[t,]-C[t,]
}





####推定
library(coda)
library(parallel)
library(R2jags)


model.file<-"cdce191210.bug"
sink(model.file)
cat("model{
	#observation
	for(t in 1:nt){
		for(i in 1:2){
			ph[t,i]<-1-exp(-exp(logc)*Eh[t,i]/A[i])
			C[t,i]~dbin(ph[t,i],M[t,i])
			Y[t,i]~dpois(M[t,i]*Em[t,i]*exp(logm))
		}
	}

	N[1,1]<-N0[1]
	N[1,2]<-N0[2]
	#process
	for(t in 1:nt){
		Ne[t,1]~dpois(N[t,1]/A[1]*P*exp(logd))
		Ne[t,2]~dpois(N[t,2]/A[2]*P*exp(logd))
		M[t,1]<-N[t,1]-Ne[t,1]+Ne[t,2]
		M[t,2]<-N[t,2]-Ne[t,2]+Ne[t,1]
		N[t+1,1:2]<-M[t,1:2]-C[t,1:2]
	}
	
	#parameter
	logm~dnorm(0,0.0001)
	logd~dnorm(0,0.0001)
	logc~dnorm(0,0.0001)
	N0[1]<-round(n0[1]*A[1])
	N0[2]<-round(n0[2]*A[2])
	log(n0[1])<-log.n0[1]
	log(n0[2])<-log.n0[2]
	log.n0[1]~dt(0,sqrt(1/5),1)
	log.n0[2]~dt(0,sqrt(1/5),1)

}
")
sink()


data<-list(C=C,
			Y=Y,
            Em=Em,
			Eh=Eh,
			A=A,
			P=P,
			nt=nt
			)

inits<-function(){
	list(logm=rnorm(1,-5,1),
			logd=rnorm(1,-5,1),
			logc=rnorm(1,-5,1),
			log.n0=rnorm(2,log(20),1)
		)
}

parameters<-c("logm","logd","logc","log.n0","N","M","Ne")


ctdce.res<-jags.parallel(model.file = model.file,
                           data = data,
                           inits = inits,
                           parameters.to.save = parameters,n.chains=4,
                           n.burnin = 100000,
                           n.iter = 200000,
                           n.thin = 100)



