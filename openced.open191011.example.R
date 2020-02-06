########Open-CE使用例
source("L:/catcheffort/openced.open191011.R",chdir=T)

#ダミーデータ生成
#設定：半径2kmの円形の捕獲事業区
nt<-50	#50日捕獲実施
Eh<-rep(20,nt)	#捕獲努力量。
Em<-rep(20,nt)	#モニタリング努力量(カメラ日など)。
A<-c(2^2*pi)	#捕獲事業区の面積

#真のパラメータ
logc<-(-4.5)	#対数捕獲効率
c<-exp(logc)	#捕獲効率
logm<-(-5)	#対数検出効率
m<-exp(logm)	#検出効率
logitd<-(-6)	#ロジット移出入率
d<-exp(logitd)	#移出入率(/km)
N0<-round(A*20)	#初期個体数(20頭/km2)

#データ生成

C<-Y<-Ne<-Ni<-M<-rep(NA,nt)
N<-rep(NA,nt+1)
N[1]<-N0
for(t in 1:nt){
	Ne[t]<-rbinom(1,N[t],d)
	Ni[t]<-rpois(1,N[t]*d)
	M[t]<-N[t]-Ne[t]+Ni[t]
	Y[t]<-rpois(1,M[t]*Em[t]/A*m)
	ph<-1-exp(-c*Eh[t]/A)
	C[t]<-rbinom(1,M[t],ph)
	N[t+1]<-M[t]-C[t]
}

#推定(overdispersion)
res1<-opence(C=C,Eh=Eh,A=A,Y=Y,Em=Em,capod=TRUE,monod=TRUE,Smax=1000)

#推定(without overdispersion)
res2<-opence(C=C,Eh=Eh,A=A,Y=Y,Em=Em,capod=FALSE,monod=FALSE,Smax=1000)

