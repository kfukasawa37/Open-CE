########Open-CEgpá
source("L:/catcheffort/openced.open191011.R",chdir=T)

#_~[f[^¶¬
#ÝèF¼a2kmÌ~`ÌßlÆæ
nt<-50	#50úßlÀ{
Eh<-rep(20,nt)	#ßlwÍÊB
Em<-rep(20,nt)	#j^OwÍÊ(JúÈÇ)B
A<-c(2^2*pi)	#ßlÆæÌÊÏ

#^Ìp[^
logc<-(-4.5)	#Îßlø¦
c<-exp(logc)	#ßlø¦
logm<-(-5)	#Îoø¦
m<-exp(logm)	#oø¦
logitd<-(-6)	#WbgÚoü¦
d<-exp(logitd)	#Úoü¦(/km)
N0<-round(A*20)	#úÂÌ(20ª/km2)

#f[^¶¬

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

#è(overdispersion)
res1<-opence(C=C,Eh=Eh,A=A,Y=Y,Em=Em,capod=TRUE,monod=TRUE,Smax=1000)

#è(without overdispersion)
res2<-opence(C=C,Eh=Eh,A=A,Y=Y,Em=Em,capod=FALSE,monod=FALSE,Smax=1000)

