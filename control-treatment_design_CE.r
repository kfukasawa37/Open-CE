####Control-removal design catch-effort model
setwd("L:\\catcheffort")	#�C�ӂ̃p�X��ݒ�



#�_�~�[�f�[�^����
#�ݒ�F���a2km�̉~�`�̕ߊl���Ƌ�̊O���ɁA���a5km�̃h�[�i�c�^�̑ΏƋ悪����A
#�ΏƋ�ł̓��j�^�����O�݂̂����{����Ɖ��肷��B
#(�`��ɂ��Ă͔C�ӂ����A�Ώۋ�Ƃ��̊O���̈ڏo���͓����������ł���Ɖ���)
nt<-50	#50���ߊl���{
Eh<-cbind(rep(20,nt),rep(0,nt))	#�ߊl�w�͗ʁB2��T�s�̍s��B1��ڂ��ߊl���Ƌ�A2��ڂ����ӂ̃o�b�t�@�]�[���B1�s1���B
Em<-cbind(rep(20,nt),rep(50,nt))	#���j�^�����O�w�͗�(�J�������Ȃ�)�B
A<-c(2^2*pi,5^2*pi-2^2*pi)	#�ߊl���Ƌ�ƃo�b�t�@�̖ʐ�
P<-2*2*pi	#�ߊl���Ƌ�̎��͒�

#�^�̃p�����[�^
logc<-(-4.5)	#�ΐ��ߊl����
c<-exp(logc)	#�ߊl����
logm<-(-5)	#�ΐ����o����
m<-exp(logm)	#���o����
logd<-(-6)	#���W�b�g�ڏo����
d<-exp(logd)	#�ڏo����(/km)
N0<-round(A*20)	#�����̐�(20��/km2)

#�f�[�^����

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





####����
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



