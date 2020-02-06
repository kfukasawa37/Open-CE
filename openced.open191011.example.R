########Open-CE�g�p��
source("L:/catcheffort/openced.open191011.R",chdir=T)

#�_�~�[�f�[�^����
#�ݒ�F���a2km�̉~�`�̕ߊl���Ƌ�
nt<-50	#50���ߊl���{
Eh<-rep(20,nt)	#�ߊl�w�͗ʁB
Em<-rep(20,nt)	#���j�^�����O�w�͗�(�J�������Ȃ�)�B
A<-c(2^2*pi)	#�ߊl���Ƌ�̖ʐ�

#�^�̃p�����[�^
logc<-(-4.5)	#�ΐ��ߊl����
c<-exp(logc)	#�ߊl����
logm<-(-5)	#�ΐ����o����
m<-exp(logm)	#���o����
logitd<-(-6)	#���W�b�g�ڏo����
d<-exp(logitd)	#�ڏo����(/km)
N0<-round(A*20)	#�����̐�(20��/km2)

#�f�[�^����

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

#����(overdispersion)
res1<-opence(C=C,Eh=Eh,A=A,Y=Y,Em=Em,capod=TRUE,monod=TRUE,Smax=1000)

#����(without overdispersion)
res2<-opence(C=C,Eh=Eh,A=A,Y=Y,Em=Em,capod=FALSE,monod=FALSE,Smax=1000)
