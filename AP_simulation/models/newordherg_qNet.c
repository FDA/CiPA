// Copyright (c) 2011-2015 by Thomas O'Hara, Yoram Rudy, 
//                            Washington University in St. Louis.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the names of the copyright holders nor the names of its
// contributors may be used to endorse or promote products derived from 
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//

// C implementation of the IKr-dynamic O'Hara-Rudy dynamic (ORd) model for
// the undiseased human ventricular action potential and calcium transient.
// This code can be compiled into a dynamically linked library (DLL) for use
// in the R software environment (https://www.R-project.org/) with the
// following R command:
//      system("R CMD SHLIB newordherg_qNet.c")
//
// The IKr-dynamic ORd model was modified from the ORd model to
// include a Markov model for IKr gating and drug binding kinetics. It is
// described in the article "Improving the in silico assessment of
// proarrhythmia risk by combining hERG (Human Ether-à-go-go-Related Gene)
// channel-drug binding kinetics and multichannel pharmacology"
// by Zhihua Li, Sara Dutta, Jiansong Sheng, Phu N. Tran, Wendy Wu, Kelly
// Chang, Thembi Mdluli, David G. Strauss, and Thomas Colatsky
// Link: http://circep.ahajournals.org/content/10/2/e004628.abstract
//
// This version of the model code includes the integration of an additional
// variable to compute qNet, a metric which is described in the article
// "Optimization of an In silico Cardiac Cell Model for Proarrhythmia Risk
// Assessment"
// by Sara Dutta, Kelly C. Chang, Kylie A. Beattie, Jiansong Sheng, Phu N.
// Tran, Wendy W. Wu, Min Wu, David G. Strauss, Thomas Colatsky, and Zhihua
// Li
// Link: https://www.frontiersin.org/article/10.3389/fphys.2017.00616
//

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <time.h>
static double parms[59];   //celltype plus 7 channels plus hERGMM parameters
#define celltype parms[0]
#define GKrfc parms[1]
#define GNaLfc parms[2]
#define GNafc parms[3]
#define GKsfc parms[4]
#define GK1fc parms[5]
#define PCafc parms[6]
#define Gtofc parms[7]
#define A1 parms[8]
#define B1 parms[9]
#define q1 parms[10]
#define A2 parms[11]
#define B2 parms[12]
#define q2 parms[13]
#define A3 parms[14]
#define B3 parms[15]
#define q3 parms[16]
#define A4 parms[17]
#define B4 parms[18]
#define q4 parms[19]
#define A11 parms[20]
#define B11 parms[21]
#define q11 parms[22]
#define A21 parms[23]
#define B21 parms[24]
#define q21 parms[25]
#define A31 parms[26]
#define B31 parms[27]
#define q31 parms[28]
#define A41 parms[29]
#define B41 parms[30]
#define q41 parms[31]
#define A51 parms[32]
#define B51 parms[33]
#define q51 parms[34]
#define A52 parms[35]
#define B52 parms[36]
#define q52 parms[37]
#define A53 parms[38]
#define B53 parms[39]
#define q53 parms[40]
#define A61 parms[41]
#define B61 parms[42]
#define q61 parms[43]
#define A62 parms[44]
#define B62 parms[45]
#define q62 parms[46]
#define A63 parms[47]
#define B63 parms[48]
#define q63 parms[49]
#define Kmax parms[50]
#define Ku parms[51]
#define n parms[52]
#define halfmax parms[53]
#define Kt parms[54]
#define Vhalf parms[55]
#define Temp parms[56]
#define ko parms[57]
#define amp parms[58]



void initmod(void (* odeparms)(int *, double *)){
int N=59;
odeparms(&N, parms);}

void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
if (ip[0] < 8 ) error("nout not enough!");   //7 currents as additional output plus dv (ydot[0])

//extracellular ionic concentrations
double nao=140.0;
double cao=1.8;

//physical constants
double R=8314.0;
double T=310.0;
double F=96485.0;

//cell geometry
double L=0.01;
double rad=0.0011;
double vcell=1000*3.14*rad*rad*L;
double Ageo=2*3.14*rad*rad+2*3.14*rad*L;
double Acap=2*Ageo;
double vmyo=0.68*vcell;
double vnsr=0.0552*vcell;
double vjsr=0.0048*vcell;
double vss=0.02*vcell;

//give names to the state vector values
double v=y[0];
double nai=y[1];
double nass=y[2];
double ki=y[3];
double kss=y[4];
double cai=y[5];
double cass=y[6];
double cansr=y[7];
double cajsr=y[8];
double m=y[9];
double hf=y[10];
double hs=y[11];
double j=y[12];
double hsp=y[13];
double jp=y[14];
double mL=y[15];
double hL=y[16];
double hLp=y[17];
double a=y[18];
double iF=y[19];
double iS=y[20];
double ap=y[21];
double iFp=y[22];
double iSp=y[23];
double d=y[24];
double ff=y[25];
double fs=y[26];
double fcaf=y[27];
double fcas=y[28];
double jca=y[29];
double nca=y[30];
double ffp=y[31];
double fcafp=y[32];
double xs1=y[33];
double xs2=y[34];
double xk1=y[35];
double Jrelnp=y[36];
double Jrelp=y[37];
double CaMKt=y[38];
double IC1=y[39];
double IC2=y[40];
double C1=y[41];
double C2=y[42];
double O=y[43];
double IO=y[44];
double IObound=y[45];
double Obound=y[46];
double Cbound=y[47];
double D=y[48];

//CaMK constants
double KmCaMK=0.15;

double aCaMK=0.05;
double bCaMK=0.00068;
double CaMKo=0.05;
double KmCaM=0.0015;
//update CaMK
double CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
double CaMKa=CaMKb+CaMKt;
ydot[38]=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt;

//reversal potentials
double ENa=(R*T/F)*log(nao/nai);
double EK=(R*T/F)*log(ko/ki);
double PKNa=0.01833;
double EKs=(R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai));

//convenient shorthand calculations
double vffrt=v*F*F/(R*T);
double vfrt=v*F/(R*T);

//calculate INa
double mss=1.0/(1.0+exp((-(v+39.57))/9.871));
double tm=1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
ydot[9]=(mss-m)/tm;
double hss=1.0/(1+exp((v+82.90)/6.086));
double thf=1.0/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
double ths=1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
double Ahf=0.99;
double Ahs=1.0-Ahf;
ydot[10]=(hss-hf)/thf;
ydot[11]=(hss-hs)/ths;
double h=Ahf*hf+Ahs*hs;
double jss=hss;
double tj=2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
ydot[12]=(jss-j)/tj;
double hssp=1.0/(1+exp((v+89.1)/6.086));
double thsp=3.0*ths;
ydot[13]=(hssp-hsp)/thsp;
double hp=Ahf*hf+Ahs*hsp;
double tjp=1.46*tj;
ydot[14]=(jss-jp)/tjp;
double GNa=75/GNafc;
double fINap=(1.0/(1.0+KmCaMK/CaMKa));
double INa=GNa*(v-ENa)*pow(m,3.0)*((1.0-fINap)*h*j+fINap*hp*jp);

//calculate INaL
double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
double tmL=tm;
ydot[15]=(mLss-mL)/tmL;
double hLss=1.0/(1.0+exp((v+87.61)/7.488));
double thL=200.0;
ydot[16]=(hLss-hL)/thL;
double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
double thLp=3.0*thL;
ydot[17]=(hLssp-hLp)/thLp;
double GNaL=0.0075/GNaLfc;
if (celltype==1)
    GNaL=GNaL*0.6;

double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
double INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

//calculate Ito
double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
ydot[18]=(ass-a)/ta;
double iss=1.0/(1.0+exp((v+43.94)/5.711));
double delta_epi=1.0;
if (celltype==1){
    delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
}else{
    delta_epi=1.0;
}
double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
tiF=tiF*delta_epi;
tiS=tiS*delta_epi;
double AiF=1.0/(1.0+exp((v-213.6)/151.2));
double AiS=1.0-AiF;
ydot[19]=(iss-iF)/tiF;
ydot[20]=(iss-iS)/tiS;
double i=AiF*iF+AiS*iS;
double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
ydot[21]=(assp-ap)/ta;
double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
double tiFp=dti_develop*dti_recover*tiF;
double tiSp=dti_develop*dti_recover*tiS;
ydot[22]=(iss-iFp)/tiFp;
ydot[23]=(iss-iSp)/tiSp;
double myip=AiF*iFp+AiS*iSp;
double Gto=0.02/Gtofc;
if (celltype==1){
    Gto=Gto*4.0;
}else if (celltype==2){
    Gto=Gto*4.0;
}
double fItop=(1.0/(1.0+KmCaMK/CaMKa));
double Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*myip);


//calculate ICaL, ICaNa, ICaK
double dss=1.0/(1.0+exp((-(v+3.940))/4.230));
double td=0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
ydot[24]=(dss-d)/td;
double fss=1.0/(1.0+exp((v+19.58)/3.696));
double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
double Aff=0.6;
double Afs=1.0-Aff;
ydot[25]=(fss-ff)/tff;
ydot[26]=(fss-fs)/tfs;
double f=Aff*ff+Afs*fs;
double fcass=fss;
double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
double Afcas=1.0-Afcaf;
ydot[27]=(fcass-fcaf)/tfcaf;
ydot[28]=(fcass-fcas)/tfcas;
double fca=Afcaf*fcaf+Afcas*fcas;
double tjca=75.0;
ydot[29]=(fcass-jca)/tjca;
double tffp=2.5*tff;
ydot[31]=(fss-ffp)/tffp;
double fp=Aff*ffp+Afs*fs;
double tfcafp=2.5*tfcaf;
ydot[32]=(fcass-fcafp)/tfcafp;
double fcap=Afcaf*fcafp+Afcas*fcas;
double Kmn=0.002;
double k2n=1000.0;
double km2n=jca*1.0;
double anca=1.0/(k2n/km2n+pow((1.0+Kmn/cass),4.0));
ydot[30]=anca*k2n-nca*km2n;
double PhiCaL=4.0*vffrt*(cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
double PhiCaNa=1.0*vffrt*(0.75*nass*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
double PhiCaK=1.0*vffrt*(0.75*kss*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
double zca=2.0;
double PCa=0.0001/PCafc;
if (celltype==1){
    PCa=PCa*1.2;
}else if (celltype==2){
    PCa=PCa*2.5;
}
double PCap=1.1*PCa;
double PCaNa=0.00125*PCa;
double PCaK=3.574e-4*PCa;
double PCaNap=0.00125*PCap;
double PCaKp=3.574e-4*PCap;
double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
double ICaL=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
double ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
double ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);


//calculate IKr
double GKr=0.046/GKrfc;
if (celltype==1){
    GKr=GKr*1.3;
}else if (celltype==2){
    GKr=GKr*0.8;
}
ydot[39] = -(A11*exp(B11*v)*IC1*exp((Temp-20)*log(q11)/10)-A21*exp(B21*v)*IC2*exp((Temp-20)*log(q21)/10))+(A51*exp(B51*v)*C1*exp((Temp-20)*log(q51)/10)-A61*exp(B61*v)*IC1*exp((Temp-20)*log(q61)/10));
ydot[40] = (A11*exp(B11*v)*IC1*exp((Temp-20)*log(q11)/10)-A21*exp(B21*v)*IC2*exp((Temp-20)*log(q21)/10))-(A3*exp(B3*y[0])*IC2*exp((Temp-20)*log(q3)/10)-A4*exp(B4*v)*IO*exp((Temp-20)*log(q4)/10))+(A52*exp(B52*v)*C2*exp((Temp-20)*log(q52)/10)-A62*exp(B62*v)*IC2*exp((Temp-20)*log(q62)/10));
ydot[41] = -(A1*exp(B1*v)*C1*exp((Temp-20)*log(q1)/10)-A2*exp(B2*v)*C2*exp((Temp-20)*log(q2)/10))-(A51*exp(B51*v)*C1*exp((Temp-20)*log(q51)/10)-A61*exp(B61*v)*IC1*exp((Temp-20)*log(q61)/10));
ydot[42] = (A1*exp(B1*v)*C1*exp((Temp-20)*log(q1)/10)-A2*exp(B2*v)*C2*exp((Temp-20)*log(q2)/10))-(A31*exp(B31*v)*C2*exp((Temp-20)*log(q31)/10)-A41*exp(B41*v)*O*exp((Temp-20)*log(q41)/10))-(A52*exp(B52*v)*C2*exp((Temp-20)*log(q52)/10)-A62*exp(B62*v)*IC2*exp((Temp-20)*log(q62)/10));
ydot[43] = (A31*exp(B31*v)*C2*exp((Temp-20)*log(q31)/10)-A41*exp(B41*v)*O*exp((Temp-20)*log(q41)/10))-(A53*exp(B53*v)*O*exp((Temp-20)*log(q53)/10)-A63*exp(B63*v)*IO*exp((Temp-20)*log(q63)/10))-(Kmax*Ku*exp(n*log(D))/(exp(n*log(D))+halfmax)*O-Ku*Obound);
ydot[44] = (A3*exp(B3*v)*IC2*exp((Temp-20)*log(q3)/10)-A4*exp(B4*v)*IO*exp((Temp-20)*log(q4)/10))+(A53*exp(B53*v)*O*exp((Temp-20)*log(q53)/10)-A63*exp(B63*v)*IO*exp((Temp-20)*log(q63)/10))-(Kmax*Ku*exp(n*log(D))/(exp(n*log(D))+halfmax)*IO-Ku*A53*exp(B53*v)*exp((Temp-20)*log(q53)/10)/(A63*exp(B63*v)*exp((Temp-20)*log(q63)/10))*IObound);
ydot[45] = (Kmax*Ku*exp(n*log(D))/(exp(n*log(D))+halfmax)*IO-Ku*A53*exp(B53*v)*exp((Temp-20)*log(q53)/10)/(A63*exp(B63*v)*exp((Temp-20)*log(q63)/10))*IObound)+(Kt/(1+exp(-(v-Vhalf)/6.789))*Cbound-Kt*IObound);
ydot[46] = (Kmax*Ku*exp(n*log(D))/(exp(n*log(D))+halfmax)*O-Ku*Obound)+(Kt/(1+exp(-(v-Vhalf)/6.789))*Cbound-Kt*Obound);
ydot[47] = -(Kt/(1+exp(-(v-Vhalf)/6.789))*Cbound-Kt*Obound)-(Kt/(1+exp(-(v-Vhalf)/6.789))*Cbound-Kt*IObound);
ydot[48] = 0;
double IKr=GKr*sqrt(ko/5.4)*O*(v-EK); 

//calculate IKs
double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
ydot[33]=(xs1ss-xs1)/txs1;
double xs2ss=xs1ss;
double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
ydot[34]=(xs2ss-xs2)/txs2;
double KsCa=1.0+0.6/(1.0+pow((3.8e-5/cai),1.4));
double GKs=0.0034/GKsfc;
if (celltype==1){
    GKs=GKs*1.4;
}
double IKs=GKs*KsCa*xs1*xs2*(v-EKs);

double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
double txk1=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
ydot[35]=(xk1ss-xk1)/txk1;
double rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
double GK1=0.1908/GK1fc;
if (celltype==1){
    GK1=GK1*1.2;
}else if (celltype==2){
    GK1=GK1*1.3;
}
double IK1=GK1*sqrt(ko)*rk1*xk1*(v-EK);


//calculate INaCa_i
double kna1=15.0;
double kna2=5.0;
double kna3=88.12;
double kasymm=12.5;
double wna=6.0e4;
double wca=6.0e4;
double wnaca=5.0e3;
double kcaon=1.5e6;
double kcaoff=5.0e3;
double qna=0.5224;
double qca=0.1670;
double hca=exp((qca*v*F)/(R*T));
double hna=exp((qna*v*F)/(R*T));
double h1=1+nai/kna3*(1+hna);
double h2=(nai*hna)/(kna3*h1);
double h3=1.0/h1;
double h4=1.0+nai/kna1*(1+nai/kna2);
double h5=nai*nai/(h4*kna1*kna2);
double h6=1.0/h4;
double h7=1.0+nao/kna3*(1.0+1.0/hna);
double h8=nao/(kna3*hna*h7);
double h9=1.0/h7;
double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
double h11=nao*nao/(h10*kna1*kna2);
double h12=1.0/h10;
double k1=h12*cao*kcaon;
double k2=kcaoff;
double k3p=h9*wca;
double k3pp=h8*wnaca;
double k3=k3p+k3pp;
double k4p=h3*wca/hca;
double k4pp=h2*wnaca;
double k4=k4p+k4pp;
double k5=kcaoff;
double k6=h6*cai*kcaon;
double k7=h5*h2*wna;
double k8=h8*h11*wna;
double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
double E1=x1/(x1+x2+x3+x4);
double E2=x2/(x1+x2+x3+x4);
double E3=x3/(x1+x2+x3+x4);
double E4=x4/(x1+x2+x3+x4);
double KmCaAct=150.0e-6;
double allo=1.0/(1.0+pow((KmCaAct/cai),2.0));
double zna=1.0;
double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
double JncxCa=E2*k2-E1*k1;
double Gncx=0.0008;
if (celltype==1){
    Gncx=Gncx*1.1;
}else if (celltype==2){
    Gncx=Gncx*1.4;
}
double INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);

//calculate INaCa_ss
h1=1+nass/kna3*(1+hna);
h2=(nass*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+nass/kna1*(1+nass/kna2);
h5=nass*nass/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+nao/kna3*(1.0+1.0/hna);
h8=nao/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
h11=nao*nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*cao*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*cass*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0/(1.0+pow((KmCaAct/cass),2.0));
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
double INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);

//calculate INaK
double k1p=949.5;
double k1m=182.4;
double k2p=687.2;
double k2m=39.4;
k3p=1899.0;
double k3m=79300.0;
k4p=639.0;
double k4m=40.0;
double Knai0=9.073;
double Knao0=27.78;
double delta=-0.1550;
double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
double Kki=0.5;
double Kko=0.3582;
double MgADP=0.05;
double MgATP=9.8;
double Kmgatp=1.698e-7;
double H=1.0e-7;
double eP=4.2;
double Khp=1.698e-7;
double Knap=224.0;
double Kxkur=292.0;
double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
double a1=(k1p*pow((nai/Knai),3.0))/(pow((1.0+nai/Knai),3.0)+pow((1.0+ki/Kki),2.0)-1.0);
double b1=k1m*MgADP;
double a2=k2p;
double b2=(k2m*pow((nao/Knao),3.0))/(pow((1.0+nao/Knao),3.0)+pow((1.0+ko/Kko),2.0)-1.0);
double a3=(k3p*pow((ko/Kko),2.0))/(pow((1.0+nao/Knao),3.0)+pow((1.0+ko/Kko),2.0)-1.0);
double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
double b4=(k4m*pow((ki/Kki),2.0))/(pow((1.0+nai/Knai),3.0)+pow((1.0+ki/Kki),2.0)-1.0);
x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
double zk=1.0;
double JnakNa=3.0*(E1*a3-E2*b3);
double JnakK=2.0*(E4*b1-E3*a1);
double Pnak=30;
if (celltype==1){
    Pnak=Pnak*0.9;
}else if (celltype==2){
    Pnak=Pnak*0.7;
}
double INaK=Pnak*(zna*JnakNa+zk*JnakK);


//calculate IKb
double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
double GKb=0.003;
if (celltype==1){
    GKb=GKb*0.6;
}
double IKb=GKb*xkb*(v-EK);

//calculate INab
double PNab=3.75e-10;
double INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

//calculate ICab
double PCab=2.5e-8;
double ICab=PCab*4.0*vffrt*(cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);

//calculate IpCa
double GpCa=0.0005;
double IpCa=GpCa*cai/(0.0005+cai);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//calculate the stimulus current, Istim
double duration=0.5;
double Istim=0.0;
if (*t<=duration){
    Istim=amp;
}else{
    Istim=0.0;
}

//update the membrane voltage
ydot[0]=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Istim);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//calculate diffusion fluxes
double JdiffNa=(nass-nai)/2.0;
double JdiffK=(kss-ki)/2.0;
double Jdiff=(cass-cai)/0.2;

//calculate ryanodione receptor calcium induced calcium release from the jsr
double bt=4.75;
double a_rel=0.5*bt;
double Jrel_inf=a_rel*(-ICaL)/(1.0+pow((1.5/cajsr),8.0));
if (celltype==2){
    Jrel_inf=Jrel_inf*1.7;
}
double tau_rel=bt/(1.0+0.0123/cajsr);

if (tau_rel<0.001){
   tau_rel=0.001;
}

ydot[36]=(Jrel_inf-Jrelnp)/tau_rel;
double btp=1.25*bt;
double a_relp=0.5*btp;
double Jrel_infp=a_relp*(-ICaL)/(1.0+pow((1.5/cajsr),8.0));
if (celltype==2){
    Jrel_infp=Jrel_infp*1.7;
}
double tau_relp=btp/(1.0+0.0123/cajsr);

if (tau_relp<0.001){
   tau_relp=0.001;
}

ydot[37]=(Jrel_infp-Jrelp)/tau_relp;
double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
double Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;

//calculate serca pump, ca uptake flux
double Jupnp=0.004375*cai/(cai+0.00092);
double Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);
if (celltype==1){
    Jupnp=Jupnp*1.3;
    Jupp=Jupp*1.3;
}
double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
double Jleak=0.0039375*cansr/15.0;
double Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

//calculate tranlocation flux
double Jtr=(cansr-cajsr)/100.0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//calcium buffer constants
double cmdnmax=0.05;
if (celltype==1){
    cmdnmax=cmdnmax*1.3;
}
double kmcmdn=0.00238;
double trpnmax=0.07;
double kmtrpn=0.0005;
double BSRmax=0.047;
double KmBSR=0.00087;
double BSLmax=1.124;
double KmBSL=0.0087;
double csqnmax=10.0;
double kmcsqn=0.8;

//update intracellular concentrations, using buffers for cai, cass, cajsr
ydot[1]=-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;
ydot[2]=-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;

ydot[3]=-(Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
ydot[4]=-(ICaK)*Acap/(F*vss)-JdiffK;

double Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow((kmcmdn+cai),2.0)+trpnmax*kmtrpn/pow((kmtrpn+cai),2.0));
ydot[5]=Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);

double Bcass=1.0/(1.0+BSRmax*KmBSR/pow((KmBSR+cass),2.0)+BSLmax*KmBSL/pow((KmBSL+cass),2.0));
ydot[6]=Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);

ydot[7]=Jup-Jtr*vjsr/vnsr;

double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow((kmcsqn+cajsr),2.0));
ydot[8]=Bcajsr*(Jtr-Jrel);

// integrate currents (qNet)
ydot[49] = INaL+ICaL+Ito+IKr+IKs+IK1;

//additional output
yout[0] = INa;
yout[1] = INaL;
yout[2] = Ito;
yout[3] = ICaL;
yout[4] = IKr;
yout[5] = IKs;
yout[6] = IK1;
yout[7] = ydot[0];
}
