// Human Ether-à-go-go-Related Gene (hERG) channel gating model with drug
// binding kinetics
//
// Written in 2016 by Zhihua Li <zhihua.li@fda.hhs.gov>
// Modified in 2017 by Kelly Chang <kelly.chang@fda.hhs.gov>
//
// LICENSE
// To the extent possible under law, the author(s) have dedicated all
// copyright and related and neighboring rights to this software to the
// public domain worldwide. This software is distributed without any
// warranty.
//
// You should have received a copy of the CC0 Public Domain Dedication along
// with this software. If not, see
// <http://creativecommons.org/publicdomain/zero/1.0/>.
// 
// DISCLAIMER
// This code does not necessarily reflect any position of the Government or the
// Food and Drug Administration.
// 
// This software and documentation (the "Software") were developed at the Food
// and Drug Administration (FDA) by employees of the Federal Government in the
// course of their official duties. Pursuant to Title 17, Section 105 of the
// United States Code, this work is not subject to copyright protection and is
// in the public domain. Permission is hereby granted, free of charge, to any
// person obtaining a copy of the Software, to deal in the Software without
// restriction, including without limitation the rights to use, copy, modify,
// merge, publish, distribute, sublicense, or sell copies of the Software or
// derivatives, and to permit persons to whom the Software is furnished to do
// so. FDA assumes no responsibility whatsoever for use by other parties of the
// Software, its source code, documentation or compiled executables, and makes
// no guarantees, expressed or implied, about its quality, reliability, or any
// other characteristic. Further, use of this code in no way implies
// endorsement by the FDA or confers any advantage in regulatory decisions.
// Although this software can be redistributed and/or modified freely, we ask
// that any derivative works bear some notice that they are derived from it,
// and any modified versions bear some notice that they have been modified.
//
// DESCRIPTION
// C implementation of the Human Ether-à-go-go-Related Gene (hERG) channel
// gating model with drug binding kinetics.
// This code can be compiled into a dynamically linked library (DLL) for use
// in the R software environment (https://www.R-project.org/) with the
// following command:
//      R>   system("R CMD SHLIB hergmod.c")
//

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <time.h>
static double parms[51];
#define A1 parms[0]
#define B1 parms[1]
#define q1 parms[2]
#define A2 parms[3]
#define B2 parms[4]
#define q2 parms[5]
#define A3 parms[6]
#define B3 parms[7]
#define q3 parms[8]
#define A4 parms[9]
#define B4 parms[10]
#define q4 parms[11]
#define A11 parms[12]
#define B11 parms[13]
#define q11 parms[14]
#define A21 parms[15]
#define B21 parms[16]
#define q21 parms[17]
#define A31 parms[18]
#define B31 parms[19]
#define q31 parms[20]
#define A41 parms[21]
#define B41 parms[22]
#define q41 parms[23]
#define A51 parms[24]
#define B51 parms[25]
#define q51 parms[26]
#define A52 parms[27]
#define B52 parms[28]
#define q52 parms[29]
#define A53 parms[30]
#define B53 parms[31]
#define q53 parms[32]
#define A61 parms[33]
#define B61 parms[34]
#define q61 parms[35]
#define A62 parms[36]
#define B62 parms[37]
#define q62 parms[38]
#define A63 parms[39]
#define B63 parms[40]
#define q63 parms[41]
#define Kmax parms[42]
#define Ku parms[43]
#define n parms[44]
#define halfmax parms[45]
#define Kt parms[46]
#define Vhalf parms[47]
#define T parms[48]
#define timeout parms[49]
#define starttime parms[50]

void initmod(void (* odeparms)(int *, double *)){
int N=51;
odeparms(&N, parms);}

void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
if (ip[0] < 0 ) error("nout not enough!");
time_t s = time(NULL);
if((int) s - (int) starttime > timeout) error("timeout!");
ydot[0] = -(A11*exp(B11*y[10])*y[0]*exp((T-20)*log(q11)/10)-A21*exp(B21*y[10])*y[1]*exp((T-20)*log(q21)/10))+(A51*exp(B51*y[10])*y[2]*exp((T-20)*log(q51)/10)-A61*exp(B61*y[10])*y[0]*exp((T-20)*log(q61)/10));
ydot[1] = (A11*exp(B11*y[10])*y[0]*exp((T-20)*log(q11)/10)-A21*exp(B21*y[10])*y[1]*exp((T-20)*log(q21)/10))-(A3*exp(B3*y[10])*y[1]*exp((T-20)*log(q3)/10)-A4*exp(B4*y[10])*y[5]*exp((T-20)*log(q4)/10))+(A52*exp(B52*y[10])*y[3]*exp((T-20)*log(q52)/10)-A62*exp(B62*y[10])*y[1]*exp((T-20)*log(q62)/10));
ydot[2] = -(A1*exp(B1*y[10])*y[2]*exp((T-20)*log(q1)/10)-A2*exp(B2*y[10])*y[3]*exp((T-20)*log(q2)/10))-(A51*exp(B51*y[10])*y[2]*exp((T-20)*log(q51)/10)-A61*exp(B61*y[10])*y[0]*exp((T-20)*log(q61)/10));
ydot[3] = (A1*exp(B1*y[10])*y[2]*exp((T-20)*log(q1)/10)-A2*exp(B2*y[10])*y[3]*exp((T-20)*log(q2)/10))-(A31*exp(B31*y[10])*y[3]*exp((T-20)*log(q31)/10)-A41*exp(B41*y[10])*y[4]*exp((T-20)*log(q41)/10))-(A52*exp(B52*y[10])*y[3]*exp((T-20)*log(q52)/10)-A62*exp(B62*y[10])*y[1]*exp((T-20)*log(q62)/10));
ydot[4] = (A31*exp(B31*y[10])*y[3]*exp((T-20)*log(q31)/10)-A41*exp(B41*y[10])*y[4]*exp((T-20)*log(q41)/10))-(A53*exp(B53*y[10])*y[4]*exp((T-20)*log(q53)/10)-A63*exp(B63*y[10])*y[5]*exp((T-20)*log(q63)/10))-(Kmax*Ku*exp(n*log(y[9]))/(exp(n*log(y[9]))+halfmax)*y[4]-Ku*y[7]);
ydot[5] = (A3*exp(B3*y[10])*y[1]*exp((T-20)*log(q3)/10)-A4*exp(B4*y[10])*y[5]*exp((T-20)*log(q4)/10))+(A53*exp(B53*y[10])*y[4]*exp((T-20)*log(q53)/10)-A63*exp(B63*y[10])*y[5]*exp((T-20)*log(q63)/10))-(Kmax*Ku*exp(n*log(y[9]))/(exp(n*log(y[9]))+halfmax)*y[5]-Ku*A53*exp(B53*y[10])*exp((T-20)*log(q53)/10)/(A63*exp(B63*y[10])*exp((T-20)*log(q63)/10))*y[6]);
ydot[6] = (Kmax*Ku*exp(n*log(y[9]))/(exp(n*log(y[9]))+halfmax)*y[5]-Ku*A53*exp(B53*y[10])*exp((T-20)*log(q53)/10)/(A63*exp(B63*y[10])*exp((T-20)*log(q63)/10))*y[6])+(Kt/(1+exp(-(y[10]-Vhalf)/6.789))*y[8]-Kt*y[6]);
ydot[7] = (Kmax*Ku*exp(n*log(y[9]))/(exp(n*log(y[9]))+halfmax)*y[4]-Ku*y[7])+(Kt/(1+exp(-(y[10]-Vhalf)/6.789))*y[8]-Kt*y[7]);
ydot[8] = -(Kt/(1+exp(-(y[10]-Vhalf)/6.789))*y[8]-Kt*y[7])-(Kt/(1+exp(-(y[10]-Vhalf)/6.789))*y[8]-Kt*y[6]);
ydot[9] = 0;
ydot[10] = 0;
}
