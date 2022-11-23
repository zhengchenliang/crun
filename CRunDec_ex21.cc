/*
  CRunDec_ex21.cc

  This file contains examples illustrating the new features
  of RunDec2.1.

  Compilation: 
  g++ -std=c++11 -o CRunDec_ex21 CRunDec_ex21.cc CRunDec.cpp
*/

#include "CRunDec.h"

#include <iostream>

#define TOP
//#define BOTTOM

using namespace std;

// Uncertainty set to 0.2%

        double xerr = 0.002;

// Numerical input adopted from arXiv:1502.01030

	CRunDec crd(5);
	#ifdef TOP
	double mtOSnum   = 173.34;
	double mtPSnum   = 171.792;
	double mtMSnum   = 163.643;
	double mt1Snum   = 172.227;
	double mtRSnum   = 171.215;
	double asMznum   = 0.1185;
	double Mznum     = 91.1876;
	double mufnum    = 20.0;
	int nlnum     = 5;
	int nloopsmin = 1;
	int nloopsnum = 4;
	double as6(double mu) {
		crd.nfMmu[0].nf   = 6;
		crd.nfMmu[0].Mth  = mtOSnum;
		crd.nfMmu[0].muth = 2*mtOSnum;
		return crd.AlL2AlH(asMznum,Mznum,crd.nfMmu,mu,nloopsnum);
	}
	double as5(double mu) {
		return crd.AlphasExact(asMznum,Mznum, mu, 5, 4);
	}
	#endif
	#ifdef BOTTOM
	double mtOSnum   = 4.78;
	double mtPSnum   = 4.483;
	double mtMSnum   = 4.163;
	double mt1Snum   = 4.670;
	double mtRSnum   = 4.365;
	double asMznum   = 0.1185;
	double Mznum     = 91.1876;
	double mufnum    = 2.0;
	int nlnum     = 4;
	int nloopsmin = 1;
	int nloopsnum = 4;
	double as6(double mu) {
		return crd.AlphasExact(asMznum,Mznum, mu, 5, 4);
	}
	double as5(double mu) {
    	double mbMSin = mtMSnum;
    	double asMzin = asMznum;
		return crd.AlphasExact(crd.DecAsDownSI(crd.AlphasExact(asMzin,Mznum, 2*mbMSin, 5, 4), mbMSin, 2*mbMSin, 4, 4), 2*mbMSin, mu, 4, 4);
	}
	#endif

int main() {

        cout << "Mass relations:" << endl;
  
	cout << "mOS2mMS: " << endl;
	for(int i = 0; i <= nloopsnum; i++) {
		if(i != 0)
			cout << crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, i) - 
					crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, i - 1) << " ";
		else
			cout << crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, i) << " ";
	}
	cout << crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, 4, 1.+xerr) - crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, 4);
	cout << endl << "Summed up: " << crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, 4) << endl;
	cout << endl;
	cout << "mMS2mOS" << endl;
	for(int i = 0; i <= nloopsnum; i++) {
		if(i != 0)
			cout << crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, i) - 
					crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, i - 1) << " ";
		else
			cout << crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, i) << " ";
	}
	cout << crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, 4, 1.+xerr) - crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, 4);
	cout << endl <<"Summed up: " << crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, 4) << endl;
	cout << endl;
	cout << "mOS2mSI" << endl;

	for(int i = 0; i <= nloopsnum; i++) {
		if(i != 0)
			cout << crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, i) - 
					crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, i - 1) << " ";
		else
			cout << crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, i) << " ";
	}
	cout << crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, 4, 1.+xerr) - crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, 4);
	cout << endl << "Summed up: " << crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, 4);
	//cout << endl << "+"<<100*xerr<<"% Error: " << crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, 4, 1-xerr) << endl;
	cout << endl;
   
   cout << endl; 
   cout << "{ctr. value, 4-loop coef. + "<<100*xerr<<"\%, 4-loop coef. - "<<100*xerr<<"\%}" << endl;
   cout << "mOS2mMS: {" << crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, 4) << ", "
        << crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, 4, 1.+xerr) << ", "
        << crd.mOS2mMS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum+1, 4, 1-xerr) << "}" << endl;
   cout << "mMS2mOS: {" << crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, 4) << ", "
        << crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, 4, 1.+xerr) << ", "
        << crd.mMS2mOS(mtMSnum, crd.mq, as6(mtMSnum), mtMSnum, nlnum+1, 4, 1-xerr) << "}" << endl;
   cout << "mOS2mSI: {" << crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, 4) << ", "
        << crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, 4, 1.+xerr) << ", "
        << crd.mOS2mSI(mtOSnum, crd.mq, as6(mtOSnum), nlnum+1, 4, 1-xerr) << "}" << endl;
   cout << endl;

	// ------------------------------------------------------------

	for(int i = nloopsmin; i <= nloopsnum; i++) {
	cout << "iloops = " << i << ", mMS2mPS = " << crd.mMS2mPS(mtMSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mOS2mPS = " << crd.mOS2mPS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mPS2mMS = " << crd.mPS2mMS(mtPSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mPS2mSI = " << crd.mPS2mSI(mtPSnum, crd.mq, &as5, mufnum, nlnum, i) << endl;
	}
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mMS2mPS = " << crd.mMS2mPS(mtMSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mPS2mMS = " << crd.mPS2mMS(mtPSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mPS2mSI = " << crd.mPS2mSI(mtPSnum, crd.mq, &as5, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << endl;


	for(int i = nloopsmin; i <= nloopsnum; i++) {
	cout << "iloops = " << i << ", mMS2m1S = " << crd.mMS2m1S(mtMSnum, crd.mq, as5(mtMSnum), mtMSnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mOS2m1S = " << crd.mOS2m1S(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", m1S2mMS = " << crd.m1S2mMS(mt1Snum, crd.mq, as5(mtMSnum), mtMSnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", m1S2mSI = " << crd.m1S2mSI(mt1Snum, crd.mq, &as5, nlnum, i) << endl;
	}
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mMS2m1S = " << crd.mMS2m1S(mtMSnum, crd.mq, as5(mtMSnum), mtMSnum, nlnum, 4, 1.+xerr) << endl;
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, m1S2mMS = " << crd.m1S2mMS(mt1Snum, crd.mq, as5(mtMSnum), mtMSnum, nlnum, 4, 1.+xerr) << endl;
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, m1S2mSI = " << crd.m1S2mSI(mt1Snum, crd.mq, &as5, nlnum, 4, 1.+xerr) << endl;
	cout << endl;

	for(int i = nloopsmin; i <= nloopsnum; i++) {
	cout << "iloops = " << i << ", mMS2mRS = " << crd.mMS2mRS(mtMSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mOS2mRS = " << crd.mOS2mRS(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mRS2mMS = " << crd.mRS2mMS(mtRSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mRS2mSI = " << crd.mRS2mSI(mtRSnum, crd.mq, &as5, mufnum, nlnum, i) << endl;
	}
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mMS2mRS = " << crd.mMS2mRS(mtMSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mRS2mMS = " << crd.mRS2mMS(mtRSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mRS2mSI = " << crd.mRS2mSI(mtRSnum, crd.mq, &as5, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << endl;

	for(int i = nloopsmin; i <= nloopsnum; i++) {
	cout << "iloops = " << i << ", mMS2mRSp = " << crd.mMS2mRSp(mtMSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mOS2mRSp = " << crd.mOS2mRSp(mtOSnum, crd.mq, as6(mtOSnum), mtOSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mRSp2mMS = " << crd.mRSp2mMS(mtRSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, i) << endl;
	cout << "iloops = " << i << ", mRSp2mSI = " << crd.mRSp2mSI(mtRSnum, crd.mq, &as5, mufnum, nlnum, i) << endl;
	}
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mMS2mRSp = " << crd.mMS2mRSp(mtMSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mRSp2mMS = " << crd.mRSp2mMS(mtRSnum, crd.mq, as5(mtMSnum), mtMSnum, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << "iloops = " << 4 << ", Error = "<<100*xerr<<"%, mRSp2mSI = " << crd.mRSp2mSI(mtRSnum, crd.mq, &as5, mufnum, nlnum, 4, 1.+xerr) << endl;
	cout << endl;

	// ------------------------------------------------------------

	// running to 5 loops:

	cout << "Running" << endl;
	double muend = 2500.;
	double mqmu0 = 2.5;
	int nf = 5;
	for(int iloops = 1; iloops<=5; iloops++) {
	double asmu = crd.AlphasExact(asMznum, Mznum, muend, nf, iloops);
  	cout << "iloops = " << iloops << ", as[" << muend << "] = " << asmu << endl;
	cout << "iloops = " << iloops << ", mq[" << muend << "] = " << crd.mMS2mMS(mqmu0, asMznum, asmu, nf, iloops) << endl;
	}

	// ------------------------------------------------------------

	// decoupling to 4 loops:
        // note: the argument 'iloops=5' 
        // refers to 5-loop running and 4-loop decoupling

	cout << "Decoupling" << endl;
	double muth   = 2500.;
	double massth = 150.;
	double alphas = 0.101;
	double massq  = 2.5;
	int nl = 5;
	for(int iloops = 1; iloops<=5; iloops++) {
  		cout << "iloops = " << iloops << ", DecAsUpOS   = " <<  crd.DecAsUpOS(alphas  ,massth ,muth ,nl ,iloops ) << endl;
  		cout << "iloops = " << iloops << ", DecAsDownOS = " <<  crd.DecAsDownOS(alphas  ,massth ,muth ,nl ,iloops ) << endl;
  		cout << "iloops = " << iloops << ", DecMqUpOS   = " <<  crd.DecMqUpOS(massq,alphas  ,massth ,muth ,nl ,iloops ) << endl;
  		cout << "iloops = " << iloops << ", DecMqDownOS = " <<  crd.DecMqDownOS(massq,alphas ,massth ,muth ,nl ,iloops ) << endl;
	}

	cout << endl;

}