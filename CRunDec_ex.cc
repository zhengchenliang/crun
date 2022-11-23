/*
  CRunDec_ex.cc

  This file contains all examples from the appendix of
  K.G. Chetyrkin, J.H. Kuhn and M.Steinhauser,
  ``RunDec: A Mathematica package for running and decoupling of the strong
  coupling and quark masses,''
  Comput.\ Phys.\ Commun.\  {\bf 133} (2000) 43
  [arXiv:hep-ph/0004189].

  Compilation: 
  g++ -c CRunDec.cpp
  g++ CRunDec.o -o CRunDec_ex CRunDec_ex.cc 
*/

#include <iostream>
#include "CRunDec.h"

using namespace std;

int main(){

  CRunDec* pObj = new CRunDec();

  // ----------

  double lamexpl = pObj -> LamExpl(asMz,Mz,5,3);
  double lamimpl = pObj -> LamImpl(asMz,Mz,5,3);
  double alslam  = pObj -> AlphasLam(0.208,Mb,5,3);
  double alsexa  = pObj -> AlphasExact(asMz,Mz,Mb,5,3);

  cout << "Running of \\alpha_s:" << endl;
  cout << "\\Lambda_expl = " << lamexpl << endl;
  cout << "\\Lambda_impl = " << lamimpl << endl;
  cout << "\\alpha_s^Lam = " << alslam  << endl;
  cout << "\\alpha_s^exa = " << alsexa  << endl << endl;

  // ----------

  double mms1    = pObj -> mOS2mMS(Mt,pObj->mq,0.107,Mt,6,3);
  double mos1    = pObj -> mMS2mOS(165,pObj->mq,0.107,Mt,6,3);
  double mms1run = pObj -> mOS2mMSrun(Mt,pObj->mq,0.107,Mt,6,3);
  double mos1run = pObj -> mMS2mOSrun(165,pObj->mq,0.107,Mt,6,3);
  double mms1it  = pObj -> mOS2mMSit(Mt,pObj->mq,0.107,Mt,6,3);

  cout << "Relations among quark masses:" << endl;
  cout << "mMS    = " << mms1    << endl;
  cout << "mOS    = " << mos1    << endl;
  cout << "mMSrun = " << mms1run << endl;
  cout << "mOSrun = " << mos1run << endl;
  cout << "mMSit  = " << mms1it  << endl << endl;

  pObj -> mq[0]=1.6;
  double msi1  = pObj -> mOS2mSI(Mb,pObj->mq,0.217,5,3);
  double mms2  = pObj -> mMS2mMS(3.85,0.217,asMz,5,4);
  double msi2  = pObj -> mMS2mSI(3.85,0.217,4.7,5,3);
  double mri1  = pObj -> mMS2mRI(2.695,asMz,5,3);
  double mms3  = pObj -> mRI2mMS(Mt,0.107,6,3);
  double mrgi1 = pObj -> mMS2mRGI(2.69,asMz,5,4);
  double mms4  = pObj -> mRGI2mMS(14.25,asMz,5,4);

  cout << "mSI(from OS)  = " << msi1 << endl;
  cout << "mMS(from MS)  = " << mms2 << endl;
  cout << "mSI(from MS)  = " << msi2 << endl;
  cout << "mRI(from MS)  = " << mri1 << endl;
  cout << "mMS(from RI)  = " << mms3 << endl;
  cout << "mRGI(from MS) = " << mrgi1 << endl;
  cout << "mMS(from RGI) = " << mms4 << endl << endl;

  // ----------

  double as6 = pObj -> DecAsUpOS(asMz,Mt,Mz,5,4);
  double as5 = pObj -> DecAsDownOS(0.105,Mt,200,6,4);
  double mb6 = pObj -> DecMqUpOS(2.7,asMz,Mt,Mz,5,4);
  double mc4 = pObj -> DecMqDownOS(0.58,asMz,4.7,Mz,5,4);

  cout << "Decoupling:" << endl;
  cout << "as6 = " << as6 << endl;
  cout << "as5 = " << as5 << endl;
  cout << "mb6 = " << mb6 << endl;
  cout << "mc4 = " << mc4 << endl << endl;

  // ----------

  pObj->nfMmu[0].nf=5;
  pObj->nfMmu[0].Mth=Mb;
  pObj->nfMmu[0].muth=5;
  pObj->nfMmu[1].nf=6;
  pObj->nfMmu[1].Mth=Mt;
  pObj->nfMmu[1].muth=200;
  double as62 = pObj -> AlL2AlH(0.338,1.6,pObj->nfMmu,500,4);

  pObj->nfMmu[0].nf=6;
  pObj->nfMmu[0].Mth=Mt;
  pObj->nfMmu[0].muth=200;
  pObj->nfMmu[1].nf=5;
  pObj->nfMmu[1].Mth=Mb;
  pObj->nfMmu[1].muth=5;
  double as4 = pObj -> AlH2AlL(0.0952,500,pObj->nfMmu,1.6,4);

  pObj->nfMmu[0].nf=5;
  pObj->nfMmu[0].Mth=Mb;
  pObj->nfMmu[0].muth=5;
  double mc42 = pObj -> mL2mH(1.2,0.403,1.2,pObj->nfMmu,Mz,4);

  pObj->nfMmu[0].nf=5;
  pObj->nfMmu[0].Mth=Mb;
  pObj->nfMmu[0].muth=5;
  double mc43 = pObj -> mH2mL(0.580,asMz,Mz,pObj->nfMmu,1.2,4);

  cout << "Running and decoupling:" << endl;
  cout << "as6(500) = " << as62 << endl;
  cout << "as4(1.6) = " << as4 << endl;
  cout << "mc4(Mz)  = " << mc42 << endl;
  cout << "mc4(1.2) = " << mc43 << endl;

  return(0);
}
