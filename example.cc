/*
  example.cc

  This file contains the examples mentioned in section 3 of 
  B. Schmidt, M. Steinhauser,
  CRunDec: a C++ package for running and decoupling of the
  strong coupling and quark masses

  Compilation: 
  g++ -c CRunDec.cpp
  g++ CRunDec.o -o example example.cc 
 */

#include <iostream>
#include <cmath>
#include "CRunDec.h"

using namespace std;

int main(){

  cout << "------------------------------------------------------------" << endl;

  // ------------------------------------------------------------

  // Running of \alpha_s from Mt to Mb, 5 active flavours:

  // Pre-defined nf=5:
  CRunDec* pObjnf5 = new CRunDec(5);
  // nf=5 has to be specified in function call:
  CRunDec* pObj = new CRunDec();

  double alpha5mb1 = pObjnf5 -> AlphasExact(0.108,Mz,Mb,4);
  int nf1 = pObjnf5 -> GetNf();
  cout << "\\alpha_s(" << Mb << ") = " << alpha5mb1;
  cout << " [with " << nf1 << " active flavours]" << endl;

  double alpha5mb2 = pObj -> AlphasExact(0.108,Mz,Mb,5,4);
  cout << "\\alpha_s(" << Mb << ") = " << alpha5mb2;
  cout << " [nf=5 explicitly specified]" << endl;

  // The number of active flavours can also be set with the help of 'SetNf'
  CRunDec* pObj2 = new CRunDec();
  pObj2 -> SetNf(5);
  double alpha5mb3 = pObj -> AlphasExact(0.108,Mz,Mb,4);
  cout << "\\alpha_s(" << Mb << ") = " << alpha5mb3;
  cout << " [nf=5 specified with SetNf()]" << endl;

  cout << "------------------------------------------------------------" << endl;

  // ------------------------------------------------------------

  // Compute \alpha_s^(5)(Mz) from \alpha_s^(3)(Mtau)
  // with 1-, 2-, 3- and 4-loop accuracy.

  CRunDec crundec;

  double as;
  for(int i=1; i<=4; i++) {

    crundec.nfMmu[0].nf   = 4;
    crundec.nfMmu[0].Mth  = Mc;
    crundec.nfMmu[0].muth = 2*Mc;
    crundec.nfMmu[1].nf   = 5;
    crundec.nfMmu[1].Mth  = Mb;
    crundec.nfMmu[1].muth = Mb;

    as = crundec.AlL2AlH(0.332, Mtau, crundec.nfMmu, Mz, i);

    cout << "\\alpha_s^{(5)}(M_Z) [from \\alpha_s^{(3)}(M_tau), "<< i 
	 << " loops] = " << as << endl;
  }

  cout << "------------------------------------------------------------" << endl;

  // ------------------------------------------------------------

  // Compute mb(mb) from mb(10 GeV)=3.610 GeV and alpha_s(Mz)=0.1189

  CRunDec cmpasmb(5);
  AsmMS asmq;

  double mu,asmu,mbmu;

  mu   = 10.; 
  mbmu = 3.610;
  asmu = cmpasmb.AlphasExact(0.1189,Mz,mu,4);
  asmq.Asexact  = asmu;
  asmq.mMSexact = mbmu;
  for (;;) {
    double mbold = asmq.mMSexact;
    asmq = cmpasmb.AsmMSrunexact(mbmu, asmu, mu, mbold, 4);
    if (abs(asmq.mMSexact - mbold) < 0.00001) break;
  }

  // Alternative: use directly mMS2mSI()
  // only mb(mb) is computed
  double mbmb = cmpasmb.mMS2mSI(mbmu,asmu,mu,5,4);

  cout << "alphas(" << mu << ")       = " << asmu << endl;
  cout << "mb(" << mu << ")           = " << mbmu << endl;
  cout << "alphas(mb)       = " << asmq.Asexact << endl;
  cout << "mb(mb)           = " << asmq.mMSexact << endl;
  cout << "mb(mb) (mMS2mSI) = " << mbmb << endl;

  cout << "------------------------------------------------------------" << endl;

  // ------------------------------------------------------------

  // Compute mtMSbar from mtOS

  CRunDec cmpmt(6);
  double  
    MtOS    = 173.2;

  // Compute \alpha_s^(6)(Mt)
  cmpmt.nfMmu[0].nf   = 6;
  cmpmt.nfMmu[0].Mth  = MtOS;
  cmpmt.nfMmu[0].muth = MtOS;

  double as6Mt = cmpmt.AlL2AlH(asMz,Mz,cmpmt.nfMmu,MtOS,4);
  double mtMt  = cmpmt.mOS2mMS(MtOS,cmpmt.mq,as6Mt,MtOS,3);
  double mtmt  = cmpmt.mMS2mSI(mtMt,as6Mt,MtOS,4);

  cout << "alphas(" << MtOS << ") = " << as6Mt << endl;
  cout << "mt(" << MtOS << ")     = " << mtMt << endl;
  cout << "mt(mt)        = " << mtmt << endl;

  cout << "------------------------------------------------------------" << endl;

  // ------------------------------------------------------------

  // On-shell charm and bottom quark mass from mc(mc) and mb(mb), resp.

  //  CRunDec* pObjnf5 = new CRunDec(5);

  double alpha5mub = pObjnf5 -> AlphasExact(asMz, Mz, mub, 4);
  double mbOS2     = pObjnf5 -> mMS2mOS(mub, pObjnf5->mq, alpha5mub, mub, 2);
  double mbOS3     = pObjnf5 -> mMS2mOS(mub, pObjnf5->mq, alpha5mub, mub, 3);
  cout << "Mb (2 loops) = " << mbOS2 << endl;
  cout << "Mb (3 loops) = " << mbOS3 << endl;

  pObjnf5 -> nfMmu[0].nf   = 4;
  pObjnf5 -> nfMmu[0].Mth  = Mc;
  pObjnf5 -> nfMmu[0].muth = 2*Mc;
  double alpha4muc = pObjnf5 -> AlH2AlL(asMz, Mz, pObjnf5 -> nfMmu, mub, 4);
  double mcOS2     = pObjnf5 -> mMS2mOS(muc, pObjnf5->mq, alpha4muc, muc, 2);
  double mcOS3     = pObjnf5 -> mMS2mOS(muc, pObjnf5->mq, alpha4muc, muc, 3);
  cout << "Mc (2 loops) = " << mcOS2 << endl;
  cout << "Mc (3 loops) = " << mcOS3 << endl;

  cout << "------------------------------------------------------------" << endl;

  // ------------------------------------------------------------

  return(0);
}
