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
  cout << "Compute mb(mb) from mb(mb)=mub GeV and alpha_s(mb)=0.223" << endl;

  CRunDec cmpasmb(5);
  AsmMS asmq;

  double mu,asmu,mbmu;

  mu   = mub; 
  mbmu = mub;
  asmu = cmpasmb.AlphasExact(0.223,mub,mub,4);
  asmq.Asexact  = asmu;
  asmq.mMSexact = mbmu;
  for (;;) {
    double mbold = asmq.mMSexact;
    asmq = cmpasmb.AsmMSrunexact(mbmu, asmu, mu, mbold, 4);
    if (abs(asmq.mMSexact - mbold) < 0.000001) break;
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
  cout << "Compute mt(mt) from mt(on shell)" << endl;

  CRunDec cmpmt(6);
  double  
    MtOS    = Mt;

  // Compute \alpha_s^(6)(Mt)
  cmpmt.nfMmu[0].nf   = 6;
  cmpmt.nfMmu[0].Mth  = MtOS;
  cmpmt.nfMmu[0].muth = MtOS;

  double as6Mt = cmpmt.AlL2AlH(asMz,Mz,cmpmt.nfMmu,MtOS,4);
  double mtMt  = cmpmt.mOS2mMS(MtOS,cmpmt.mq,as6Mt,MtOS,4);
  double mtmt  = cmpmt.mMS2mSI(mtMt,as6Mt,MtOS,4);

  cout << "alphas(" << MtOS << ") = " << as6Mt << endl;
  cout << "mt(" << MtOS << ")     = " << mtMt << endl;
  cout << "mt(mt)         = " << mtmt << endl;

  cout << "------------------------------------------------------------" << endl;
  cout << "Compute mb(Mt) and mb(mt(mt)) from mb(mb)" << endl;

  CRunDec cmpasmbtt(5);
  AsmMS asmqbt, asmqbtt;

  // asmqbtt = cmpasmbtt.AsmMSrunexact(asmq.mMSexact, asmq.Asexact, asmq.mMSexact, mbmb, 4);
  asmqbt = cmpasmbtt.AsmMSrunexact(mbmu, asmu, mu, MtOS, 4);
  asmqbtt = cmpasmbtt.AsmMSrunexact(mbmu, asmu, mu, mtmt, 4);
  cout << "mb(" << MtOS << ")       = " << asmqbt.mMSexact << endl;
  cout << "mb(mt(mt))       = " << asmqbtt.mMSexact << endl;

  cout << "------------------------------------------------------------" << endl;
  cout << "compute mc(on shell) & mb(on shell) from mc(mc) and mb(mb)" << endl;

  CRunDec* pObjnf5 = new CRunDec(5);

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
