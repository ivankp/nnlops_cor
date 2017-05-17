#include <iostream>

#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "timed_counter.hh"

using std::cout;
using std::cerr;
using std::endl;

const double b_N_j_30[] = { 0,1,2,3,4 };
const double b_pT_yy[] = { 0,20,30,45,60,80,120,170,220,350 };
const double b_pT_j1_30[] = { 30,40,55,75,95,120,170,400 };
const double b_m_jj_30[] = { 0,200,500,1000 };
const double b_Dphi_j_j_30[] = { 0,1.0472,2.0944,3.1416 };
const double b_Dphi_j_j_30_signed[] = { -3.1416,-1.5708,0,1.5708,3.1416 };

#define HIST(NAME) \
  TH1D *h_##NAME = new TH1D(#NAME,"", \
    std::extent<decltype(b_##NAME)>::value-1,b_##NAME);

int main(int argc, char* argv[]) {
  if (argc<3) {
    cout << "usage: " << argv[0] << " out.root in1.root ..." << endl;
    return 1;
  }

  TChain chain("tree");
  for (int i=2; i<argc; ++i)
    if (!chain.Add(argv[i],0)) return 1;

  TTreeReader reader(&chain);

  TTreeReaderValue<Bool_t>  _isFiducial(reader,"isFiducial");
  TTreeReaderValue<Int_t>   _N_j_30(reader,"N_j_30");
  TTreeReaderValue<Float_t> _pT_yy(reader,"pT_yy");
  TTreeReaderValue<Float_t> _pT_j1_30(reader,"pT_j1_30");
  TTreeReaderValue<Float_t> _m_jj_30(reader,"m_jj_30");
  TTreeReaderValue<Float_t> _Dphi_j_j_30(reader,"Dphi_j_j_30");
  TTreeReaderValue<Float_t> _Dphi_j_j_30_signed(reader,"Dphi_j_j_30_signed");

  TTreeReaderValue<Double_t> _w_nominal(reader,"w_nominal");

  TFile f(argv[1],"recreate");
  HIST(N_j_30)
  HIST(pT_yy)
  HIST(pT_j1_30)
  HIST(m_jj_30)
  HIST(Dphi_j_j_30)
  HIST(Dphi_j_j_30_signed)

  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    if (!*_isFiducial) continue;

    const double w_nominal = *_w_nominal;
    const Int_t N_j_30 = *_N_j_30;

    h_N_j_30->Fill(N_j_30,w_nominal);
    h_pT_yy->Fill(*_pT_yy*1e-3,w_nominal);

    if (N_j_30 < 1) continue; // 1 jet ==============================

    h_pT_j1_30->Fill(*_pT_j1_30*1e-3,w_nominal);

    if (N_j_30 < 2) continue; // 2 jets =============================

    h_m_jj_30->Fill(*_m_jj_30*1e-3,w_nominal);
    h_Dphi_j_j_30->Fill(*_Dphi_j_j_30,w_nominal);
    h_Dphi_j_j_30_signed->Fill(*_Dphi_j_j_30_signed,w_nominal);

  }

  f.Write();

  return 0;
}
