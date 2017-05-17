#include <iostream>
#include <vector>

#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "timed_counter.hh"
#include "binner.hh"

template <typename... T> struct bad_type;

struct hist_bin {
  unsigned n = 0;
  std::vector<double> w;
  void operator++() {
    if (!n) w = weights;
    ++n;
    for (size_t i=weights.size(); i; ) { --i;
      w[i] += weights[i];
    }
  }
  static std::vector<double> weights;
};
std::vector<double> hist_bin::weights;

template <typename A>
using hist = ivanp::binner<hist_bin, std::tuple<
  ivanp::axis_spec<ivanp::container_axis<std::vector<A>>>>>;

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char* argv[]) {
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
  TTreeReaderArray<double> _w_pdf4lhc_unc(reader,"w_pdf4lhc_unc");
  TTreeReaderArray<double> _w_nnpdf30_unc(reader,"w_nnpdf30_unc");
  TTreeReaderArray<double> _w_qcd(reader,"w_qcd");
  TTreeReaderArray<double> _w_qcd_nnlops(reader,"w_qcd_nnlops");

  hist<Int_t>    h_N_j_30({ 0,1,2,3,4 });
  hist<Double_t> h_pT_yy({ 0.,20.,30.,45.,60.,80.,120.,170.,220.,350. });
  hist<Double_t> h_pT_j1_30({ 30.,40.,55.,75.,95.,120.,170.,400. });
  hist<Double_t> h_m_jj_30({ 0.,200.,500.,1000. });
  hist<Double_t> h_Dphi_j_j_30({ 0.,1.0472,2.0944,3.1416 });
  hist<Double_t> h_Dphi_j_j_30_signed({ -3.1416,-1.5708,0.,1.5708,3.1416 });

  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    if (!*_isFiducial) continue;

    // Read the weights
    hist_bin::weights.clear();
    hist_bin::weights.push_back(*_w_nominal);
    for (const double w : _w_pdf4lhc_unc) hist_bin::weights.push_back(w);
    for (const double w : _w_nnpdf30_unc) hist_bin::weights.push_back(w);
    for (const double w : _w_qcd)         hist_bin::weights.push_back(w);
    for (const double w : _w_qcd_nnlops)  hist_bin::weights.push_back(w);

    const Int_t N_j_30 = *_N_j_30;

    h_N_j_30(N_j_30);
    h_pT_yy(*_pT_yy*1e-3);

    if (N_j_30 < 1) continue; // 1 jet ==============================

    h_pT_j1_30(*_pT_j1_30*1e-3);

    if (N_j_30 < 2) continue; // 2 jets =============================

    h_m_jj_30(*_m_jj_30*1e-3);
    h_Dphi_j_j_30(*_Dphi_j_j_30);
    h_Dphi_j_j_30_signed(*_Dphi_j_j_30_signed);
  }

  TFile f(argv[1],"recreate");

  TH1D *out_h_pT_yy = new TH1D("pT_yy","",
    h_pT_yy.axis().nbins(), h_pT_yy.axis().edges().data() );

  unsigned n = 0;
  for (const auto& b : h_pT_yy.bins()) {
    static int i = 0;
    if (!b.n) continue;
    (*out_h_pT_yy)[i] = b.w.at(0);
    n += b.n;
    ++i;
  }
  out_h_pT_yy->SetEntries(n);

  f.Write();

  return 0;
}
