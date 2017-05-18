#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TParameter.h>

#include "timed_counter.hh"
#include "catstr.hh"
#include "binner.hh"
#include "mat.hh"

using std::cout;
using std::cerr;
using std::endl;

template <typename T>
const auto& get_param(TFile* f, const char* name) {
  return *dynamic_cast<const TParameter<T>*>(f->Get(name));
}

struct hist_bin {
  unsigned n = 0;
  std::vector<double> w;
  void operator++() {
    if (!n) w = weights;
    else for (size_t i=weights.size(); i; ) {
      --i;
      w[i] += weights[i];
    }
    ++n;
  }
  static std::vector<double> weights;
};
std::vector<double> hist_bin::weights;

using hist = ivanp::binner<hist_bin, std::tuple<
  ivanp::axis_spec<ivanp::container_axis<std::vector<double>>,0,0>>>;

TH1D* th1(const ivanp::named_ptr<hist>& hp, unsigned w) {
  const auto& hax = hp->axis();
  TH1D *h = new TH1D(hp.name.c_str(),"", hax.nbins(), 0, 1 );
    // hp->axis().nbins(), hp->axis().edges().data() );
  TAxis *ax = h->GetXaxis();

  unsigned n = 0;
  const auto& bins = hp->bins();
  for (unsigned i=1; i<=bins.size(); ++i) {
    const auto& b = bins[i-1];
    if (!b.n) continue;
    (*h)[i] = b.w.at(w);
    ax->SetBinLabel(i,cat('[',hax.lower(i),',',hax.upper(i),')').c_str());
    n += b.n;
  }
  h->SetEntries(n);

  return h;
}

sym_mat<double> cov(const hist& h, unsigned i) {
  return { h.bins(), [i](const hist_bin& bin){
      return bin.n ? bin.w.at(i)-bin.w[0] : 0.; // err_i
    }
  };
}

int main(int argc, char* argv[]) {
  if (argc<3) {
    cout << "usage: " << argv[0] << " out.root in1.root ..." << endl;
    return 1;
  }

  TChain chain("tree");
  for (int i=2; i<argc; ++i)
    if (!chain.Add(argv[i],0)) return 1;

  cout << '\n';
  const auto& xs_br_fe =
    get_param<double>(chain.GetFile(),"crossSectionBRfilterEff");
  cout << xs_br_fe.GetName() << " = " << xs_br_fe.GetVal() <<'\n'<< endl;

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

  hist h_pT_yy("pT_yy",{0.,20.,30.,45.,60.,80.,120.,170.,220.,350.,99999.});
  hist h_N_j_30("N_j_30",{ 0.,1.,2.,3.,9999. });
  hist h_m_jj_30("m_jj_30",{0.,170.,500.,1500.,99999.});
  hist h_Dphi_j_j_30("Dphi_j_j_30",{0.,1.0472,2.0944,3.15});
  hist h_Dphi_j_j_30_signed("Dphi_j_j_30_signed",{-3.15,-1.570796,0.,1.570796,3.15});
  hist h_pT_j1_30("pT_j1_30",{30.,55.,75.,120.,350.,99999.});

  std::vector<double> total_weight;

  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) { // LOOP
    // Read the weights ---------------------------------------------
    hist_bin::weights.clear();
    hist_bin::weights.push_back(*_w_nominal);
    for (const double w : _w_pdf4lhc_unc) hist_bin::weights.push_back(w);
    for (const double w : _w_nnpdf30_unc) hist_bin::weights.push_back(w);
    for (const double w : _w_qcd)         hist_bin::weights.push_back(w);
    for (const double w : _w_qcd_nnlops)  hist_bin::weights.push_back(w);
    // --------------------------------------------------------------

    if (total_weight.size()==0) total_weight = hist_bin::weights;
    else for (unsigned i=total_weight.size(); i; ) { --i;
      total_weight[i] += hist_bin::weights[i];
    }

    if (!*_isFiducial) continue;

    const Int_t N_j_30 = *_N_j_30;

    h_N_j_30(N_j_30);
    h_pT_yy(*_pT_yy*1e-3);
    h_pT_j1_30(*_pT_j1_30*1e-3);
    h_m_jj_30(*_m_jj_30*1e-3);
    h_Dphi_j_j_30(*_Dphi_j_j_30);
    h_Dphi_j_j_30_signed(*_Dphi_j_j_30_signed);

  } // end event loop
  cout << '\n';

  cout << "Total nominal weight: " << total_weight[0] << endl;
  cout << '\n';

  // Output =========================================================

  std::vector<std::vector<double>> bins(total_weight.size());

  TFile f(argv[1],"recreate");

  for (const auto& h : hist::all) {
    cout << "\033[0;1m" << h.name << "\033[0m" << endl;
    // scale to cross section
    for (auto& b : h->bins()) {
      for (unsigned i=total_weight.size(); i; ) { --i;
        b.w[i] *= (xs_br_fe.GetVal() / total_weight[i]);
      }
    }

    // Save nominal histograms
    th1(h,0);

    // construct covariance matrix for this histogram
    auto m_cov = cov(*h,1);
    for (unsigned i=2, n=_w_pdf4lhc_unc.GetSize()+1; i<n; ++i)
      m_cov += cov(*h,i);

    // construct correlation matrix for this histogram
    const auto m_cor = cor(m_cov);

    cout << "\nstandard deviations\n";
    std::ostream sci(cout.rdbuf());
    sci << std::scientific << std::setprecision(3);
    for (double s : m_cor.stdev) sci << s << endl;
    cout << "\ncorrelation matrix\n";
    cout << m_cor.cor << endl;

    // join histograms
    for (unsigned i=1; i<total_weight.size(); ++i)
      for (auto& b : h->bins())
        bins[i].push_back(b.w[i]-b.w[0]);
  }

  f.Write();

  return 0;
}
