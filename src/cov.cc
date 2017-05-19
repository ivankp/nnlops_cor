#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TParameter.h>
#include <TH1.h>
#include <TH2.h>

#include "timed_counter.hh"
#include "catstr.hh"
#include "binner.hh"
#include "mat.hh"

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using ivanp::cat;

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
sym_mat<double> cov(const std::vector<std::vector<double>>& h, unsigned i) {
  return { h, [i](const auto& bin){ return bin.at(i)-bin[0]; } };
}
template <typename H>
sym_mat<double> cov(const H& h, unsigned i, unsigned n) {
  auto m = cov(h,i);
  for ( ++i; i<n; ++i) m += cov(h,i);
  return m;
}

template <typename M, typename L>
TH2D* mat_hist(const M& m, const char* name, const char* title, L labels,
  std::pair<double,double> range = {0,0}
) {
  const unsigned n = m.size();
  TH2D *h = new TH2D(name,title, n,0,1, n,0,1);

  for (unsigned j=0; j<n; ++j)
    for (unsigned i=0; i<n; ++i)
      h->SetBinContent(j+1,i+1,m(i,j));

  TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis();
  for (unsigned i=1; i<=n; ++i) {
    const std::string label = labels(i);
    ax->SetBinLabel(i,label.c_str());
    ay->SetBinLabel(i,label.c_str());
  }

  if (range.first!=range.second) {
    h->SetMinimum(range.first);
    h->SetMaximum(range.second);
  }
  h->SetOption("colz");

  return h;
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

    // Fill histograms
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

  std::vector<std::vector<double>> all_bins;
  std::vector<std::string> all_bin_labels;

  TFile f(argv[1],"recreate");

  for (const auto& h : hist::all) { // loop over histograms
    cout << h.name << endl;
    const auto& axis = h->axis();

    for (unsigned bi=0; bi<axis.nbins(); ++bi) {
      auto& bin = h->bins()[bi];
      all_bins.emplace_back();
      for (unsigned wi=0; wi<total_weight.size(); ++wi) {
        auto& w = bin.w[wi];
        // scale to cross section
        w *= (xs_br_fe.GetVal() / total_weight[wi]);
        // join histograms
        all_bins.back().push_back(w);
      }
      // make bin labels
      all_bin_labels.emplace_back(cat(
        h.name," [",axis.lower(bi+1),',',axis.upper(bi+1),')'));
    }

    // Save nominal histograms
    th1(h,0);

    for (const auto* w : {
      &_w_pdf4lhc_unc, &_w_nnpdf30_unc, &_w_qcd, &_w_qcd_nnlops
    }) {
      static unsigned si = 1;

      // construct correlation matrix for this histogram
      const auto m_cor = cor(cov( *h, si, w->GetSize()+si ));

      // Save correlation matrix as a histogram
      mat_hist(m_cor.cor,
        cat("cor",w->GetBranchName()+1,'_',h.name).c_str(),
        (h.name+" correlation matrix").c_str(),
        [&axis](unsigned i){
          return cat('[',axis.lower(i),',',axis.upper(i),')');
        }, {-1,1}
      );

      if (w==&_w_qcd_nnlops) si = 1;
      else si += w->GetSize();
    }
  } // end loop over histograms

  for (const auto* w : {
    &_w_pdf4lhc_unc, &_w_nnpdf30_unc, &_w_qcd, &_w_qcd_nnlops
  }) {
    static unsigned si = 1;

    // construct the big correlation matrix
    const auto cor_all = cor(cov( all_bins, si, w->GetSize()+si ));

    mat_hist(cor_all.cor,
      cat("cor",w->GetBranchName()+1,"_all").c_str(),
      "correlation matrix",
      [&](unsigned i){ return all_bin_labels.at(i-1); }, {-1,1}
    );

    if (w==&_w_qcd_nnlops) si = 1;
    else si += w->GetSize();
  }

  f.Write();

  return 0;
}
