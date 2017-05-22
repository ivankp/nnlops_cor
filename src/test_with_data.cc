#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

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

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#include "mat.hh"

using std::cout;
using std::cerr;
using std::endl;
using ivanp::cat;

template <typename T>
auto& get_param(TFile* f, const char* name) {
  return *dynamic_cast<TParameter<T>*>(f->Get(name));
}

struct hist_bin {
  unsigned n = 0;
  std::vector<double> w;
  void operator++() {
    if (!n) w = weights;
    else for (unsigned i=weights.size(); i; ) {
      --i;
      w[i] += weights[i];
    }
    ++n;
  }
  static std::vector<double> weights, total_weights;
};
std::vector<double> hist_bin::weights, hist_bin::total_weights;

using hist = ivanp::binner<hist_bin, std::tuple<
  ivanp::axis_spec<ivanp::container_axis<std::vector<double>>,0,0>>>;

template <typename A>
std::string bin_label(const A& axis, unsigned i) {
  return cat('[',axis.lower(i),',',axis.upper(i),')');
}

TH1D* th1(const ivanp::named_ptr<hist>& hp, unsigned w, const char* title="") {
  const auto& hax = hp->axis();
  TH1D *h = new TH1D(hp.name.c_str(),title, hax.nbins(), 0, 1 );
  TAxis *ax = h->GetXaxis();

  unsigned n = 0;
  const auto& bins = hp->bins();
  for (unsigned i=1; i<=bins.size(); ++i) {
    const auto& b = bins[i-1];
    if (!b.n) continue;
    (*h)[i] = b.w.at(w);
    ax->SetBinLabel(i,bin_label(hax,i).c_str());
    n += b.n;
  }
  h->SetEntries(n);

  return h;
}

template <typename L>
TH1D* th1(const std::vector<double>& v,
  const char* name, const char* title, L labels
) {
  TH1D *h = new TH1D(name,title,v.size(),0,1);
  TAxis *ax = h->GetXaxis();

  for (unsigned i=1; i<=v.size(); ++i) {
    (*h)[i] = v[i-1];
    ax->SetBinLabel(i,labels(i).c_str());
  }
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
sym_mat<double> cov(const H& h, unsigned i, unsigned end) {
  auto m = cov(h,i);
  for ( ++i; i<end; ++i) m += cov(h,i);
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
  if (argc<2) {
    cout << "usage: " << argv[0] << " out.root in1.root ..." << endl;
    return 1;
  }

  TChain chain("tree");
  for (int i=1; i<argc; ++i)
    if (!chain.Add(argv[i],0)) return 1;

  cout << '\n';
  auto& xs_br_fe =
    get_param<double>(chain.GetFile(),"crossSectionBRfilterEff");
  if (xs_br_fe.GetVal() < 0.) xs_br_fe.SetVal(0.1101404);
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
  std::array<TTreeReaderArray<double>,1> _weights {{
    {reader,"w_pdf4lhc_unc"}
    // {reader,"w_nnpdf30_unc"},
    // {reader,"w_qcd"}
    // {reader,"w_qcd_nnlops"}
  }};

  hist h_pT_yy("pT_yy",{0.,20.,30.,45.,60.,80.,120.,170.,220.,350.,99999.});
  hist h_N_j_30("N_j_30",{ 0.,1.,2.,3.,9999. });
  // hist h_m_jj_30("m_jj_30",{0.,170.,500.,1500.,99999.});
  // hist h_Dphi_j_j_30("Dphi_j_j_30",{0.,1.0472,2.0944,3.15});
  // hist h_Dphi_j_j_30_signed("Dphi_j_j_30_signed",{-3.15,-1.570796,0.,1.570796,3.15});
  // hist h_pT_j1_30("pT_j1_30",{30.,55.,75.,120.,350.,99999.});

  unsigned n_bad_weights = 0, n_events_with_bad_weights = 0,
           n_fiducial = 0;

  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) { // LOOP
    // if (ent==100000) break; // TEST
    // Read the weights ---------------------------------------------
    hist_bin::weights.clear();
    hist_bin::weights.push_back(*_w_nominal);
    for (auto& _w : _weights)
      for (const double w : _w)
        hist_bin::weights.push_back(w);

    // handle bad weights
    bool had_bad_weight = false;
    for (double& w : hist_bin::weights) {
      if (!std::isfinite(w) || std::abs(w) > 150.) {
        // cout << "w = " << w << endl;
        w = 1.;
        ++n_bad_weights;
        if (!had_bad_weight) had_bad_weight = true;
      }
    }
    if (had_bad_weight) ++n_events_with_bad_weights;

    // add the event's weights to total
    if (hist_bin::total_weights.size()==0) {
      hist_bin::total_weights = hist_bin::weights;
    } else for (unsigned i=hist_bin::weights.size(); i; ) { --i;
      hist_bin::total_weights[i] += hist_bin::weights[i];
    }
    // --------------------------------------------------------------

    // TEST( hist_bin::weights[0] )

    if (!*_isFiducial) continue;
    ++n_fiducial;

    // Fill histograms
    h_N_j_30(*_N_j_30);
    h_pT_yy(*_pT_yy*1e-3);
    // h_pT_j1_30(*_pT_j1_30*1e-3);
    // h_m_jj_30(*_m_jj_30*1e-3);
    // h_Dphi_j_j_30(*_Dphi_j_j_30);
    // h_Dphi_j_j_30_signed(*_Dphi_j_j_30_signed);

    // if (n_fiducial==5) break; // TEST

  } // end event loop
  cout << '\n';

  cout << "N fiducial events: " << n_fiducial << endl;
  cout << "Total nominal weight: " << hist_bin::total_weights[0] << endl;
  cout << '\n';

  TEST( n_bad_weights )
  TEST( n_events_with_bad_weights )

  cout << '\n';

  // Output =========================================================

  // for (int wi : {0,1,2})
  //   TEST( (xs_br_fe.GetVal() / hist_bin::total_weights[wi]) );

  // for (auto w : hist_bin::total_weights)
  //   cout << w << endl;
  // cout << endl;

  // TEST( _weights[0].GetSize() )
  // TEST( _weights[1].GetSize() )

  for (const auto& h : hist::all) { // loop over histograms
    cout << h.name << endl;
    const auto& axis = h->axis();

    for (unsigned bi=0; bi<axis.nbins(); ++bi) {
      auto& bin = h->bins()[bi];
      if (bin.w.size()==0) // assign zeros to empty bin
        bin.w.assign(hist_bin::weights.size(),0.);
      for (unsigned wi=0; wi<hist_bin::weights.size(); ++wi) {
        auto& w = bin.w[wi];
        // scale to cross section
        // if (wi<2) cout << ' ' << w;
        // const double f = (xs_br_fe.GetVal() / hist_bin::total_weights[wi]);
        // TEST( f )
        // w *= f;
        // w *= (xs_br_fe.GetVal() / hist_bin::total_weights[wi]);
        w *= (0.1101404*36100.0 / hist_bin::total_weights[wi]);
        // if (wi<2) cout << ' ' << w;
      }
      // cout << endl;
    }
    cout << endl;

    TEST( h->bins()[0].w.size() )
    for (const auto& b : h->bins()) {
      for (auto w : b.w) cout << ' ' << w;
      cout << endl;
    }
    cout << endl;

    unsigned i1 = 1;
    for (auto& w : _weights) {
      cout << w.GetBranchName() << endl;
      // construct correlation matrix for this histogram
      const auto m_cov = cov( *h, i1, w.GetSize()+i1 );
      // const auto m_cov = cov( *h, 1, 30+1 );
      cout << m_cov << endl;
      const auto m_cor = cor(m_cov);
      cout << m_cor.cor << endl;

      i1 += w.GetSize();
      break;
    }
  } // end loop over histograms

  return 0;
}
