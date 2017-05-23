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
#include "mat.hh"

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using ivanp::cat;

template <typename T>
auto& get_param(TFile* f, const char* name) {
  return *dynamic_cast<TParameter<T>*>(f->Get(name));
}

struct hist_bin {
  unsigned n = 0;
  std::vector<double> w;
  static std::vector<double> weights, total_weights;

  void operator++() {
    if (!n) w = weights;
    else for (size_t i=weights.size(); i; ) {
      --i;
      w[i] += weights[i];
    }
    ++n;
  }
  double min(unsigned first, unsigned last) const {
    return *std::min_element(w.begin()+first,w.begin()+last);
  }
  double max(unsigned first, unsigned last) const {
    return *std::max_element(w.begin()+first,w.begin()+last);
  }

  template <typename T>
  auto operator[](T&& f) const -> decltype(f(*this)) { return f(*this); }
  template <typename T>
  inline std::enable_if_t<std::is_integral<T>::value,double>
  operator[](T i) const { return w[i]; }
  template <typename T>
  inline std::enable_if_t<std::is_integral<T>::value,double>&
  operator[](T i) { return w[i]; }
};
std::vector<double> hist_bin::weights, hist_bin::total_weights;

using hist = ivanp::binner<hist_bin, std::tuple<
  ivanp::axis_spec<ivanp::container_axis<std::vector<double>>,0,0>>>;

template <typename A>
std::string bin_label(const A& axis, unsigned i) {
  return cat('[',axis.lower(i),',',axis.upper(i),')');
}

template <typename WeightPred>
TH1D* th1(const ivanp::named_ptr<hist>& hp,
  WeightPred predicate,
  const std::string& name, const std::string& title
) {
  const auto& hax = hp->axis();
  TH1D *h = new TH1D(name.c_str(), title.c_str(), hax.nbins(), 0, hax.nbins());
  TAxis *ax = h->GetXaxis();

  unsigned n = 0;
  const auto& bins = hp->bins();
  for (unsigned i=1; i<=bins.size(); ++i) {
    const auto& b = bins[i-1];
    if (!b.n) continue;
    (*h)[i] = b[predicate];
    ax->SetBinLabel(i,bin_label(hax,i).c_str());
    n += b.n;
  }
  h->SetEntries(n);

  return h;
}

template <typename L>
TH1D* th1(const std::vector<double>& v,
  const std::string& name, const std::string& title, L labels
) {
  TH1D *h = new TH1D(name.c_str(),title.c_str(),v.size(),0,v.size());
  TAxis *ax = h->GetXaxis();

  for (unsigned i=1; i<=v.size(); ++i) {
    (*h)[i] = v[i-1];
    ax->SetBinLabel(i,labels(i).c_str());
  }
  return h;
}

sym_mat<double> cov(const hist& h, unsigned i) {
  return { h.bins(), [i](const hist_bin& bin){
      return bin.n ? bin[i]-bin[0] : 0.; // err_i
    }
  };
}
sym_mat<double> cov(const std::vector<std::vector<double>>& h, unsigned i) {
  return { h, [i](const auto& bin){ return bin[i]-bin[0]; } };
}
template <typename H>
sym_mat<double> cov(const H& h, unsigned i, unsigned n) {
  auto m = cov(h,i);
  for ( ++i; i<n; ++i) m += cov(h,i);
  return m;
}

template <typename M, typename L>
TH2D* mat_hist(const M& m,
  const std::string& name, const std::string& title, L labels,
  std::pair<double,double> range = {0,0}
) {
  const unsigned n = m.size();
  TH2D *h = new TH2D(name.c_str(),title.c_str(), n,0,n, n,0,n);

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
  h->SetOption("COLZTEXT");

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
  std::array<TTreeReaderArray<double>,4> _weights {{
    {reader,"w_pdf4lhc_unc"},
    {reader,"w_nnpdf30_unc"},
    {reader,"w_qcd"},
    {reader,"w_qcd_nnlops"}
  }};

  hist h_pT_yy("pT_yy",{0.,20.,30.,45.,60.,80.,120.,170.,220.,350.});
  hist h_N_j_30("N_j_30",{ 0.,1.,2.,3.,9999. });
  hist h_m_jj_30("m_jj_30",{0.,170.,500.,1500.});
  hist h_Dphi_j_j_30("Dphi_j_j_30",{0.,1.0472,2.0944,3.15});
  hist h_Dphi_j_j_30_signed("Dphi_j_j_30_signed",{-3.15,-1.570796,0.,1.570796,3.15});
  hist h_pT_j1_30("pT_j1_30",{30.,55.,75.,120.,350.});

  unsigned n_fiducial = 0;

  struct {
    unsigned weights = 0, events = 0, nominal = 0;
  } n_inf, n_large;

  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) { // LOOP
    // Read the weights ---------------------------------------------
    hist_bin::weights.clear(); // current event weights
    hist_bin::weights.push_back(*_w_nominal);
    for (auto& _w : _weights)
      for (const double w : _w)
        hist_bin::weights.push_back(w);

    // handle bad weights
    bool first_bad_weight = true, is_nom_w = true;
    for (double& w : hist_bin::weights) {
      bool is_inf = false;
      if ((is_inf = !std::isfinite(w)) || std::abs(w) > 150.) {
        w = 1.;
        if (is_inf) { // count bad weights
          ++n_inf.weights;
          if (first_bad_weight) ++n_inf.events;
          if (is_nom_w) ++n_inf.nominal;
        } else {
          ++n_large.weights;
          if (first_bad_weight) ++n_large.events;
          if (is_nom_w) ++n_large.nominal;
        }
        first_bad_weight = false;
      }
      is_nom_w = false;
    }

    // add the event's weights to total
    if (hist_bin::total_weights.size()==0) {
      hist_bin::total_weights = hist_bin::weights;
    } else for (unsigned i=hist_bin::weights.size(); i; ) { --i;
      hist_bin::total_weights[i] += hist_bin::weights[i];
    }
    // --------------------------------------------------------------

    if (!*_isFiducial) continue;
    ++n_fiducial;

    // Fill histograms
    h_N_j_30(*_N_j_30);
    h_pT_yy(*_pT_yy*1e-3);
    h_pT_j1_30(*_pT_j1_30*1e-3);
    h_m_jj_30(*_m_jj_30*1e-3);
    h_Dphi_j_j_30(*_Dphi_j_j_30);
    h_Dphi_j_j_30_signed(*_Dphi_j_j_30_signed);

  } // end event loop
  cout << '\n';

  cout << "Total nominal weight: " << hist_bin::total_weights[0] << '\n';
  cout << "Fiducial events: " << n_fiducial << '\n';

  cout << "\nBad weights\n          "
          "weights   events  nominal";
  cout << "\ninf,nan "
       << setw(9) << n_inf.weights
       << setw(9) << n_inf.events
       << setw(9) << n_inf.nominal;
  cout << "\nlarge   "
       << setw(9) << n_large.weights
       << setw(9) << n_large.events
       << setw(9) << n_large.nominal;
  cout <<'\n'<< endl;

  // Output =========================================================

  std::vector<std::vector<double>> all_bins;
  std::vector<std::string> all_bin_labels;

  TFile f(argv[1],"recreate");

  f.mkdir("vars")->cd(); // write variable-specific histograms in a directory

  for (const auto& h : hist::all) { // loop over histograms
    cout << h.name << endl;
    gDirectory->mkdir(h.name.c_str())->cd();
    const auto& axis = h->axis();

    for (unsigned bi=0; bi<axis.nbins(); ++bi) {
      auto& bin = h->bins()[bi];
      all_bins.emplace_back();
      for (unsigned wi=0; wi<hist_bin::total_weights.size(); ++wi) {
        auto& w = bin[wi];
        // scale to cross section
        w *= (xs_br_fe.GetVal() / hist_bin::total_weights[wi]);
        // join histograms
        all_bins.back().push_back(w);
      }
      // make bin labels
      all_bin_labels.emplace_back(h.name+bin_label(axis,bi+1));
    }

    // Save nominal histograms
    th1(h,0,"nominal",(h.name+" nominal distribution"));

    unsigned i1 = 1, i2 = 1; // weight indices
    for (const auto& w : _weights) {
      i2 += w.GetSize();
      const char *name = w.GetBranchName()+1;

      if ( strlen(name)>=4 && !memcmp(name+1,"qcd",3) ) { // QCD weights

        // QCD envelope histograms
        th1(h,[=](const hist_bin& b){ return b.min(i1,i2); },
            cat("down",name), cat(h.name,' ',name+1," down variation"));
        th1(h,[=](const hist_bin& b){ return b.max(i1,i2); },
            cat("up",name), cat(h.name,' ',name+1," up variation"));

        // Symmetrized QCD uncertainties
        th1(h,[=](const hist_bin& b){
              return std::max( b[0]-b.min(i1,i2), b.max(i1,i2)-b[0] ); },
            cat("err",name),
            cat(h.name,' ',name+1," symmetrized uncertainties"));

      } else { // PDF weights

        // construct correlation matrix for this histogram
        const auto m_cor = cor(cov( *h, i1, i2 ));

        th1(m_cor.err, cat("err",name), cat(h.name,' ',name+1," errors"),
          [&axis](unsigned i){ return bin_label(axis,i); }
        );

        // Save correlation matrix as a histogram
        mat_hist(m_cor.cor, cat("cor",name),
          cat(h.name,' ',name+1," correlation matrix"),
          [&axis](unsigned i){ return bin_label(axis,i); }, {-1,1}
        );

      }
      i1 = i2;
    }

    gDirectory->cd("..");
  } // end loop over histograms

  f.cd(); // cd back to the file
  f.mkdir("error_sources")->cd();

  std::vector<sym_mat<double>> big_cov;
  big_cov.reserve(4);

  for (const auto& w : _weights) {
    static unsigned i1 = 1, i2 = 1; // weight indices
    i2 += w.GetSize();
    const char *name = w.GetBranchName()+1;

    if ( strlen(name)>=4 && !memcmp(name+1,"qcd",3) ) { // QCD weights

      // Symmetrized QCD uncertainties
      std::vector<double> symm_err(all_bins.size(),0.);
      for (unsigned i=0; i<all_bins.size(); ++i) {
        const auto& b = all_bins[i];
        const auto begin = b.begin();
        symm_err[i] = std::max(
          b[0] - *std::min_element(begin+i1,begin+i2),
          *std::max_element(begin+i1,begin+i2) - b[0] );
      }

      th1(symm_err,
        cat("err",name), cat(name+1," errors"),
        [&](unsigned i){ return all_bin_labels.at(i-1); }
      );

      big_cov.emplace_back(symm_err);

    } else { // PDF weights

      // construct the big correlation matrix
      auto cov_all = cov( all_bins, i1, i2 );
      const auto cor_all = cor(cov_all);

      big_cov.emplace_back(std::move(cov_all));

      th1(cor_all.err,
        cat("err",name), cat(name+1," errors"),
        [&](unsigned i){ return all_bin_labels.at(i-1); }
      );

      mat_hist(cor_all.cor,
        cat("cor",name), cat(name+1," correlation matrix"),
        [&](unsigned i){ return all_bin_labels.at(i-1); }, {-1,1}
      );

    }
    i1 = i2;
  }

  f.cd(); // cd back to the file

  mat_hist(cor(big_cov[0]+big_cov[2]).cor,
    "cor_pdf4lhc_qcd",
    cat(_weights[0].GetBranchName()+2,'+',_weights[2].GetBranchName()+2,
        " correlation matrix"),
    [&](unsigned i){ return all_bin_labels.at(i-1); }, {-1,1}
  );

  mat_hist(cor(big_cov[1]+big_cov[3]).cor,
    "cor_nnpdf30_qcd_nnlops",
    cat(_weights[1].GetBranchName()+2,'+',_weights[3].GetBranchName()+2,
        " correlation matrix"),
    [&](unsigned i){ return all_bin_labels.at(i-1); }, {-1,1}
  );

  f.Write();
  cout <<'\n'<< f.GetName() << " \033[32m✔\033[0m" << endl;

  return 0;
}
