#include <iostream>
#include <vector>
#include <stdexcept>

#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TMatrixD.h>
#include <TAxis.h>

#include "catstr.hh"
#include "cstr.hh"

using std::cout;
using std::cerr;
using std::endl;
using ivanp::cat;

TH1 *h1 = nullptr;

void do_hist(TObject* obj, const char* name=nullptr) {
  if (!obj) throw std::runtime_error("null TObject pointer");

  TH1 *h = dynamic_cast<TH1*>(obj);
  cout << h->GetName() << endl;
  h->Write(name);

  if (!h1) h1 = h;
}

void do_mat(TObject* obj, const char* name=nullptr) {
  if (!obj) throw std::runtime_error("null TObject pointer");

  TH2 *h = dynamic_cast<TH2*>(obj);
  cout << h->GetName() << endl;
  const unsigned n = h->GetNbinsX();

  TMatrixD m(n,n);
  for (unsigned i=0; i<n; ++i)
    for (unsigned j=0; j<n; ++j)
      m[i][j] = h->GetBinContent(i,j);
  m.Write(name ? name : cat("m_",h->GetName()).c_str());
}

int main(int argc, char* argv[]) {
  if (argc!=3 && argc!=5) {
    cout << "usage: " << argv[0] << " in.root out.root [PDF QCD]" << endl;
    return 1;
  }

  TFile fin(argv[1],"read");
  if (fin.IsZombie()) return 1;
  TFile fout(argv[2],"recreate");
  if (fout.IsZombie()) return 1;

  cout << argv[1] << " -> " << argv[2] << endl;

  TDirectory *dir = fin.GetDirectory("error_sources");
  
  if (argc==3) {
    TIter next(dir->GetListOfKeys());
    for (TKey *key; (key=static_cast<TKey*>(next())); ) {
      if (starts_with(key->GetName(),"err_")) // error histogram
        do_hist(key->ReadObj());
      else if (starts_with(key->GetName(),"cor_")) // correlation matrix
        do_mat(key->ReadObj());
    }
  } else {
    do_hist( dir->Get(cat("err_",argv[3]).c_str()), "h_PDF_SUMCOV_SD" );
    do_mat ( dir->Get(cat("cor_",argv[3]).c_str()), "m_PDF_CORR" );

    do_hist( dir->Get(cat("err_",argv[4]).c_str()), "h_QCD_SUMCOV_SD" );
    do_mat ( dir->Get(cat("cor_",argv[4]).c_str()), "m_QCD_CORR" );
  }

  TAxis *ax = h1->GetXaxis();
  std::vector<std::pair<std::string,int>> n_var_bins;
  for (int i=1, n=ax->GetNbins(); i<=n; ++i) {
    const char* label = ax->GetBinLabel(i);
    std::string var(label,strchr(label,' ')-label);

    auto it = std::find_if(n_var_bins.rbegin(), n_var_bins.rend(),
      [&](const auto& p){ return p.first == var; }
    );
    if (it==n_var_bins.rend()) n_var_bins.emplace_back(std::move(var),1);
    else ++it->second;
  }

  TH1I *hn = new TH1I("varName_and_Nbins","",n_var_bins.size(),0,n_var_bins.size());
  for (const auto& b : n_var_bins) {
    static int i = 1;
    hn->SetBinContent(i,b.second);
    hn->GetXaxis()->SetBinLabel(i,b.first.c_str());
    ++i;
  }

  fout.Write();

  return 0;
}
