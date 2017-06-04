#include <iostream>
#include <vector>
#include <stdexcept>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>

#include "catstr.hh"
#include "math.hh"

using std::cout;
using std::cerr;
using std::endl;
using ivanp::cat;
using namespace ivanp::math;

template <typename T>
auto* get(TDirectory* dir, const char* namecycle) {
  auto* obj = dir->Get(namecycle);
  if (!obj) throw std::runtime_error(cat(
    "no object ",namecycle, " in directory ",dir->GetName()));
  return dynamic_cast<T*>(obj);
}

void style(TH1* h, Color_t c) {
  h->SetLineWidth(2);
  h->SetLineColor(c);
  h->SetMarkerSize(0);
  h->SetMarkerColor(c);
  h->GetXaxis()->SetLabelSize(0.025);
}

int main(int argc, char* argv[]) {
  if (argc!=3) {
    cout << "usage: " << argv[0] << " in.root out.pdf" << endl;
    return 1;
  }

  TFile f(argv[1]);
  auto* dir = get<TDirectory>(&f,"error_sources");

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2f");

  TCanvas canv;
  canv.SaveAs(cat(argv[2],'[').c_str());

  // ================================================================

  TH1 *xsec = get<TH1>(dir,"nominal_xsec");
  xsec->Scale(1e3);
  xsec->SetYTitle("#sigma  [fb]");
  style(xsec,602);
  xsec->Draw();
  canv.SaveAs(argv[2]);

  // ----------------------------------------------------------------

  TLegend leg(0.7,0.7,0.88,0.9);
  leg.SetBorderSize(0);
  leg.SetFillColorAlpha(0,0);

  TH1 *qcd_err = get<TH1>(dir,"err_ggf_qcd_2017");
  qcd_err->Scale(1e3);
  qcd_err->SetTitle("Absolute uncertainties;;#Delta#sigma  [fb]");
  qcd_err->GetYaxis()->SetRangeUser(0,2.5);
  style(qcd_err,602);
  qcd_err->Draw();
  leg.AddEntry(qcd_err,"QCD");

  TH1 *pdf_err = get<TH1>(dir,"err_pdf4lhc_unc");
  pdf_err->Scale(1e3);
  style(pdf_err,46);
  pdf_err->Draw("same");
  leg.AddEntry(pdf_err,"PDF");

  TH1 *sum_err = static_cast<TH1*>(qcd_err->Clone("err_sum"));
  for (int i=1, n=sum_err->GetNbinsX(); i<=n; ++i)
    sum_err->SetBinContent(i,
      std::sqrt(sq( qcd_err->GetBinContent(i), pdf_err->GetBinContent(i) )) );
  style(sum_err,8);
  sum_err->Draw("same");
  leg.AddEntry(sum_err,"Total");

  leg.Draw();
  canv.SaveAs(argv[2]);

  // ----------------------------------------------------------------

  qcd_err->SetTitle("Fractional uncertainties;;#Delta#sigma / #sigma");
  qcd_err->GetYaxis()->SetRangeUser(0,0.5);
  qcd_err->Divide(xsec);
  qcd_err->Draw();

  pdf_err->Divide(xsec);
  pdf_err->Draw("same");

  sum_err->Divide(xsec);
  sum_err->Draw("same");

  leg.Draw();
  canv.SaveAs(argv[2]);

  // ================================================================

  TH2 *qcd_cor = get<TH2>(dir,"cor_ggf_qcd_2017");
  qcd_cor->GetXaxis()->SetLabelSize(0.025);
  qcd_cor->GetYaxis()->SetLabelSize(0.025);
  qcd_cor->Draw();
  canv.SaveAs(argv[2]);

  TH2 *pdf_cor = get<TH2>(dir,"cor_pdf4lhc_unc");
  pdf_cor->GetXaxis()->SetLabelSize(0.025);
  pdf_cor->GetYaxis()->SetLabelSize(0.025);
  pdf_cor->Draw();
  canv.SaveAs(argv[2]);

  // ================================================================

  canv.SaveAs(cat(argv[2],']').c_str());
}
