Compilation: `make`

Produce correlation matrices and error histograms:
`./cor output.root input1.root ...`

`cor` creates a root file with a lot of things in it and `TH2D`s for matrices.

`convert` reads that output and makes a root file with the following structure:

```
TH1D h_PDF_SUMCOV_SD
TMatrixT<double> m_PDF_CORR
TH1D h_QCD_SUMCOV_SD
TMatrixT<double> m_QCD_CORR
TH1I varName_and_Nbins
```

Usage example:
```
./cor nnlops_correlations.root user*.root
./convert nnlops_correlations.root theory_unc.root pdf4lhc_unc ggf_qcd_2017
```

