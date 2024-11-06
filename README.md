# ScaleAndSmearing
Scale and smearing corrections for ECAL on J/psi

gROOT->SetMacroPath("/Users/elenadesantis/Documents/CERN2024/ScaleAndSmearing");


TFile *filehisto = TFile::Open("outputHistograms_DATA_partF.root");
TH2D *h2d = (TH2D*)filehisto->Get("invMass_vsrun_check_yz");
int nbin = h2d->GetNbinsY();
std::cout << nbin << std::endl;