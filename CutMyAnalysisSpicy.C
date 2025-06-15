#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooCrystalBall.h>
#include <RooFitResult.h>
#include <iostream>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TF1.h>
#include <RooPolynomial.h>
#include <TColor.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....Funzione globale per scritta Dati....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void WriteData(){
// Aggiungi una scritta in testo italic al bordo esterno della cornice
    TLatex latex;
    latex.SetTextSize(0.03); // Dimensione del testo
    latex.SetTextFont(42);   // Tipo di font
    latex.SetTextAlign(11);  // Allineamento: 11 = sinistra, alto
    latex.SetTextAngle(0);   // Angolo del testo
    latex.DrawLatexNDC(0.17, 0.92, "#it{2022 Double electron Parking - 17.7819 fb^{-1}}"); // Posizione e testo

}

void SetPlotMargins(TPad* pad, Double_t leftMargin = 0.15, Double_t rightMargin = 0.05) {
    pad->SetLeftMargin(leftMargin);
    pad->SetRightMargin(rightMargin);
}


void CutMyAnalysisSpicy()
{
    ///////////////Definisco i colori
    int bluCMS = TColor::GetColor("#5790fc");
    int rosaCMS = TColor::GetColor("#964a8b");
    int rossoCMS = TColor::GetColor("#e42536");
    int gialloCMS = TColor::GetColor("#ffa90e");
    int violaCMS = TColor::GetColor("#7a21dd");
    int mattoneCMS = TColor::GetColor("#bd1f01");
    ///////////////
    // Apri i file ROOT
    TFile *dati_cuts = TFile::Open("cutcomparison_DATA_partF.root");
    TFile *dati = TFile::Open("outputHistograms_DATA_partF.root");
    //TFile *MC_cuts = TFile::Open("cutcomparison_MC_LLcut.root");
    //TFile *MC = TFile::Open("outputHistograms_MC_LLcut.root");

    // Estrai gli istogrammi dal file dati_cuts
    TH1F *h_invMass_ECAL_ele = (TH1F *)dati_cuts->Get("h_invMass_withcuts");
    TH1F *h_invMass_ECAL_ele_nocuts = (TH1F *)dati_cuts->Get("h_invMass_nocuts");

    /*// Massa invariante dal MC
    TH1F *h_invMass_ECAL_ele_MC = (TH1F *)MC_cuts->Get("h_invMass_ECAL_ele");
    TH1F *h_invMass_ECAL_ele_MC_nocuts = (TH1F *)MC_cuts->Get("h_invMass_ECAL_ele_nocuts");

    // Estrai gli istogrammi dell'ID dal file dati
    TH1F *h_LowE_pfmvaIdEle = (TH1F *)dati->Get("h_LowE_pfmvaIdEle");
    TH1F *h_LowE_pfmvaIdEle_sign = (TH1F *)dati->Get("h_LowE_pfmvaIdEle_sign");
    TH1F *h_LowE_pfmvaIdEle_back = (TH1F *)dati->Get("h_LowE_pfmvaIdEle_back");

    // Estrai gli istogrammi dell'ID dal file del MC
    TH1F *h_LowE_pfmvaIdEle_MC = (TH1F *)MC->Get("h_LowE_pfmvaIdEle");
    TH1F *h_LowE_pfmvaIdEle_sign_MC = (TH1F *)MC->Get("h_LowE_pfmvaIdEle_sign");
    TH1F *h_LowE_pfmvaIdEle_back_MC = (TH1F *)MC->Get("h_LowE_pfmvaIdEle_back");


    //// Istogrammi massa invariante ECAL, Raw e con info tracciatore SENZA TAGLI (situazione iniziale)
    // Dati
    TH1F *h_invMass_data_nocuts = (TH1F *)dati_cuts->Get("h_invMass_nocuts");
    TH1F *h_invMass_RawSC_data_nocuts = (TH1F *)dati_cuts->Get("h_invMass_rawSC_nocuts");
    // MC
    TH1F *h_invMass_MC_nocuts = (TH1F *)MC_cuts->Get("h_invMass_nocuts");
    TH1F *h_invMass_RawSC_MC_nocuts = (TH1F *)MC_cuts->Get("h_invMass_rawSC_nocuts");

    //Distribuzione eta
    TH1F *h_eta_mc = (TH1F *)MC->Get("h_etaEle");
    TH1F *h_eta_dati = (TH1F *)dati->Get("h_etaEle");*/

    // Prima canvas: Sovrapposizione h_invMass_ECAL_ele e h_invMass_ECAL_ele_nocuts
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass Comparison", 800, 600);

    h_invMass_ECAL_ele->SetXTitle("m(e^{+}e^{-}) [GeV]");
    h_invMass_ECAL_ele->SetTitle("");
    h_invMass_ECAL_ele->SetLineColor(kBlack);
    h_invMass_ECAL_ele->SetFillColorAlpha(rossoCMS, 1);
    h_invMass_ECAL_ele->SetLineWidth(2);
    h_invMass_ECAL_ele->SetStats(kFALSE);

    h_invMass_ECAL_ele_nocuts->SetXTitle("m(e^{+}e^{-}) [GeV]");
    h_invMass_ECAL_ele_nocuts->SetTitle("");
    //h_invMass_ECAL_ele_nocuts->Scale(1.0 / h_invMass_ECAL_ele_nocuts->Integral());
    h_invMass_ECAL_ele_nocuts->SetLineColor(kBlack);
    h_invMass_ECAL_ele_nocuts->SetFillColorAlpha(gialloCMS, 1);
    h_invMass_ECAL_ele_nocuts->SetLineWidth(2);
    h_invMass_ECAL_ele_nocuts->SetStats(kFALSE);

    h_invMass_ECAL_ele_nocuts->Draw("HIST");
    h_invMass_ECAL_ele->Draw("HIST SAME");
    WriteData();

    // Aggiungi una legenda
    TLegend *leg1 = new TLegend(0.55, 0.77, 0.9, 0.9);
    // leg1->AddEntry(h_invMass_ECAL_ele_MC, "MC (#it{B #rightarrow K J/#psi})", "f");
    leg1->AddEntry(h_invMass_ECAL_ele, "With cuts on mva ID", "f");
    leg1->AddEntry(h_invMass_ECAL_ele_nocuts, "Without cuts", "f");
    leg1->SetTextSize(0.03);
    leg1->Draw();

    c1->SaveAs("PlotConID2022/InvariantMassComparison.png");

   /* // Seconda canvas: Sovrapposizione h_LowE_pfmvaIdEle, h_LowE_pfmvaIdEle_sign, h_LowE_pfmvaIdEle_back
    TCanvas *c2 = new TCanvas("c2", "Multivariate ID Comparison", 1000, 600);

    // Normalizza gli istogrammi
    h_LowE_pfmvaIdEle->Scale(1.0 / h_LowE_pfmvaIdEle->Integral());
    h_LowE_pfmvaIdEle_sign->Scale(1.0 / h_LowE_pfmvaIdEle_sign->Integral());
    h_LowE_pfmvaIdEle_back->Scale(1.0 / h_LowE_pfmvaIdEle_back->Integral());

    h_LowE_pfmvaIdEle_back->SetXTitle("lowE_MVAID");
    h_LowE_pfmvaIdEle_sign->SetXTitle("lowE_MVAID");
    h_LowE_pfmvaIdEle->SetXTitle("lowE_MVAID");
    h_LowE_pfmvaIdEle_back->SetTitle("Multivariate ID distribution - data");
    h_LowE_pfmvaIdEle_back->GetXaxis()->SetRangeUser(-12, 8);
    h_LowE_pfmvaIdEle_back->SetLineColor(gialloCMS);
    h_LowE_pfmvaIdEle_back->SetStats(kFALSE);
    // h_LowE_pfmvaIdEle_back->SetFillColorAlpha(gialloCMS, 0.5);
    h_LowE_pfmvaIdEle_back->SetLineWidth(5);

    h_LowE_pfmvaIdEle->SetLineColor(violaCMS);
    h_LowE_pfmvaIdEle->SetStats(kFALSE);
    // h_LowE_pfmvaIdEle->SetFillColorAlpha(bluCMS, 0.5);
    h_LowE_pfmvaIdEle->SetLineWidth(5);

    h_LowE_pfmvaIdEle_sign->SetLineColor(rossoCMS);
    h_LowE_pfmvaIdEle_sign->SetStats(kFALSE);
    // h_LowE_pfmvaIdEle_sign->SetFillColorAlpha(rossoCMS, 0.5);
    h_LowE_pfmvaIdEle_sign->SetLineWidth(5);

    h_LowE_pfmvaIdEle_sign->Draw("HIST");
    h_LowE_pfmvaIdEle->Draw("HIST SAME");
    h_LowE_pfmvaIdEle_back->Draw("HIST SAME");

    // ** Crea una copia dell'istogramma del segnale per colorare solo x > -1 **
    TH1F *h_LowE_pfmvaIdEle_sign_filled = (TH1F*)h_LowE_pfmvaIdEle_sign->Clone("h_LowE_pfmvaIdEle_sign_filled");

    // Cicla sui bin e mantieni solo quelli per cui x > -1
    for (int i = 1; i <= h_LowE_pfmvaIdEle_sign_filled->GetNbinsX(); ++i) {
    if (h_LowE_pfmvaIdEle_sign_filled->GetBinCenter(i) <= -1) {
        h_LowE_pfmvaIdEle_sign_filled->SetBinContent(i, 0);  // Imposta a 0 il contenuto per x <= -1
    }
    }

    // Imposta colore e riempi solo per x > -1

    h_LowE_pfmvaIdEle_sign_filled->SetFillColorAlpha(rossoCMS, 0.5);  // Colore di riempimento con trasparenza
    h_LowE_pfmvaIdEle_sign_filled->SetLineWidth(0);
    h_LowE_pfmvaIdEle_sign_filled->Draw("HIST SAME");

    // Aggiungi una legenda
    TLegend *leg2 = new TLegend(0.7, 0.77, 0.9, 0.9);
    leg2->AddEntry(h_LowE_pfmvaIdEle, "Full distribution", "l");
    leg2->AddEntry(h_LowE_pfmvaIdEle_sign, "Signal region", "l"); //(2.5 GeV < m < 3.2 GeV)
    leg2->AddEntry(h_LowE_pfmvaIdEle_back, "Sidebands", "l");
    leg2->SetTextSize(0.03); 
    leg2->Draw();

    c2->SaveAs("Plot Presentazione/MultivariateIDComparison_data.png");

    /// Plotto la distribuzione dell'ID nel MC
    TCanvas *c3 = new TCanvas("c3", "Multivariate ID - Monte Carlo", 800, 600);

    h_LowE_pfmvaIdEle_MC->Scale(1.0 / h_LowE_pfmvaIdEle_MC->Integral());
    h_LowE_pfmvaIdEle_MC->SetLineColor(bluCMS);
    h_LowE_pfmvaIdEle_MC->SetXTitle("lowE_MVAID");
    h_LowE_pfmvaIdEle_MC->SetTitle("Multivariate ID distribution - MC");
    h_LowE_pfmvaIdEle_MC->SetStats(kFALSE);
    // h_LowE_pfmvaIdEle->SetFillColorAlpha(bluCMS, 0.5);
    h_LowE_pfmvaIdEle_MC->SetLineWidth(5);
    h_LowE_pfmvaIdEle_MC->Draw("HIST");

    // ** Crea una copia dell'istogramma del segnale per colorare solo x > -1 **
    TH1F *h_LowE_pfmvaIdEle_MC_filled = (TH1F*) h_LowE_pfmvaIdEle_MC->Clone("h_LowE_pfmvaIdEle_MC_filled");

    // Cicla sui bin e mantieni solo quelli per cui x > -1
    for (int i = 1; i <= h_LowE_pfmvaIdEle_MC_filled->GetNbinsX(); ++i) {
    if (h_LowE_pfmvaIdEle_MC_filled->GetBinCenter(i) <= -1) {
        h_LowE_pfmvaIdEle_MC_filled->SetBinContent(i, 0);  // Imposta a 0 il contenuto per x <= -1
    }
    }

    // Imposta colore e riempi solo per x > -1

    h_LowE_pfmvaIdEle_MC_filled->SetFillColorAlpha(bluCMS, 0.5);  // Colore di riempimento con trasparenza
    h_LowE_pfmvaIdEle_MC_filled->SetLineWidth(0);
    h_LowE_pfmvaIdEle_MC_filled->Draw("HIST SAME");

    c3->SaveAs("Plot Presentazione/MultivariateIDComparison_MC.png");

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    // plot varie masse invarianti
    TCanvas *c4 = new TCanvas("c4", "Invariant masses", 2000, 600);
    c4->Divide(3, 1);
    // Normalizza gli istogrammi
    h_invMass_RawSC_data_nocuts->Scale(1.0 / h_invMass_RawSC_data_nocuts->Integral());
    h_invMass_data_nocuts->Scale(1.0 / h_invMass_data_nocuts->Integral());
    h_invMass_RawSC_MC_nocuts->Scale(1.0 / h_invMass_RawSC_MC_nocuts->Integral());
    h_invMass_MC_nocuts->Scale(1.0 / h_invMass_MC_nocuts->Integral());
    h_invMass_ECAL_ele_MC_nocuts->Scale(1.0 / h_invMass_ECAL_ele_MC_nocuts->Integral());
    //Massa raw
    c4->cd(1);
    SetPlotMargins((TPad*)c4->cd(1)); 
    h_invMass_RawSC_MC_nocuts->SetXTitle("m(e^{+}e^{-}) [GeV]");
    h_invMass_RawSC_MC_nocuts->SetTitle("Raw SC invariant mass");
    h_invMass_RawSC_MC_nocuts->SetLineColor(kBlack);
    h_invMass_RawSC_MC_nocuts->SetFillColorAlpha(bluCMS, 0.6);
    h_invMass_RawSC_MC_nocuts->SetLineWidth(2);
    h_invMass_RawSC_MC_nocuts->SetStats(kFALSE);
    h_invMass_RawSC_MC_nocuts->Draw("HIST");

    h_invMass_RawSC_data_nocuts->SetXTitle("m(e^{+}e^{-}) [GeV]");
    h_invMass_RawSC_data_nocuts->SetTitle("");
    h_invMass_RawSC_data_nocuts->SetLineColor(kBlack);
    h_invMass_RawSC_data_nocuts->SetFillColorAlpha(rossoCMS, 0.7);
    h_invMass_RawSC_data_nocuts->SetLineWidth(2);
    h_invMass_RawSC_data_nocuts->SetStats(kFALSE);
    h_invMass_RawSC_data_nocuts->Draw("HIST SAME");

    // Aggiungi una legenda
    TLegend *leg_raw = new TLegend(0.7, 0.77, 0.95, 0.9);
    leg_raw->AddEntry(h_invMass_RawSC_MC_nocuts, "MC", "f");
    leg_raw->AddEntry(h_invMass_RawSC_data_nocuts, "Data", "f");
    leg_raw->SetTextSize(0.03);
    leg_raw->Draw();

    //Massa ECAL
    c4->cd(2);
    SetPlotMargins((TPad*)c4->cd(2)); 
    h_invMass_ECAL_ele_MC_nocuts->SetXTitle("m(e^{+}e^{-}) [GeV]");
    h_invMass_ECAL_ele_MC_nocuts->SetTitle("ECAL invariant mass");
    h_invMass_ECAL_ele_MC_nocuts->SetLineColor(kBlack);
    h_invMass_ECAL_ele_MC_nocuts->SetFillColorAlpha(bluCMS, 0.6);
    h_invMass_ECAL_ele_MC_nocuts->SetLineWidth(2);
    h_invMass_ECAL_ele_MC_nocuts->SetStats(kFALSE);
    h_invMass_ECAL_ele_MC_nocuts->Draw("HIST");

    h_invMass_ECAL_ele_nocuts->SetFillColorAlpha(rossoCMS, 0.6);
    h_invMass_ECAL_ele_nocuts->Draw("HIST SAME");

    TLegend *leg_ECAL = new TLegend(0.7, 0.77, 0.95, 0.9);
    leg_ECAL->AddEntry(h_invMass_ECAL_ele_MC_nocuts, "MC", "f");
    leg_ECAL->AddEntry(h_invMass_ECAL_ele_nocuts, "Data", "f");
    leg_ECAL->SetTextSize(0.03);
    leg_ECAL->Draw();

    //Massa con track info
    c4->cd(3);
    SetPlotMargins((TPad*)c4->cd(3)); 
    h_invMass_MC_nocuts->SetXTitle("m(e^{+}e^{-}) [GeV]");
    h_invMass_MC_nocuts->SetTitle("Invariant mass with track info");
    h_invMass_MC_nocuts->SetLineColor(kBlack);
    h_invMass_MC_nocuts->SetFillColorAlpha(bluCMS, 0.6);
    h_invMass_MC_nocuts->SetLineWidth(2);
    h_invMass_MC_nocuts->SetStats(kFALSE);
    h_invMass_MC_nocuts->Draw("HIST");

    h_invMass_data_nocuts->SetXTitle("m(e^{+}e^{-}) [GeV]");
    h_invMass_data_nocuts->SetTitle("");
    h_invMass_data_nocuts->SetLineColor(kBlack);
    h_invMass_data_nocuts->SetFillColorAlpha(rossoCMS, 0.7);
    h_invMass_data_nocuts->SetLineWidth(2);
    h_invMass_data_nocuts->SetStats(kFALSE);
    h_invMass_data_nocuts->Draw("HIST SAME");

    // Aggiungi una legenda
    TLegend *leg_track = new TLegend(0.7, 0.77, 0.95, 0.9);
    leg_track->AddEntry(h_invMass_MC_nocuts, "MC", "f");
    leg_track->AddEntry(h_invMass_data_nocuts, "Data", "f");
    leg_track->SetTextSize(0.03);
    leg_track->Draw();

    c4->SaveAs("Plot Presentazione/invMasses_datavsMC.png");

/////////////////////////////////////////////////// Distribuzione eta
TCanvas *c5 = new TCanvas("c5", "eta distribution", 600, 600);
gPad->SetMargin(0.15, 0.1, 0.1, 0.1); 
    // Normalizza gli istogrammi
    h_eta_dati->Scale(1.0 / h_eta_dati->Integral());
    h_eta_mc->Scale(1.0 / h_eta_mc->Integral());
    h_eta_mc->SetXTitle("#eta");
    h_eta_mc->SetYTitle("A.U.");
    h_eta_mc->SetTitle("Raw SC invariant mass");
    h_eta_mc->SetLineColor(kBlack);
    h_eta_mc->SetFillColorAlpha(bluCMS, 0.4);
    h_eta_mc->SetLineWidth(2);
    h_eta_mc->SetStats(kFALSE);
    h_eta_mc->Draw("HIST");

    h_eta_dati->SetXTitle("#eta");
    h_eta_dati->SetTitle("");
    h_eta_dati->SetLineColor(kBlack);
    h_eta_dati->SetFillColorAlpha(rossoCMS, 0.4);
    h_eta_dati->SetLineWidth(2);
    h_eta_dati->SetStats(kFALSE);
    h_eta_dati->Draw("HIST SAME");

    // Aggiungi una legenda
    TLegend *leg_eta = new TLegend(0.7, 0.77, 0.9, 0.9);
    leg_eta->AddEntry(h_eta_mc, "MC", "f");
    leg_eta->AddEntry(h_eta_dati, "Data", "f");
    leg_eta->SetTextSize(0.03);
    leg_eta->Draw();
    c5->SaveAs("Plot Presentazione/eta.png");

*/

delete c1;
/*delete c2;
delete c3;
delete c4;
delete c5;*/



dati_cuts->Close();
//MC_cuts->Close();
dati->Close();
//MC->Close();
delete dati_cuts;
//delete MC_cuts;

}
