// Macro to compare projections from h_invMass_ECAL_corrected_scandsm and h_invMass_ECAL_corrected
// and plot the ratio of their invariant mass projections for each pt and run number bin
#include <TFile.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <iostream>

#define NbinsPt 6
#define NbinsRun 9

void PlotRatioScandsmVsCorrected() {


    // Open the ROOT file for MC and extract h_scale_vs_pt
    TFile *fileMC = TFile::Open("outputHistograms_MC.root");
    if (!fileMC || fileMC->IsZombie()) {
        std::cerr << "Error: Could not open outputHistograms_MC.root" << std::endl;
    } 

    TH2D *h_invmass_vs_pt_MC = (TH2D*)fileMC->Get("h_scale_vs_pt");


    // Open the ROOT file for data
    TFile *filedata = TFile::Open("outputHistograms_DATA_partF.root");
    if (!filedata || filedata->IsZombie()) {
        std::cerr << "Error: Could not open outputHistograms_DATA_partF.root" << std::endl;
        return;
    }

    // Get the histograms
    TH3D *h_corr = (TH3D*)filedata->Get("h_invMass_ECAL_corrected");
    TH3D *h_scandsm = (TH3D*)filedata->Get("h_invMass_ECAL_corrected_scandsm");
    if (!h_corr || !h_scandsm) {
        std::cerr << "Error: Could not find one or both TH3D histograms in filedata." << std::endl;
        return;
    }

    // Bin definitions (update if needed)
    double runBins[] = {356309, 356900, 357538, 357732, 360000, 360400, 361000, 361600, 362200, 362760}; // 9 bins
    double bincenters[NbinsPt] = {5.5, 8, 10, 12.5, 17, 30};
    double binhalfwidths[NbinsPt] = {1.5, 1, 1, 1.5, 3, 10};

    // Loop over pt bins
    for (int i = 0; i < NbinsPt; ++i) {
        int binLow = i+1;
        int binHigh = i+1;

        // Project to 2D for this pt bin
        h_corr->GetXaxis()->SetRange(binLow, binHigh);
        TH2D *h2D_corr = (TH2D*)h_corr->Project3D("yz");
        h2D_corr->SetName(Form("h2D_corrscale_ptbin%d", i+1));

        h_scandsm->GetXaxis()->SetRange(binLow, binHigh);
        TH2D *h2D_scandsm = (TH2D*)h_scandsm->Project3D("yz");
        h2D_scandsm->SetName(Form("h2D_scandsm_ptbin%d", i+1));

        // Extract the ProjectionY of h_invmass_vs_pt_MC for bin i+1 (pt bin)
        TH1D *h1_invmass_vs_pt_MC = nullptr;
        h1_invmass_vs_pt_MC = h_invmass_vs_pt_MC->ProjectionY(Form("h1_invmass_vs_pt_MC_ptbin%d", i+1), i+1, i+1);

        // Loop over run number bins
        for (int j = 0; j < NbinsRun; ++j) {
            // Project to 1D invariant mass for this run bin
            TH1D *h1_corr = (TH1D*)h2D_corr->ProjectionX(Form("h1_corrscale_pt%d_run%d", i+1, j+1), j+1, j+1);
            TH1D *h1_scandsm = (TH1D*)h2D_scandsm->ProjectionX(Form("h1_scandsm_pt%d_run%d", i+1, j+1), j+1, j+1);

            // Check for valid histograms
            if (!h1_corr || !h1_scandsm) continue;
            if (h1_corr->GetEntries() == 0 || h1_scandsm->GetEntries() == 0) continue;

            // Apply mass window cut: only keep bins between 2.6 and 3.4 GeV
            int bin_min = h1_corr->GetXaxis()->FindBin(2.6);
            int bin_max = h1_corr->GetXaxis()->FindBin(3.4);
            for (int b = 1; b <= h1_corr->GetNbinsX(); ++b) {
                if (b < bin_min || b > bin_max) {
                    h1_corr->SetBinContent(b, 0);
                    h1_corr->SetBinError(b, 0);
                    h1_scandsm->SetBinContent(b, 0);
                    h1_scandsm->SetBinError(b, 0);
                    h1_invmass_vs_pt_MC->SetBinContent(b, 0);
                    h1_invmass_vs_pt_MC->SetBinError(b, 0);
                }
            }

            // Normalize histograms to unit area (integral)
            double int_corr = h1_corr->Integral();
            double int_scandsm = h1_scandsm->Integral();
            double int_mc = h1_invmass_vs_pt_MC->Integral();
            if (int_corr > 0) h1_corr->Scale(1.0 / int_corr);
            if (int_scandsm > 0) h1_scandsm->Scale(1.0 / int_scandsm);
            if (int_mc > 0) h1_invmass_vs_pt_MC->Scale(1.0 / int_mc);

            // Ratio before smearing corrections (corr / MC)
            TH1D *h_ratio = (TH1D*)h1_corr->Clone(Form("h_ratio_pt%d_run%d", i+1, j+1));
            h_ratio->SetTitle(Form("Ratio: corrected/MC, Pt bin %d, Run bin %d", i+1, j+1));
            h_ratio->Divide(h1_invmass_vs_pt_MC);
            h_ratio->GetYaxis()->SetTitle("Ratio corrected/MC");
            h_ratio->GetXaxis()->SetTitle("Invariant Mass [GeV]");
            h_ratio->SetMarkerStyle(20);
            h_ratio->SetMarkerColor(kBlue);
            h_ratio->SetLineColor(kBlue);
            h_ratio->SetLineWidth(2);

            // Ratio after smearing corrections (scandsm / MC)
            TH1D *h_ratio2 = (TH1D*)h1_scandsm->Clone(Form("h_ratio2_pt%d_run%d", i+1, j+1));
            h_ratio2->SetTitle(Form("Ratio: scandsm/MC, Pt bin %d, Run bin %d", i+1, j+1));
            h_ratio2->Divide(h1_invmass_vs_pt_MC);
            h_ratio2->GetYaxis()->SetTitle("Ratio scandsm/MC");
            h_ratio2->GetXaxis()->SetTitle("Invariant Mass [GeV]");
            h_ratio2->SetMarkerStyle(21);
            h_ratio2->SetMarkerColor(kRed);
            h_ratio2->SetLineColor(kRed);
            h_ratio2->SetLineWidth(2);

            // Draw both ratios on the same canvas with points and error bars, y axis centered around 1
            TCanvas *c = new TCanvas(Form("c_ratio_pt%d_run%d", i+1, j+1), Form("Ratio Pt%d Run%d", i+1, j+1), 800, 600);
            gStyle->SetOptStat(0);
            h_ratio->SetMinimum(0.7);

            // Set axis ranges
            h_ratio->GetXaxis()->SetRangeUser(2.6, 3.4);
            h_ratio2->GetXaxis()->SetRangeUser(2.6, 3.4);
            h_ratio->SetMinimum(0.0);
            h_ratio->SetMaximum(2.0);
            h_ratio2->SetMinimum(0.0);
            h_ratio2->SetMaximum(2.0);

            h_ratio->Draw("P E");
            h_ratio2->Draw("P E SAME");

            // Add legend
            TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
            leg->AddEntry(h_ratio, "Corrected/Data / MC", "p");
            leg->AddEntry(h_ratio2, "Scandsm/Data / MC", "p");
            leg->Draw();

            // Draw a horizontal line at y=1
            TLine *line = new TLine(2.6, 1.0, 3.4, 1.0);
            line->SetLineStyle(2);
            line->SetLineColor(kBlack);
            line->Draw();

            c->SaveAs(Form("PlotConID2022/RatioPostSmearingCorr/Pt_bin%d/Ratio_Pt%d_Run%d.png", i+1, i+1, j+1));
            delete leg;
            delete line;
            delete c;
            delete h_ratio;
            delete h_ratio2;
            delete h1_corr;
            delete h1_scandsm;
        }
        delete h2D_corr;
        delete h2D_scandsm;
    }
    filedata->Close();
}
