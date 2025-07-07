// PlotScaleandSmearingVsPtRunN.C
// Macro to read and plot scale and smearing histograms vs pT and run number

#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TColor.h>
#include <TPaletteAxis.h>

void PlotScaleandSmearingVsPtRunN() {
    // Set better color palette for 2D plots
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    
    // Open the file containing the histograms
    TFile* inputFile = TFile::Open("scale_corrections.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open corrections file!" << std::endl;
        return;
    }
    
    // Retrieve the histograms
    TH2D* h_scale = (TH2D*)inputFile->Get("h_scale");
    TH2D* h_smearing = (TH2D*)inputFile->Get("h_smearing");
    TH2D* h_corr_1ele = (TH2D*)inputFile->Get("h_corr_1ele");
    
    if (!h_scale || !h_smearing || !h_corr_1ele) {
        std::cerr << "Error: Could not find required histograms in the file!" << std::endl;
        inputFile->Close();
        return;
    }
    
    // Create canvas for scale plot
    TCanvas* c_scale = new TCanvas("c_scale", "Scale vs pT and Run Number", 900, 700);
    c_scale->SetRightMargin(0.15); // Make room for color palette
    c_scale->SetLeftMargin(0.13);
    c_scale->SetBottomMargin(0.12);
    
    // Draw scale histogram with COLZ option
    h_scale->SetTitle("");
    h_scale->GetXaxis()->SetTitle("p_{T} [GeV]");
    h_scale->GetYaxis()->SetTitle("Run Number");
    h_scale->GetZaxis()->SetTitle("1 - m^{(J/#psi)}_{data}/m^{(J/#psi)}_{MC}");
    h_scale->Draw("COLZ");
    
    // Add CMS Preliminary label
    TLatex cmsLabel;
    cmsLabel.SetNDC();
    cmsLabel.SetTextFont(62);
    cmsLabel.SetTextSize(0.04);
    cmsLabel.SetTextAlign(11);
    cmsLabel.DrawLatex(0.19, 0.91, "CMS");
    
    // Add Preliminary in italics
    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextFont(52);
    prelimLabel.SetTextSize(0.04);
    prelimLabel.SetTextAlign(11);
    prelimLabel.DrawLatex(0.26, 0.91, "Preliminary");
    
    // Draw luminosity/energy information
    TLatex lumiLabel;
    lumiLabel.SetNDC();
    lumiLabel.SetTextFont(42);
    lumiLabel.SetTextSize(0.04);
    lumiLabel.SetTextAlign(31);
    lumiLabel.DrawLatex(0.85, 0.91, "38.01 fb^{-1} (2022, 13.7 TeV)");

    // Save scale plot
    c_scale->SaveAs("ScaleVsPtRunN.png");
    
    // Create canvas for smearing plot
    TCanvas* c_smearing = new TCanvas("c_smearing", "Smearing vs pT and Run Number", 900, 700);
    c_smearing->SetRightMargin(0.15); // Make room for color palette
    c_smearing->SetLeftMargin(0.13);
    c_smearing->SetBottomMargin(0.12);
    
    // Draw smearing histogram with COLZ option
    h_smearing->SetTitle("");
    h_smearing->GetXaxis()->SetTitle("p_{T} [GeV]");
    h_smearing->GetYaxis()->SetTitle("Run Number");
    h_smearing->GetZaxis()->SetTitle("#sigma_{data}/#sigma_{MC}");
    h_smearing->Draw("COLZ");
    
    // Add CMS Preliminary label
    cmsLabel.DrawLatex(0.19, 0.91, "CMS");
    prelimLabel.DrawLatex(0.26, 0.91, "Preliminary");
    lumiLabel.DrawLatex(0.85, 0.91, "38.01 fb^{-1} (2022, 13.7 TeV)");

    // Save smearing plot
    c_smearing->SaveAs("SmearingVsPtRunN.png");
    
    // Create canvas for single electron correction plot
    TCanvas* c_corr = new TCanvas("c_corr", "Single Electron Correction vs pT and Run Number", 900, 700);
    c_corr->SetRightMargin(0.15); // Make room for color palette
    c_corr->SetLeftMargin(0.13);
    c_corr->SetBottomMargin(0.12);
    
    // Draw correction histogram with COLZ option
    h_corr_1ele->SetTitle("");
    h_corr_1ele->GetXaxis()->SetTitle("p_{T} [GeV]");
    h_corr_1ele->GetYaxis()->SetTitle("Run Number");
    h_corr_1ele->GetZaxis()->SetTitle("m^{(J/#psi)}_{MC} / m^{(J/#psi)}_{data}");
    h_corr_1ele->Draw("COLZ");
    
    // Add CMS Preliminary label
    cmsLabel.DrawLatex(0.19, 0.91, "CMS");
    prelimLabel.DrawLatex(0.26, 0.91, "Preliminary");
    lumiLabel.DrawLatex(0.85, 0.91, "38.01 fb^{-1} (2022, 13.7 TeV)");

    // Save correction plot
    c_corr->SaveAs("CorrVsPtRunN.png");
    
    // Clean up
    inputFile->Close();
    delete c_scale;
    delete c_smearing;
    delete c_corr;
    
    std::cout << "Scale, smearing, and correction plots created successfully!" << std::endl;
}