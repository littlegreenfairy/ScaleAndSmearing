// ResolutionAnalysis.C - A macro to analyze and visualize electron energy resolution
// from true resolution histograms
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooCrystalBall.h>
#include <RooFitResult.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>
#include <iostream>


void ResolutionAnalysis() {
    // Open the file containing the resolution histograms
    TFile* inputFile = TFile::Open("outputHistograms_MC.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file!" << std::endl;
        return;
    }
    
    // Retrieve the resolution histograms
    TH2D* h_truerawResolution_vs_pt = (TH2D*)inputFile->Get("h_truerawResolution_vs_pt");
    TH2D* h_truecorrectedSCResolution_vs_pt = (TH2D*)inputFile->Get("h_truecorrectedSCResolution_vs_pt");
    TH2D* h_trueecaltrkregResolution_vs_pt = (TH2D*)inputFile->Get("h_trueecaltrkregResolution_vs_pt");
    
    if (!h_truerawResolution_vs_pt || !h_truecorrectedSCResolution_vs_pt || !h_trueecaltrkregResolution_vs_pt) {
        std::cerr << "Error: Could not find resolution histograms in the file!" << std::endl;
        inputFile->Close();
        return;
    }
    
    // Set nicer titles for the histograms
    h_truerawResolution_vs_pt->SetTitle("Raw Supercluster Energy Resolution");
    h_truecorrectedSCResolution_vs_pt->SetTitle("Corrected Supercluster Energy Resolution");
    h_trueecaltrkregResolution_vs_pt->SetTitle("ECAL-Track Regression Energy Resolution");
    
    // Create resolution distribution histograms for specific pT bins
    double ptBins[] = {4, 7, 9, 11, 14, 20, 40};
    int nPtBins = sizeof(ptBins)/sizeof(double) - 1;
    
    // Create canvas with pads for resolution distributions
    TCanvas* c_resDistributions = new TCanvas("c_resDistributions", "Resolution Distributions", 2000, 900);
    c_resDistributions->Divide(3, 2);
    
    // Draw resolution distributions for each pT bin
    for (int i = 0; i < nPtBins && i < 6; i++) { // Show up to 6 bins
        c_resDistributions->cd(i+1);
        
        // Get the bin numbers for the lower and upper pT edges
        int binLow = h_truerawResolution_vs_pt->GetXaxis()->FindBin(ptBins[i] + 0.1);
        int binHigh = h_truerawResolution_vs_pt->GetXaxis()->FindBin(ptBins[i+1] - 0.1);
        
        // Project the 2D histograms onto the Y axis (resolution) for the current pT bin
        TH1D* h_rawRes = h_truerawResolution_vs_pt->ProjectionY(
            Form("h_rawRes_bin%d", i), binLow, binHigh);
        TH1D* h_correctedRes = h_truecorrectedSCResolution_vs_pt->ProjectionY(
            Form("h_correctedRes_bin%d", i), binLow, binHigh);
        TH1D* h_ecaltrkregRes = h_trueecaltrkregResolution_vs_pt->ProjectionY(
            Form("h_ecaltrkregRes_bin%d", i), binLow, binHigh);
        
        // Set titles and styles
        h_rawRes->SetTitle(Form("Resolution for %.1f < p_{T} < %.1f GeV", ptBins[i], ptBins[i+1]));
        h_rawRes->GetXaxis()->SetTitle("(E_{gen} - E_{reco})/E_{gen}");
        h_rawRes->GetYaxis()->SetTitle("Entries");
        h_rawRes->SetLineColor(kRed);
        h_rawRes->SetLineWidth(2);
        h_rawRes->SetStats(0);
        
        h_correctedRes->SetLineColor(kBlue);
        h_correctedRes->SetLineWidth(2);
        
        h_ecaltrkregRes->SetLineColor(kGreen+2);
        h_ecaltrkregRes->SetLineWidth(2);
        
        // Find the maximum value to properly scale the Y-axis
        double maxValue = std::max({h_rawRes->GetMaximum(), 
                                  h_correctedRes->GetMaximum()});
        h_rawRes->SetMaximum(maxValue * 1.1);
        
        // Draw the histograms
        h_rawRes->Draw("HIST");
        h_correctedRes->Draw("HIST SAME");
        //h_ecaltrkregRes->Draw("HIST SAME");
        
        // Add legend to each pad
        TLegend* padLeg = new TLegend(0.65, 0.7, 0.89, 0.89);
        padLeg->SetBorderSize(0);
        padLeg->SetFillStyle(0);
        padLeg->AddEntry(h_rawRes, "Raw SC", "l");
        padLeg->AddEntry(h_correctedRes, "Corrected SC", "l");
        //padLeg->AddEntry(h_ecaltrkregRes, "ECAL-Track Reg.", "l");
        padLeg->Draw();
    }
    
    // Create output file to save results
    TFile* outputFile = new TFile("ResolutionAnalysis.root", "RECREATE");
    
    // Save canvas
    c_resDistributions->Write();
    
    // Save canvas as image
    c_resDistributions->SaveAs("PlotConID2022/TrueResolutionfromEGMmcprompt/ResolutionDistributions.png");
    
    // Clean up
    outputFile->Close();
    inputFile->Close();
    delete c_resDistributions;
    
    std::cout << "Resolution analysis completed. Results saved to ResolutionAnalysis.root" << std::endl;
    std::cout << "Distribution plot saved as PNG image." << std::endl;
}