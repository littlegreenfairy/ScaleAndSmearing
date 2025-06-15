#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooCrystalBall.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <iostream>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TF1.h>
#include <TColor.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>

void plotCorrectionsVsRunNumber() {
    // Open the ROOT file
    TFile *f = TFile::Open("scale_corrections.root");

    if (!f || f->IsZombie()) {
        std::cerr << "Failed to open file!" << std::endl;
        return;
    }

    // Get the histograms from the file
    TH2D *h_corr_1ele = (TH2D*)f->Get("h_corr_1ele");
    TH1D *h_corr_1ele_inclusiveRun = (TH1D*)f->Get("h_corr_1ele_inclusiveRun");

    if (!h_corr_1ele || !h_corr_1ele_inclusiveRun) {
        std::cerr << "Failed to retrieve histograms!" << std::endl;
        return;
    }

    // Define the number of pt bins and pt bin edges
    const int nPtBins = 6;
    double Ptbins[7] = {4, 7, 9, 11, 14, 20, 40}; // Example bin edges for Pt

    // Loop over each pt bin to extract projections and create canvases
    for (int i = 0; i < nPtBins; ++i) {
        // Define the range for the Pt bin (using Ptbins)
        double ptLow = Ptbins[i];
        double ptHigh = Ptbins[i + 1];

        // Extract the projection for this pt bin
        TH1D *h_proj = h_corr_1ele->ProjectionY(Form("h_proj_pt%d", i+1), i+1, i+1);
        h_proj->SetTitle(Form("Corrections vs Run number for %.1f < Pt < %.1f GeV", ptLow, ptHigh));
        h_proj->GetXaxis()->SetTitle("Run Number");
        h_proj->GetYaxis()->SetTitle("Correction");
        
        // Set stats=false
        h_proj->SetStats(kFALSE);
        
        // Set marker style
        h_proj->SetMarkerStyle(21);
        h_proj->SetMarkerSize(0.8);

        // Create a new canvas for this pt bin
        TCanvas *c = new TCanvas(Form("c_pt%d", i), Form("Canvas for Pt bin %.1f < Pt < %.1f GeV", ptLow, ptHigh), 800, 600);
        h_proj->Draw("P");  // Use P option to draw with markers

        // Get the inclusive correction for this pt bin from h_corr_1ele_inclusiveRun
        int bin = h_corr_1ele_inclusiveRun->FindBin((ptLow + ptHigh) / 2.0);
        double inclusiveCorrection = h_corr_1ele_inclusiveRun->GetBinContent(bin);
        double inclusiveCorrectionError = h_corr_1ele_inclusiveRun->GetBinError(bin);

        // Create a line representing the inclusive correction for this pt bin
        TLine *line = new TLine(h_proj->GetXaxis()->GetXmin(), inclusiveCorrection, h_proj->GetXaxis()->GetXmax(), inclusiveCorrection);
        line->SetLineColor(kRed);
        line->SetLineStyle(kDashed);
        line->SetLineWidth(2);
        line->Draw("same");

        // Create a shaded band representing +/- 1 sigma
        TBox *box = new TBox(h_proj->GetXaxis()->GetXmin(), inclusiveCorrection - inclusiveCorrectionError,
                             h_proj->GetXaxis()->GetXmax(), inclusiveCorrection + inclusiveCorrectionError);
        box->SetFillColorAlpha(kRed, 0.3);
        box->Draw("same");
        
        // Draw histogram points again to make sure they're visible on top of the band
        h_proj->Draw("P SAME");

        // Add a legend
        TLegend *legend = new TLegend(0.65, 0.15, 0.89, 0.30);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(h_proj, "Run-by-run corrections", "p");
        legend->AddEntry(line, "Inclusive correction", "l");
        legend->AddEntry(box, "Inclusive correction \pm 1#sigma", "f");
        legend->Draw();

        // Save the canvas as a PNG (optional)
        c->SaveAs(Form("PlotConID2022/ScalevsRunN/corrections_vs_run_pt%d.png", i+1));

        // Cleanup
        delete c;
    }

    // Close the file
    f->Close();
}
