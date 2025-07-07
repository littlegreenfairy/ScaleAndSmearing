#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TColor.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooCrystalBall.h>
#include <RooFitResult.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TF1.h>

using namespace std;

#define Nbins 6

void plotMassProjByPtBins() {
    // Open the file containing histograms
    TFile* file = TFile::Open("outputHistograms_MC.root", "READ");
    if (!file || file->IsZombie()) {
        cout << "Error opening file!" << endl;
        return;
    }
    
    // Count how many pT bins we have by checking for projections
    int nbinsPt = 0;
    while (true) {
        TH1D* h = (TH1D*)file->Get(Form("proj_bin_%d", nbinsPt + 1));
        if (!h) break;
        nbinsPt++;
    }
    
    if (nbinsPt == 0) {
        cout << "No projections found in the file!" << endl;
        file->Close();
        return;
    }
    
    cout << "Found " << nbinsPt << " pT bins" << endl;
    
    // Define y-axis maximum values for each pT bin
    double yMaxValues[] = {600, 400, 150, 110, 70, 25};
    
    // Fit parameters for regressed SC
    double means_reg[Nbins] = {3.035, 3.02, 3.04, 3.0, 3.0743, 3.0916};
    double sigmas_reg[Nbins] = {0.16, 0.1602, 0.13, 0.14, 0.1152, 0.12};
    double alphaLs_reg[Nbins] = {0.92, 0.92, 0.92, 0.92, 0.92, 0.92};
    double nLs_reg[Nbins] = {3.2, 3.2, 3.2, 3.2, 3.2, 3.2};
    double alphaRs_reg[Nbins] = {1.67, 1.67, 1.67, 1.67, 1.67, 1.67};
    double nRs_reg[Nbins] = {4.59, 4.59, 4.59, 4.59, 4.59, 4.59};
    double meanInf_reg[Nbins] = {3.01, 2.9, 3.0, 2.85, 2.9, 2.9};
    double meanSup_reg[Nbins] = {3.2, 3.2, 3.2, 3.2, 3.2, 3.2};
    
    // Fit parameters for raw SC (adjust these values as needed)
    double means_raw[Nbins] = {2.4, 2.4, 2.5, 2.6, 2.7, 2.8};
    double sigmas_raw[Nbins] = {0.2, 0.15, 0.15, 0.12, 0.1, 0.2};
    double alphaLs_raw[Nbins] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    double nLs_raw[Nbins] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
    double alphaRs_raw[Nbins] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
    double nRs_raw[Nbins] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
    double meanInf_raw[Nbins] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    double meanSup_raw[Nbins] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
    
    // Fit ranges for regressed SC
    double massMin_reg[Nbins] = {2.2, 2.2, 2.3, 2.1, 2.65, 2.65};
    double massMax_reg[Nbins] = {3.9, 3.9, 3.8, 3.8, 3.5, 3.5};
    
    // Fit ranges for raw SC (typically wider due to poorer resolution)
    double massMin_raw[Nbins] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5}; // Adjusted for better fit
    double massMax_raw[Nbins] = {4.6, 4.6, 4.6, 4.6, 4.6, 4.6}; // Adjusted for better fit
    
    // pT bin centers and half-widths for plotting sigma vs pT
    double ptBinCenters[Nbins] = {5.5, 8.0, 10.5, 12.5, 17.0, 30.0}; // GeV
    double ptBinHalfWidths[Nbins] = {1.5, 1.0, 1.5, 1.5, 3.0, 10.0}; // GeV
    
    // Arrays to store fit results
    double sigmas_fitted_reg[Nbins], sigmaErrors_reg[Nbins];
    double sigmas_fitted_raw[Nbins], sigmaErrors_raw[Nbins];
    
    // Calculate canvas layout
    int nCols = ceil(sqrt(nbinsPt));
    int nRows = ceil(double(nbinsPt) / nCols);
    
    // Create a canvas
    TCanvas* c = new TCanvas("c_mass_proj", "Mass Projections by pT Bin", 2000, 1000);
    c->Divide(nCols, nRows);
    
    // Loop over pT bins
    for (int i = 1; i <= nbinsPt; i++) {
        // Get projections
        TH1D* proj = (TH1D*)file->Get(Form("proj_bin_%d", i));
        TH1D* proj_rawSC = (TH1D*)file->Get(Form("RawSC_bin_%d", i));
        
        if (!proj || !proj_rawSC) {
            cout << "Missing projection for bin " << i << endl;
            continue;
        }
        
        // Rebin histograms by a factor of 2
        //proj->Rebin(2);
        //proj_rawSC->Rebin(2);
        
        proj->SetStats(kFALSE);
        proj_rawSC->SetStats(kFALSE);
        
        // Set different colors and styles
        proj->SetLineColor(kBlue);
        proj->SetLineWidth(2);
        proj_rawSC->SetLineColor(kRed);
        proj_rawSC->SetLineWidth(2);
        
        // Move to appropriate pad
        c->cd(i);
        
        // Fit regressed SC histogram
        RooRealVar mass_reg("mass_reg", "m(e^{+}e^{-})", massMin_reg[i-1], massMax_reg[i-1]);
        RooDataHist data_reg("data_reg", "Dataset from histogram", mass_reg, proj);
        
        RooRealVar cb_mean_reg(Form("cb_mean_reg_%d", i), "Mean of CB", means_reg[i-1], meanInf_reg[i-1], meanSup_reg[i-1]);
        RooRealVar cb_sigma_reg(Form("cb_sigma_reg_%d", i), "Sigma of CB", sigmas_reg[i-1], 0.1 * sigmas_reg[i-1], 2 * sigmas_reg[i-1]);
        RooRealVar cb_alphaL_reg(Form("cb_alphaL_reg_%d", i), "AlphaL of CB", alphaLs_reg[i-1], 0, 5.0);
        RooRealVar cb_nL_reg(Form("cb_nL_reg_%d", i), "nL of CB", nLs_reg[i-1], 100, 200.0);
        RooRealVar cb_alphaR_reg(Form("cb_alphaR_reg_%d", i), "AlphaR of CB", alphaRs_reg[i-1], 0, 25.0);
        RooRealVar cb_nR_reg(Form("cb_nR_reg_%d", i), "nR of CB", nRs_reg[i-1], 0.1, 200.0);
        
        RooCrystalBall cb_reg(Form("cb_reg_%d", i), "Asymmetric Crystal Ball", mass_reg, cb_mean_reg, cb_sigma_reg, cb_alphaL_reg, cb_nL_reg, cb_alphaR_reg, cb_nR_reg);
        
        // Fit raw SC histogram
        RooRealVar mass_raw("mass_raw", "m(e^{+}e^{-})", massMin_raw[i-1], massMax_raw[i-1]);
        RooDataHist data_raw("data_raw", "Dataset from histogram", mass_raw, proj_rawSC);
        
        RooRealVar cb_mean_raw(Form("cb_mean_raw_%d", i), "Mean of CB", means_raw[i-1], meanInf_raw[i-1], meanSup_raw[i-1]);
        RooRealVar cb_sigma_raw(Form("cb_sigma_raw_%d", i), "Sigma of CB", sigmas_raw[i-1], 0.1 * sigmas_raw[i-1], 2 * sigmas_raw[i-1]);
        RooRealVar cb_alphaL_raw(Form("cb_alphaL_raw_%d", i), "AlphaL of CB", alphaLs_raw[i-1], 0.1, 5.0);
        RooRealVar cb_nL_raw(Form("cb_nL_raw_%d", i), "nL of CB", nLs_raw[i-1], 100, 200.0);
        RooRealVar cb_alphaR_raw(Form("cb_alphaR_raw_%d", i), "AlphaR of CB", alphaRs_raw[i-1], 0.1, 25.0);
        RooRealVar cb_nR_raw(Form("cb_nR_raw_%d", i), "nR of CB", nRs_raw[i-1], 0.1, 200.0);
        
        RooCrystalBall cb_raw(Form("cb_raw_%d", i), "Asymmetric Crystal Ball", mass_raw, cb_mean_raw, cb_sigma_raw, cb_alphaL_raw, cb_nL_raw, cb_alphaR_raw, cb_nR_raw);
        
        // Create the frame for the plot
        double frameMin = min(massMin_reg[i-1], massMin_raw[i-1]);
        double frameMax = max(massMax_reg[i-1], massMax_raw[i-1]);
        
        // Create TGraphErrors for data points with error bars
        int nBins_reg = proj->GetNbinsX();
        int nBins_raw = proj_rawSC->GetNbinsX();
        
        vector<double> x_reg, y_reg, ex_reg, ey_reg;
        vector<double> x_raw, y_raw, ex_raw, ey_raw;
        
        // Fill regressed data points
        for (int bin = 1; bin <= nBins_reg; bin++) {
            double x = proj->GetBinCenter(bin);
            if (x >= frameMin && x <= frameMax && proj->GetBinContent(bin) > 0) {
                x_reg.push_back(x);
                y_reg.push_back(proj->GetBinContent(bin));
                ex_reg.push_back(proj->GetBinWidth(bin) / 2.0);
                ey_reg.push_back(proj->GetBinError(bin));
            }
        }
        
        // Fill raw data points
        for (int bin = 1; bin <= nBins_raw; bin++) {
            double x = proj_rawSC->GetBinCenter(bin);
            if (x >= frameMin && x <= frameMax && proj_rawSC->GetBinContent(bin) > 0) {
                x_raw.push_back(x);
                y_raw.push_back(proj_rawSC->GetBinContent(bin));
                ex_raw.push_back(proj_rawSC->GetBinWidth(bin) / 2.0);
                ey_raw.push_back(proj_rawSC->GetBinError(bin));
            }
        }
        
        // Create TGraphErrors
        TGraphErrors* gr_reg = new TGraphErrors(x_reg.size(), &x_reg[0], &y_reg[0], &ex_reg[0], &ey_reg[0]);
        TGraphErrors* gr_raw = new TGraphErrors(x_raw.size(), &x_raw[0], &y_raw[0], &ex_raw[0], &ey_raw[0]);
        
        // Set styles
        gr_reg->SetMarkerStyle(20);
        gr_reg->SetMarkerColor(kBlue);
        gr_reg->SetLineColor(kBlue);
        
        gr_raw->SetMarkerStyle(22);
        gr_raw->SetMarkerColor(kRed);
        gr_raw->SetLineColor(kRed);
        
        // Set up the plot range
        int maxIndex = min(i-1, 5);
        double yMax = yMaxValues[maxIndex];
        
        // Create a histogram for axes
        TH1F* hFrame = new TH1F(Form("hFrame_%d", i), "", 100, frameMin, frameMax);
        hFrame->SetMaximum(yMax);
        hFrame->SetMinimum(0);
        hFrame->GetXaxis()->SetTitle("m(ee) [GeV]");
        hFrame->GetYaxis()->SetTitle("Entries");
        hFrame->SetStats(kFALSE);
        
        // Draw the frame and data points
        hFrame->Draw();
        gr_raw->Draw("P SAME");
        gr_reg->Draw("P SAME");
        
        // Perform fits and get results
        RooFitResult *fit_result_reg = cb_reg.fitTo(data_reg, RooFit::Save());
        RooFitResult *fit_result_raw = cb_raw.fitTo(data_raw, RooFit::Save());

        // Create TF1 functions using explicit parameter values instead of lambda
        // Get fitted parameter values
        double mean_reg_val = cb_mean_reg.getVal();
        double sigma_reg_val = cb_sigma_reg.getVal();
        double alphaL_reg_val = cb_alphaL_reg.getVal();
        double nL_reg_val = cb_nL_reg.getVal();
        double alphaR_reg_val = cb_alphaR_reg.getVal();
        double nR_reg_val = cb_nR_reg.getVal();
        
        double mean_raw_val = cb_mean_raw.getVal();
        double sigma_raw_val = cb_sigma_raw.getVal();
        double alphaL_raw_val = cb_alphaL_raw.getVal();
        double nL_raw_val = cb_nL_raw.getVal();
        double alphaR_raw_val = cb_alphaR_raw.getVal();
        double nR_raw_val = cb_nR_raw.getVal();
        
        // Create TF1 functions for the fit curves using fitted parameters
        TF1* f_reg = new TF1(Form("f_reg_%d", i), "gaus", massMin_reg[i-1], massMax_reg[i-1]);
        f_reg->SetParameters(proj->GetMaximum(), mean_reg_val, sigma_reg_val);
        f_reg->SetLineColor(kBlue+2);
        f_reg->SetLineWidth(3);
        f_reg->SetLineStyle(2);
        
        TF1* f_raw = new TF1(Form("f_raw_%d", i), "gaus", massMin_raw[i-1], massMax_raw[i-1]);
        f_raw->SetParameters(proj_rawSC->GetMaximum(), mean_raw_val, sigma_raw_val);
        f_raw->SetLineColor(kRed+2);
        f_raw->SetLineWidth(3);
        f_raw->SetLineStyle(2);
        
        // Draw the fit functions
        f_reg->Draw("SAME");
        f_raw->Draw("SAME");
        
        // Print fit results and store sigma values
        if (fit_result_reg && fit_result_reg->status() == 0) {
            cout << "Regressed fit converged for bin " << i << endl;
            cout << "Mean: " << cb_mean_reg.getVal() << " ± " << cb_mean_reg.getError() << endl;
            cout << "Sigma: " << cb_sigma_reg.getVal() << " ± " << cb_sigma_reg.getError() << endl;
            sigmas_fitted_reg[i-1] = cb_sigma_reg.getVal();
            sigmaErrors_reg[i-1] = cb_sigma_reg.getError();
        } else {
            sigmas_fitted_reg[i-1] = 0;
            sigmaErrors_reg[i-1] = 0;
        }
      
        if (fit_result_raw && fit_result_raw->status() == 0) {
            cout << "Raw fit converged for bin " << i << endl;
            cout << "Mean: " << cb_mean_raw.getVal() << " ± " << cb_mean_raw.getError() << endl;
            cout << "Sigma: " << cb_sigma_raw.getVal() << " ± " << cb_sigma_raw.getError() << endl;
            sigmas_fitted_raw[i-1] = cb_sigma_raw.getVal();
            sigmaErrors_raw[i-1] = cb_sigma_raw.getError();
        } else {
            sigmas_fitted_raw[i-1] = 0;
            sigmaErrors_raw[i-1] = 0;
        }

        
        // Add fit parameter box
        TPaveText* fitBox = new TPaveText(0.62, 0.15, 0.98, 0.89, "NDC");
        fitBox->SetFillColor(0);
        fitBox->SetBorderSize(1);
        fitBox->SetTextSize(0.025);
        fitBox->SetTextAlign(12); // Left aligned
        
        // Add regressed fit parameters
        if (fit_result_reg && fit_result_reg->status() == 0) {
            fitBox->AddText(Form("Regressed: (converged)"));
            fitBox->AddText(Form("#mu = %.4f #pm %.4f", cb_mean_reg.getVal(), cb_mean_reg.getError()));
            fitBox->AddText(Form("#sigma = %.4f #pm %.4f", cb_sigma_reg.getVal(), cb_sigma_reg.getError()));
            fitBox->AddText(Form("#alpha_{L} = %.3f #pm %.3f", cb_alphaL_reg.getVal(), cb_alphaL_reg.getError()));
            fitBox->AddText(Form("n_{L} = %.2f #pm %.2f", cb_nL_reg.getVal(), cb_nL_reg.getError()));
            fitBox->AddText(Form("#alpha_{R} = %.3f #pm %.3f", cb_alphaR_reg.getVal(), cb_alphaR_reg.getError()));
            fitBox->AddText(Form("n_{R} = %.2f #pm %.2f", cb_nR_reg.getVal(), cb_nR_reg.getError()));
        } else {
            fitBox->AddText(Form("Regressed: (failed)"));
            fitBox->AddText(Form("#mu = %.4f", cb_mean_reg.getVal()));
            fitBox->AddText(Form("#sigma = %.4f", cb_sigma_reg.getVal()));
            fitBox->AddText(Form("#alpha_{L} = %.3f", cb_alphaL_reg.getVal()));
            fitBox->AddText(Form("n_{L} = %.2f", cb_nL_reg.getVal()));
            fitBox->AddText(Form("#alpha_{R} = %.3f", cb_alphaR_reg.getVal()));
            fitBox->AddText(Form("n_{R} = %.2f", cb_nR_reg.getVal()));
        }
        
        fitBox->AddText(""); // Empty line for separation
        
        // Add raw fit parameters
        if (fit_result_raw && fit_result_raw->status() == 0) {
            fitBox->AddText(Form("Raw: (converged)"));
            fitBox->AddText(Form("#mu = %.4f #pm %.4f", cb_mean_raw.getVal(), cb_mean_raw.getError()));
            fitBox->AddText(Form("#sigma = %.4f #pm %.4f", cb_sigma_raw.getVal(), cb_sigma_raw.getError()));
            fitBox->AddText(Form("#alpha_{L} = %.3f #pm %.3f", cb_alphaL_raw.getVal(), cb_alphaL_raw.getError()));
            fitBox->AddText(Form("n_{L} = %.2f #pm %.2f", cb_nL_raw.getVal(), cb_nL_raw.getError()));
            fitBox->AddText(Form("#alpha_{R} = %.3f #pm %.3f", cb_alphaR_raw.getVal(), cb_alphaR_raw.getError()));
            fitBox->AddText(Form("n_{R} = %.2f #pm %.2f", cb_nR_raw.getVal(), cb_nR_raw.getError()));
        } else {
            fitBox->AddText(Form("Raw: (failed)"));
            fitBox->AddText(Form("#mu = %.4f", cb_mean_raw.getVal()));
            fitBox->AddText(Form("#sigma = %.4f", cb_sigma_raw.getVal()));
            fitBox->AddText(Form("#alpha_{L} = %.3f", cb_alphaL_raw.getVal()));
            fitBox->AddText(Form("n_{L} = %.2f", cb_nL_raw.getVal()));
            fitBox->AddText(Form("#alpha_{R} = %.3f", cb_alphaR_raw.getVal()));
            fitBox->AddText(Form("n_{R} = %.2f", cb_nR_raw.getVal()));
        }
        
        fitBox->Draw();
        
        // Add a legend
        TLegend* leg = new TLegend(0.15, 0.75, 0.45, 0.89);
        
        leg->AddEntry(gr_reg, "Regressed Mass", "p");
        leg->AddEntry(gr_raw, "Raw Mass", "p");
        leg->AddEntry(f_reg, "Regressed Fit", "l");
        leg->AddEntry(f_raw, "Raw Fit", "l");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->Draw();
        
        // Clean up RooPlot frames
    }
    
    // Save the canvas
    c->Update(); // Force canvas update before saving
    c->SaveAs("mass_projections_with_fits.png");
    
    // Create sigma vs pT plot
    TCanvas* c_sigma = new TCanvas("c_sigma", "Sigma vs pT", 800, 600);
    
    // Create TGraphErrors for sigma vs pT
    TGraphErrors* gr_sigma_reg = new TGraphErrors(nbinsPt, ptBinCenters, sigmas_fitted_reg, ptBinHalfWidths, sigmaErrors_reg);
    TGraphErrors* gr_sigma_raw = new TGraphErrors(nbinsPt, ptBinCenters, sigmas_fitted_raw, ptBinHalfWidths, sigmaErrors_raw);
    
    // Set styles
    gr_sigma_reg->SetMarkerStyle(20);
    gr_sigma_reg->SetMarkerColor(kBlue);
    gr_sigma_reg->SetLineColor(kBlue);
    gr_sigma_reg->SetMarkerSize(1.2);
    
    gr_sigma_raw->SetMarkerStyle(22);
    gr_sigma_raw->SetMarkerColor(kRed);
    gr_sigma_raw->SetLineColor(kRed);
    gr_sigma_raw->SetMarkerSize(1.2);
    
    // Set titles and labels
    gr_sigma_reg->SetTitle("Peak Width vs p_{T};p_{T} [GeV];#sigma [GeV]");
    gr_sigma_reg->GetXaxis()->SetTitleSize(0.045);
    gr_sigma_reg->GetYaxis()->SetTitleSize(0.045);
    gr_sigma_reg->GetXaxis()->SetLabelSize(0.04);
    gr_sigma_reg->GetYaxis()->SetLabelSize(0.04);
    
    // Draw the plots
    gr_sigma_reg->Draw("AP");
    gr_sigma_raw->Draw("P SAME");
    
    // Add legend
    TLegend* leg_sigma = new TLegend(0.65, 0.75, 0.89, 0.89);
    leg_sigma->AddEntry(gr_sigma_reg, "Regressed Mass", "p");
    leg_sigma->AddEntry(gr_sigma_raw, "Raw Mass", "p");
    leg_sigma->SetBorderSize(0);
    leg_sigma->SetFillStyle(0);
    leg_sigma->Draw();
    
    // Save the sigma vs pT plot
    c_sigma->SaveAs("sigma_vs_pt.png");
    
    cout << "Sigma vs pT plot saved as sigma_vs_pt.png" << endl;
    
    // Clean up
    file->Close();
    cout << "Plot saved as mass_projections_with_fits.png" << endl;
}