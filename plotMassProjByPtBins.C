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
#include <cmath>

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
    double yMaxValues[] = {600, 400, 300, 220, 140, 25};
    
    // Define rebin factors for each pT bin (regressed and raw can have different rebinning)
    int rebinFactor_reg[Nbins] = {1,1,2,2,2,2}; // Customize rebin factor for regressed histograms
    int rebinFactor_raw[Nbins] = {1,1,2,2,2,2}; // Customize rebin factor for raw histograms
    
    // Fit parameters for regressed SC
    double means_reg[Nbins] = {3.035, 3.02, 3.04, 3.0, 3.0743, 3.0916};
    double meanInf_reg[Nbins] = {3.01, 2.9, 3.0, 2.85, 2.9, 2.9};
    double meanSup_reg[Nbins] = {3.2, 3.2, 3.2, 3.2, 3.2, 3.2};
    
    double sigmas_reg[Nbins] = {0.16, 0.1602, 0.13, 0.14, 0.1152, 0.12};
    double sigmaInf_reg[Nbins] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    double sigmaSup_reg[Nbins] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3};
    
    double alphaLs_reg[Nbins] = {0.92, 0.92, 0.92, 0.92, 0.92, 0.92};
    double alphaLInf_reg[Nbins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double alphaLSup_reg[Nbins] = {15.0, 15.0, 15.0, 15.0, 15.0, 15.0};
    
    double nLs_reg[Nbins] = {16.2, 16.2, 16.2, 16.2, 16.2, 16.2};
    double nLInf_reg[Nbins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double nLSup_reg[Nbins] = {200.0, 200.0, 200.0, 200.0, 200.0, 200.0};
    
    double alphaRs_reg[Nbins] = {1.67, 1.67, 1.67, 1.67, 1.67, 1.67};
    double alphaRInf_reg[Nbins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double alphaRSup_reg[Nbins] = {15.0, 15.0, 15.0, 15.0, 15.0, 15.0};
    
    double nRs_reg[Nbins] = {4.59, 4.59, 4.59, 4.59, 4.59, 4.59};
    double nRInf_reg[Nbins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double nRSup_reg[Nbins] = {100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
    
    // Fit parameters for raw SC (adjust these values as needed)
    double means_raw[Nbins] = {2.4, 2.4, 2.5, 2.6, 2.7, 2.8};
    double meanInf_raw[Nbins] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    double meanSup_raw[Nbins] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
    
    double sigmas_raw[Nbins] = {0.2, 0.15, 0.15, 0.12, 0.1, 0.1};
    double sigmaInf_raw[Nbins] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    double sigmaSup_raw[Nbins] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    
    double alphaLs_raw[Nbins] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    double alphaLInf_raw[Nbins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double alphaLSup_raw[Nbins] = {15.0, 15.0, 15.0, 15.0, 15.0, 15.0};
    
    double nLs_raw[Nbins] = {16.0, 16.0, 16.0, 16.0, 16.0, 16.0};
    double nLInf_raw[Nbins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double nLSup_raw[Nbins] = {200.0, 200.0, 200.0, 200.0, 200.0, 200.0};
    
    double alphaRs_raw[Nbins] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
    double alphaRInf_raw[Nbins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double alphaRSup_raw[Nbins] = {15.0, 15.0, 15.0, 15.0, 15.0, 15.0};
    
    double nRs_raw[Nbins] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
    double nRInf_raw[Nbins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double nRSup_raw[Nbins] = {100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
    
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
    double means_fitted_reg[Nbins], meanErrors_reg[Nbins];
    double means_fitted_raw[Nbins], meanErrors_raw[Nbins];
    double resolution_reg[Nbins], resolutionErrors_reg[Nbins];
    double resolution_raw[Nbins], resolutionErrors_raw[Nbins];
    
    // Create individual canvases for each pt bin
    vector<TCanvas*> canvases;
    
    // Create summary canvas for all bins
    TCanvas* c_summary = new TCanvas("c_summary", "Mass Projections - All pT Bins", 1800, 1200);
    c_summary->Divide(3, 2); // 3 columns x 2 rows for 6 bins
    
    // Loop over pT bins
    for (int i = 1; i <= nbinsPt; i++) {
        // Create individual canvas for this pt bin
        TCanvas* c = new TCanvas(Form("c_mass_proj_bin_%d", i), Form("Mass Projections pT Bin %d", i), 1200, 800);
        canvases.push_back(c);
        
        // Switch to individual canvas for fitting and plotting
        c->cd();
        
        // Get projections
        TH1D* proj = (TH1D*)file->Get(Form("proj_bin_%d", i));
        TH1D* proj_rawSC = (TH1D*)file->Get(Form("RawSC_bin_%d", i));
        
        if (!proj || !proj_rawSC) {
            cout << "Missing projection for bin " << i << endl;
            continue;
        }
        
        // Apply customizable rebin factors
        if (rebinFactor_reg[i-1] > 1) {
            proj->Rebin(rebinFactor_reg[i-1]);
            cout << "Rebinned regressed histogram by factor " << rebinFactor_reg[i-1] << " for bin " << i << endl;
        }
        if (rebinFactor_raw[i-1] > 1) {
            proj_rawSC->Rebin(rebinFactor_raw[i-1]);
            cout << "Rebinned raw histogram by factor " << rebinFactor_raw[i-1] << " for bin " << i << endl;
        }
        
        proj->SetStats(kFALSE);
        proj_rawSC->SetStats(kFALSE);
        
        // Set different colors and styles
        proj->SetLineColor(kBlue);
        proj->SetLineWidth(2);
        proj_rawSC->SetLineColor(kRed);
        proj_rawSC->SetLineWidth(2);
        
        // Fit regressed SC histogram with customizable parameter ranges
        RooRealVar mass_reg("mass_reg", "m(e^{+}e^{-})", massMin_reg[i-1], massMax_reg[i-1]);
        RooDataHist data_reg("data_reg", "Dataset from histogram", mass_reg, proj);
        
        RooRealVar cb_mean_reg(Form("cb_mean_reg_%d", i), "Mean of CB", means_reg[i-1], meanInf_reg[i-1], meanSup_reg[i-1]);
        RooRealVar cb_sigma_reg(Form("cb_sigma_reg_%d", i), "Sigma of CB", sigmas_reg[i-1], sigmaInf_reg[i-1], sigmaSup_reg[i-1]);
        RooRealVar cb_alphaL_reg(Form("cb_alphaL_reg_%d", i), "AlphaL of CB", alphaLs_reg[i-1], alphaLInf_reg[i-1], alphaLSup_reg[i-1]);
        RooRealVar cb_nL_reg(Form("cb_nL_reg_%d", i), "nL of CB", nLs_reg[i-1], nLInf_reg[i-1], nLSup_reg[i-1]);
        RooRealVar cb_alphaR_reg(Form("cb_alphaR_reg_%d", i), "AlphaR of CB", alphaRs_reg[i-1], alphaRInf_reg[i-1], alphaRSup_reg[i-1]);
        RooRealVar cb_nR_reg(Form("cb_nR_reg_%d", i), "nR of CB", nRs_reg[i-1], nRInf_reg[i-1], nRSup_reg[i-1]);
        
        RooCrystalBall cb_reg(Form("cb_reg_%d", i), "Asymmetric Crystal Ball", mass_reg, cb_mean_reg, cb_sigma_reg, cb_alphaL_reg, cb_nL_reg, cb_alphaR_reg, cb_nR_reg);
        
        // Fit raw SC histogram with customizable parameter ranges
        RooRealVar mass_raw("mass_raw", "m(e^{+}e^{-})", massMin_raw[i-1], massMax_raw[i-1]);
        RooDataHist data_raw("data_raw", "Dataset from histogram", mass_raw, proj_rawSC);
        
        RooRealVar cb_mean_raw(Form("cb_mean_raw_%d", i), "Mean of CB", means_raw[i-1], meanInf_raw[i-1], meanSup_raw[i-1]);
        RooRealVar cb_sigma_raw(Form("cb_sigma_raw_%d", i), "Sigma of CB", sigmas_raw[i-1], sigmaInf_raw[i-1], sigmaSup_raw[i-1]);
        RooRealVar cb_alphaL_raw(Form("cb_alphaL_raw_%d", i), "AlphaL of CB", alphaLs_raw[i-1], alphaLInf_raw[i-1], alphaLSup_raw[i-1]);
        RooRealVar cb_nL_raw(Form("cb_nL_raw_%d", i), "nL of CB", nLs_raw[i-1], nLInf_raw[i-1], nLSup_raw[i-1]);
        RooRealVar cb_alphaR_raw(Form("cb_alphaR_raw_%d", i), "AlphaR of CB", alphaRs_raw[i-1], alphaRInf_raw[i-1], alphaRSup_raw[i-1]);
        RooRealVar cb_nR_raw(Form("cb_nR_raw_%d", i), "nR of CB", nRs_raw[i-1], nRInf_raw[i-1], nRSup_raw[i-1]);
        
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
        RooFitResult *fit_result_reg = cb_reg.fitTo(data_reg, RooFit::Save(), RooFit::PrintLevel(-1));
        RooFitResult *fit_result_raw = cb_raw.fitTo(data_raw, RooFit::Save(), RooFit::PrintLevel(-1));

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
        
        // Create TF1 functions for the Crystal Ball fit curves using fitted parameters
        // Crystal Ball function: A * exp(-0.5*((x-mean)/sigma)^2) for |t| < alpha
        //                       A * (n/alpha)^n * exp(-0.5*alpha^2) * (n/alpha - alpha - t)^(-n) for t >= alpha
        // where t = (x-mean)/sigma
        
        // Define Crystal Ball function for regressed fit
        TF1* f_reg = new TF1(Form("f_reg_%d", i), 
            [=](double *x, double *par) {
                double t = (x[0] - par[1]) / par[2]; // (x - mean) / sigma
                double abs_alpha_L = TMath::Abs(par[3]); // |alphaL|
                double abs_alpha_R = TMath::Abs(par[5]); // |alphaR|
                double n_L = par[4]; // nL
                double n_R = par[6]; // nR
                
                if (t < -abs_alpha_L) {
                    // Left tail
                    double A_L = TMath::Power(n_L/abs_alpha_L, n_L) * TMath::Exp(-0.5*abs_alpha_L*abs_alpha_L);
                    double B_L = n_L/abs_alpha_L - abs_alpha_L;
                    return par[0] * A_L * TMath::Power(B_L - t, -n_L);
                } else if (t > abs_alpha_R) {
                    // Right tail
                    double A_R = TMath::Power(n_R/abs_alpha_R, n_R) * TMath::Exp(-0.5*abs_alpha_R*abs_alpha_R);
                    double B_R = n_R/abs_alpha_R - abs_alpha_R;
                    return par[0] * A_R * TMath::Power(B_R + t, -n_R);
                } else {
                    // Gaussian core
                    return par[0] * TMath::Exp(-0.5*t*t);
                }
            }, massMin_reg[i-1], massMax_reg[i-1], 7);
        
        f_reg->SetParameters(proj->GetMaximum(), mean_reg_val, sigma_reg_val, alphaL_reg_val, nL_reg_val, alphaR_reg_val, nR_reg_val);
        f_reg->SetLineColor(kBlue+2);
        f_reg->SetLineWidth(3);
        f_reg->SetLineStyle(2);
        
        // Define Crystal Ball function for raw fit
        TF1* f_raw = new TF1(Form("f_raw_%d", i), 
            [=](double *x, double *par) {
                double t = (x[0] - par[1]) / par[2]; // (x - mean) / sigma
                double abs_alpha_L = TMath::Abs(par[3]); // |alphaL|
                double abs_alpha_R = TMath::Abs(par[5]); // |alphaR|
                double n_L = par[4]; // nL
                double n_R = par[6]; // nR
                
                if (t < -abs_alpha_L) {
                    // Left tail
                    double A_L = TMath::Power(n_L/abs_alpha_L, n_L) * TMath::Exp(-0.5*abs_alpha_L*abs_alpha_L);
                    double B_L = n_L/abs_alpha_L - abs_alpha_L;
                    return par[0] * A_L * TMath::Power(B_L - t, -n_L);
                } else if (t > abs_alpha_R) {
                    // Right tail
                    double A_R = TMath::Power(n_R/abs_alpha_R, n_R) * TMath::Exp(-0.5*abs_alpha_R*abs_alpha_R);
                    double B_R = n_R/abs_alpha_R - abs_alpha_R;
                    return par[0] * A_R * TMath::Power(B_R + t, -n_R);
                } else {
                    // Gaussian core
                    return par[0] * TMath::Exp(-0.5*t*t);
                }
            }, massMin_raw[i-1], massMax_raw[i-1], 7);
        
        f_raw->SetParameters(proj_rawSC->GetMaximum(), mean_raw_val, sigma_raw_val, alphaL_raw_val, nL_raw_val, alphaR_raw_val, nR_raw_val);
        f_raw->SetLineColor(kRed+2);
        f_raw->SetLineWidth(3);
        f_raw->SetLineStyle(2);
        
        // Draw the fit functions
        f_reg->Draw("SAME");
        f_raw->Draw("SAME");
        
        // Print fit results and store sigma values
        cout << "=== BIN " << i << " FIT RESULTS ===" << endl;
        if (fit_result_reg && fit_result_reg->status() == 0) {
            cout << "Regressed fit CONVERGED for bin " << i << endl;
            cout << "  Mean: " << cb_mean_reg.getVal() << " ± " << cb_mean_reg.getError() << endl;
            cout << "  Sigma: " << cb_sigma_reg.getVal() << " ± " << cb_sigma_reg.getError() << endl;
            cout << "  Fit status: " << fit_result_reg->status() << " (0 = success)" << endl;
            sigmas_fitted_reg[i-1] = cb_sigma_reg.getVal();
            sigmaErrors_reg[i-1] = cb_sigma_reg.getError();
            means_fitted_reg[i-1] = cb_mean_reg.getVal();
            meanErrors_reg[i-1] = cb_mean_reg.getError();
            // Calculate resolution = sigma/mu and its error
            resolution_reg[i-1] = cb_sigma_reg.getVal() / cb_mean_reg.getVal();
            // Error propagation: d(sigma/mu) = sqrt((dsigma/mu)^2 + (sigma*dmu/mu^2)^2)
            double rel_sigma_err = cb_sigma_reg.getError() / cb_sigma_reg.getVal();
            double rel_mean_err = cb_mean_reg.getError() / cb_mean_reg.getVal();
            resolutionErrors_reg[i-1] = resolution_reg[i-1] * sqrt(rel_sigma_err*rel_sigma_err + rel_mean_err*rel_mean_err);
        } else {
            cout << "Regressed fit FAILED for bin " << i << endl;
            cout << "  Mean: " << cb_mean_reg.getVal() << " ± " << cb_mean_reg.getError() << endl;
            cout << "  Sigma: " << cb_sigma_reg.getVal() << " ± " << cb_sigma_reg.getError() << endl;
            if (fit_result_reg) {
                cout << "  Fit status: " << fit_result_reg->status() << " (0 = success, others = failure)" << endl;
            } else {
                cout << "  Fit result is NULL!" << endl;
            }
            sigmas_fitted_reg[i-1] = cb_sigma_reg.getVal(); // Store value even if fit failed
            sigmaErrors_reg[i-1] = cb_sigma_reg.getError(); // Store error even if fit failed
            means_fitted_reg[i-1] = cb_mean_reg.getVal();
            meanErrors_reg[i-1] = cb_mean_reg.getError();
            // Calculate resolution even for failed fits
            resolution_reg[i-1] = cb_sigma_reg.getVal() / cb_mean_reg.getVal();
            double rel_sigma_err = cb_sigma_reg.getError() / cb_sigma_reg.getVal();
            double rel_mean_err = cb_mean_reg.getError() / cb_mean_reg.getVal();
            resolutionErrors_reg[i-1] = resolution_reg[i-1] * sqrt(rel_sigma_err*rel_sigma_err + rel_mean_err*rel_mean_err);
        }
      
        if (fit_result_raw && fit_result_raw->status() == 0) {
            cout << "Raw fit CONVERGED for bin " << i << endl;
            cout << "  Mean: " << cb_mean_raw.getVal() << " ± " << cb_mean_raw.getError() << endl;
            cout << "  Sigma: " << cb_sigma_raw.getVal() << " ± " << cb_sigma_raw.getError() << endl;
            cout << "  Fit status: " << fit_result_raw->status() << " (0 = success)" << endl;
            sigmas_fitted_raw[i-1] = cb_sigma_raw.getVal();
            sigmaErrors_raw[i-1] = cb_sigma_raw.getError();
            means_fitted_raw[i-1] = cb_mean_raw.getVal();
            meanErrors_raw[i-1] = cb_mean_raw.getError();
            // Calculate resolution = sigma/mu and its error
            resolution_raw[i-1] = cb_sigma_raw.getVal() / cb_mean_raw.getVal();
            // Error propagation: d(sigma/mu) = sqrt((dsigma/mu)^2 + (sigma*dmu/mu^2)^2)
            double rel_sigma_err = cb_sigma_raw.getError() / cb_sigma_raw.getVal();
            double rel_mean_err = cb_mean_raw.getError() / cb_mean_raw.getVal();
            resolutionErrors_raw[i-1] = resolution_raw[i-1] * sqrt(rel_sigma_err*rel_sigma_err + rel_mean_err*rel_mean_err);
        } else {
            cout << "Raw fit FAILED for bin " << i << endl;
            cout << "  Mean: " << cb_mean_raw.getVal() << " ± " << cb_mean_raw.getError() << endl;
            cout << "  Sigma: " << cb_sigma_raw.getVal() << " ± " << cb_sigma_raw.getError() << endl;
            if (fit_result_raw) {
                cout << "  Fit status: " << fit_result_raw->status() << " (0 = success, others = failure)" << endl;
            } else {
                cout << "  Fit result is NULL!" << endl;
            }
            sigmas_fitted_raw[i-1] = cb_sigma_raw.getVal(); // Store value even if fit failed
            sigmaErrors_raw[i-1] = cb_sigma_raw.getError(); // Store error even if fit failed
            means_fitted_raw[i-1] = cb_mean_raw.getVal();
            meanErrors_raw[i-1] = cb_mean_raw.getError();
            // Calculate resolution even for failed fits
            resolution_raw[i-1] = cb_sigma_raw.getVal() / cb_mean_raw.getVal();
            double rel_sigma_err = cb_sigma_raw.getError() / cb_sigma_raw.getVal();
            double rel_mean_err = cb_mean_raw.getError() / cb_mean_raw.getVal();
            resolutionErrors_raw[i-1] = resolution_raw[i-1] * sqrt(rel_sigma_err*rel_sigma_err + rel_mean_err*rel_mean_err);
        }
        cout << "=================================" << endl;

        
        // Add fit parameter box with both regressed and raw results
        TPaveText* fitBox = new TPaveText(0.62, 0.7, 0.9, 0.9, "NDC");
        fitBox->SetFillColor(0);
        fitBox->SetBorderSize(1);
        fitBox->SetTextSize(0.025);
        fitBox->SetTextAlign(12); // Left aligned
        
        // Add regressed fit parameters
        if (fit_result_reg && fit_result_reg->status() == 0) {
            fitBox->AddText("#color[4]{Regressed:}");
            fitBox->AddText(Form("#color[4]{#mu = %.4f #pm %.4f}", cb_mean_reg.getVal(), cb_mean_reg.getError()));
            fitBox->AddText(Form("#color[4]{#sigma = %.4f #pm %.4f}", cb_sigma_reg.getVal(), cb_sigma_reg.getError()));
        } else {
            fitBox->AddText("#color[4]{Regressed:}");
            fitBox->AddText(Form("#color[4]{#mu = %.4f #pm %.4f}", cb_mean_reg.getVal(), cb_mean_reg.getError()));
            fitBox->AddText(Form("#color[4]{#sigma = %.4f #pm %.4f}", cb_sigma_reg.getVal(), cb_sigma_reg.getError()));
        }
        
        fitBox->AddText(""); // Empty line for separation
        
        // Add raw fit parameters
        if (fit_result_raw && fit_result_raw->status() == 0) {
            fitBox->AddText("#color[2]{Raw:}");
            fitBox->AddText(Form("#color[2]{#mu = %.4f #pm %.4f}", cb_mean_raw.getVal(), cb_mean_raw.getError()));
            fitBox->AddText(Form("#color[2]{#sigma = %.4f #pm %.4f}", cb_sigma_raw.getVal(), cb_sigma_raw.getError()));
        } else {
            fitBox->AddText("#color[2]{Raw:}");
            fitBox->AddText(Form("#color[2]{#mu = %.4f #pm %.4f}", cb_mean_raw.getVal(), cb_mean_raw.getError()));
            fitBox->AddText(Form("#color[2]{#sigma = %.4f #pm %.4f}", cb_sigma_raw.getVal(), cb_sigma_raw.getError()));
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
        
        // Save individual canvas
        c->Update();
        c->SaveAs(Form("mass_projection_bin_%d.png", i));
        cout << "Saved: mass_projection_bin_" << i << ".png" << endl;
        
        // Now create the same plot on the summary canvas
        c_summary->cd(i); // Switch to the i-th pad on summary canvas
        
        // Clone the frame histogram for summary canvas
        TH1F* hFrame_summary = (TH1F*)hFrame->Clone(Form("hFrame_summary_%d", i));
        hFrame_summary->SetTitle(Form("p_{T} bin %d", i));
        hFrame_summary->Draw();
        
        // Clone and draw the data points
        TGraphErrors* gr_raw_summary = (TGraphErrors*)gr_raw->Clone(Form("gr_raw_summary_%d", i));
        TGraphErrors* gr_reg_summary = (TGraphErrors*)gr_reg->Clone(Form("gr_reg_summary_%d", i));
        gr_raw_summary->Draw("P SAME");
        gr_reg_summary->Draw("P SAME");
        
        // Clone and draw the fit functions
        TF1* f_reg_summary = (TF1*)f_reg->Clone(Form("f_reg_summary_%d", i));
        TF1* f_raw_summary = (TF1*)f_raw->Clone(Form("f_raw_summary_%d", i));
        f_reg_summary->Draw("SAME");
        f_raw_summary->Draw("SAME");
        
        // Clone and draw the fit parameter box
        TPaveText* fitBox_summary = (TPaveText*)fitBox->Clone(Form("fitBox_summary_%d", i));
        fitBox_summary->SetTextSize(0.04); // Slightly larger text for summary canvas
        fitBox_summary->Draw();
        
        // Add a smaller legend for summary canvas
        TLegend* leg_summary = new TLegend(0.15, 0.75, 0.45, 0.89);
        leg_summary->AddEntry(gr_reg_summary, "Regressed", "p");
        leg_summary->AddEntry(gr_raw_summary, "Raw", "p");
        leg_summary->AddEntry(f_reg_summary, "Reg. Fit", "l");
        leg_summary->AddEntry(f_raw_summary, "Raw Fit", "l");
        leg_summary->SetBorderSize(0);
        leg_summary->SetFillStyle(0);
        leg_summary->SetTextSize(0.05);
        leg_summary->Draw();
        
        // Update the summary pad
        c_summary->cd(i)->Update();
        
        // Delete the canvas to free memory
        delete c;
        canvases[i-1] = nullptr; // Set pointer to null for safety
        
        // Clean up RooPlot frames
    }
    
    // Save summary canvas with all bins
    c_summary->Update();
    c_summary->SaveAs("mass_projections_all_bins_summary.png");
    cout << "Saved summary canvas: mass_projections_all_bins_summary.png" << endl;
    
    // Delete the summary canvas
    delete c_summary;
    
    // Note: Individual canvases have been deleted after saving
    // Summary canvas created and saved with all 6 bins
    cout << "All individual canvases have been saved and deleted" << endl;
    
    // Create resolution vs pT plot
    TCanvas* c_resolution = new TCanvas("c_resolution", "Resolution vs pT", 800, 600);
    
    // Create TGraphErrors for resolution vs pT
    TGraphErrors* gr_resolution_reg = new TGraphErrors(nbinsPt, ptBinCenters, resolution_reg, ptBinHalfWidths, resolutionErrors_reg);
    TGraphErrors* gr_resolution_raw = new TGraphErrors(nbinsPt, ptBinCenters, resolution_raw, ptBinHalfWidths, resolutionErrors_raw);
    
    // Set styles
    gr_resolution_reg->SetMarkerStyle(20);
    gr_resolution_reg->SetMarkerColor(kBlue);
    gr_resolution_reg->SetLineColor(kBlue);
    gr_resolution_reg->SetMarkerSize(1.2);
    
    gr_resolution_raw->SetMarkerStyle(22);
    gr_resolution_raw->SetMarkerColor(kRed);
    gr_resolution_raw->SetLineColor(kRed);
    gr_resolution_raw->SetMarkerSize(1.2);
    
    // Set titles and labels
    gr_resolution_reg->SetTitle("Mass Resolution vs p_{T};p_{T} [GeV];#sigma/#mu");
    gr_resolution_reg->GetXaxis()->SetTitleSize(0.045);
    gr_resolution_reg->GetYaxis()->SetTitleSize(0.045);
    gr_resolution_reg->GetXaxis()->SetLabelSize(0.04);
    gr_resolution_reg->GetYaxis()->SetLabelSize(0.04);
    
    // Draw the plots
    gr_resolution_reg->Draw("AP");
    
    // Set y-axis range for better visibility
    gr_resolution_reg->GetYaxis()->SetRangeUser(0.0, 0.1); // Adjust range as needed
    
    gr_resolution_raw->Draw("P SAME");
    
    // Add legend
    TLegend* leg_resolution = new TLegend(0.65, 0.75, 0.89, 0.89);
    leg_resolution->AddEntry(gr_resolution_reg, "Regressed Mass", "p");
    leg_resolution->AddEntry(gr_resolution_raw, "Raw Mass", "p");
    leg_resolution->SetBorderSize(0);
    leg_resolution->SetFillStyle(0);
    leg_resolution->Draw();
    
    // Save the resolution vs pT plot
    c_resolution->SaveAs("resolution_vs_pt.png");
    
    cout << "Resolution vs pT plot saved as resolution_vs_pt.png" << endl;
    
    // Delete the resolution canvas
    delete c_resolution;
    
    // Clean up
    file->Close();
    cout << "Individual plots saved as mass_projection_bin_X.png" << endl;
    cout << "Resolution vs pT plot saved as resolution_vs_pt.png" << endl;
    cout << "All canvases have been properly deleted after saving" << endl;
}