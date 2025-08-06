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
#include <TSystem.h>

#define Nbins 6
#define Nsigma 1.0  // Change this to adjust constraint width (1.0 = 1 sigma, 2.0 = 2 sigma, etc.)


// Define CMS colors
Color_t gialloCMS = TColor::GetColor("#ffcc00");
Color_t violaCMS = TColor::GetColor("#660099");
Color_t rossoCMS = TColor::GetColor("#cc0000");
int bluCMS = TColor::GetColor("#5790FC");

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....Funzione globale per scritta MC....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void WriteSimulation(){
    // Draw "CMS" in bold
    TLatex cmsLabel;
    cmsLabel.SetNDC();
    cmsLabel.SetTextFont(62);  // Helvetica Bold
    cmsLabel.SetTextSize(0.05);
    cmsLabel.SetTextAlign(11);
    cmsLabel.DrawLatex(0.12, 0.93, "CMS");
    
    // Draw "Simulation" in italic
    TLatex simLabel;
    simLabel.SetNDC();
    simLabel.SetTextFont(52);  // Italic font
    simLabel.SetTextSize(0.05);
    simLabel.SetTextAlign(11);
    simLabel.DrawLatex(0.21, 0.93, "Simulation");
}

void WritePreliminary(){
    // Add CMS Preliminary label
    TLatex cmsLabel;
    cmsLabel.SetNDC();
    cmsLabel.SetTextFont(62);
    cmsLabel.SetTextSize(0.05);
    cmsLabel.SetTextAlign(11);
    cmsLabel.DrawLatex(0.14, 0.92, "CMS");

    // Add Preliminary in italics
    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextFont(52);
    prelimLabel.SetTextSize(0.05);
    prelimLabel.SetTextAlign(11);
    prelimLabel.DrawLatex(0.22, 0.92, "Preliminary");

    // Draw luminosity/energy information
    TLatex lumiLabel;
    lumiLabel.SetNDC();
    lumiLabel.SetTextFont(42);
    lumiLabel.SetTextSize(0.045);
    lumiLabel.SetTextAlign(31);
    lumiLabel.DrawLatex(0.9, 0.92, "38.01 fb^{-1} (2022, 13.7 TeV)");
}

void PlotBackgroundFit(RooRealVar& mass, RooDataHist& data, RooAddPdf& background, int i) {
    // Create a canvas to draw the background fit
    TCanvas* cBackgroundFit = new TCanvas(Form("cBackgroundFit_%d", i+1), "Background Fit", 800, 600);
    gPad->SetLeftMargin(0.13);

    // Create a frame for the mass variable to hold the plot
    RooPlot* frame = mass.frame();
    
    // Plot the data points on the frame
    data.plotOn(frame, RooFit::Name("Data"), RooFit::MarkerStyle(20), RooFit::MarkerSize(0.7));

    // Plot the background fit function
    background.plotOn(frame, RooFit::Name("BackgroundFit"), RooFit::LineColor(kMagenta), RooFit::LineStyle(kDashed));

    // Add labels and titles
    frame->SetTitle(Form("Background Fit for Raw SC Pt Bin %d", i+1));
    frame->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle("Events");

    // Draw the frame on the canvas
    frame->Draw();

    // Add legend for clarity
    TLegend* legend = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend->AddEntry(frame->findObject("Data"), "Data", "p");
    legend->AddEntry(frame->findObject("BackgroundFit"), "Background Fit", "l");
    legend->Draw();

    // Extract parameters from the model for printing
    RooRealVar* mu = (RooRealVar*)background.getVariables()->find(Form("gauss_mu_%d", i+1));
    RooRealVar* sigma = (RooRealVar*)background.getVariables()->find(Form("gauss_sigma_%d", i+1));

    // Optional: add text annotations for fit parameters
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.15, 0.85, Form("Background Fit for Pt bin %d", i+1));
    if(mu && sigma) {
        latex.DrawLatex(0.15, 0.8, Form("#mu = %.4f +/- %.4f, #sigma = %.4f +/- %.4f", mu->getVal(), mu->getError(), sigma->getVal(), sigma->getError()));
    }

    WritePreliminary();

    // Save the canvas as a file
    cBackgroundFit->SaveAs(Form("PlotConID2022/RawScaleandSmearing/Pt_bin%d/BackgroundFit_RawSC_Pt%d.png", i+1, i+1));
    delete cBackgroundFit;
}

void PlotDataFit(RooRealVar& mass, RooDataHist& data, RooAddPdf& model, RooAddPdf& background, RooCrystalBall& crystal, double leftlim, double rightlim, int i, int fitstatus) {

    // Create canvas for plotting
    TCanvas *c = new TCanvas(Form("c_bin_%d", i+1), Form("Bin %d", i+1), 950, 700);
    gPad->SetLeftMargin(0.13);

    // Frame for mass with restricted range
    RooPlot* frame = mass.frame(RooFit::Range(leftlim, rightlim));
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("Invariant Mass [GeV]");
    frame->GetYaxis()->SetTitle("Events");

    // Plot the data points
    data.plotOn(frame, RooFit::Name("Data"));

    // Plot the background component of the model with dotted style and CMS yellow color
    model.plotOn(frame, RooFit::Components(background), RooFit::LineStyle(kDotted), RooFit::LineColor(gialloCMS), RooFit::LineWidth(5), RooFit::Name("Background"));

    // Plot the signal (Crystal Ball) component with dashed style and CMS purple color
    model.plotOn(frame, RooFit::Components(crystal), RooFit::LineStyle(kDashed), RooFit::LineColor(violaCMS), RooFit::LineWidth(5), RooFit::Name("CrystalBall"));

    // Plot the full model with solid line and CMS red color
    model.plotOn(frame, RooFit::LineColor(rossoCMS), RooFit::LineWidth(5), RooFit::Name("CombinedFit"));

    // Draw frame on canvas
    frame->Draw();

    // Add legend to clarify components
    TLegend *legend = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend->SetBorderSize(0);
    legend->AddEntry(frame->findObject("Data"), "Data", "p");
    legend->AddEntry(frame->findObject("Background"), "Background Fit", "l");
    legend->AddEntry(frame->findObject("CrystalBall"), "Signal (Crystal Ball)", "l");
    legend->AddEntry(frame->findObject("CombinedFit"), "Combined Fit", "l");
    legend->Draw();

    // Extract parameters from the model for printing
    RooRealVar* mu = (RooRealVar*)model.getVariables()->find(Form("mu_cb_%d", i+1));
    RooRealVar* sigma = (RooRealVar*)model.getVariables()->find(Form("sigma_cb_%d", i+1));
    double chi2 = frame->chiSquare();

    // Optional: Add annotations for fit information
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    if(mu && sigma) {
        latex.DrawLatex(0.17, 0.8, Form("#mu = %.4f #pm %.4f GeV", mu->getVal(), mu->getError()));
        latex.DrawLatex(0.17, 0.75, Form("#sigma = %.4f #pm %.4f GeV", sigma->getVal(), sigma->getError())); 
    }

    WritePreliminary();

    // Save plot as image
    c->SaveAs(Form("PlotConID2022/RawScaleandSmearing/Pt_bin%d/DataFit_RawSC_Pt%d.png", i+1, i+1));
    delete c;
}

void FitRawScaleAndSmearing() {
    
    // Create directories for output plots
    for(int i = 0; i < Nbins; i++) {
        gSystem->Exec(Form("mkdir -p PlotConID2022/RawScaleandSmearing/Pt_bin%d", i+1));
    }

    // Open MC and Data files
    TFile *fileMC = TFile::Open("outputHistograms_MC.root");
    TFile *fileData = TFile::Open("outputHistograms_DATA_partF.root");

    // Get the 2D histograms
    TH2D* h2d_MC = (TH2D*)fileMC->Get("h_rawSC_vs_pt_reweighted");
    TH2D* h2d_Data = (TH2D*)fileData->Get("h_rawSC_vs_pt");
    
    if (!h2d_MC) {
        std::cerr << "MC histogram h_rawSC_vs_pt_reweighted not found!" << std::endl;
        return;
    }
    if (!h2d_Data) {
        std::cerr << "Data histogram h_rawSC_vs_pt not found!" << std::endl;
        return;
    }

    // Arrays for MC fit parameters
    double means[Nbins], sigmas[Nbins], alphaLs[Nbins], nLs[Nbins], alphaRs[Nbins], nRs[Nbins];
    double meanInf[Nbins], meanSup[Nbins];
    double mu_mc[Nbins], inc_mu_mc[Nbins], val_sigma_mc[Nbins], inc_sigma_mc[Nbins];
    double val_nL[Nbins], inc_nL[Nbins], val_alphL[Nbins], inc_alphL[Nbins];
    double val_nR[Nbins], inc_nR[Nbins], val_alphR[Nbins], inc_alphR[Nbins];

    // Initialize MC fit parameters
    for (int i = 0; i < Nbins; i++) {
        means[i] = 2.5;
        sigmas[i] = 0.14;
        alphaLs[i] = 0.92;
        nLs[i] = 16.2;
        alphaRs[i] = 1.67;
        nRs[i] = 4.59;
        meanInf[i] = means[i] - 2 * sigmas[i];
        meanSup[i] = means[i] + 3 * sigmas[i];
    }

    //customized initialization for monte carlo
    means[0] = 2.38;
    sigmas[0] = 0.177; 
    meanInf[0] = 2.2;
    meanSup[0] = 2.5;      

    means[1] = 2.44;
    sigmas[1] = 0.17;

    means[2] = 2.53;
    sigmas[2] = 0.1416;

    means[3] = 2.617;
    sigmas[3] = 0.11;

    means[4] = 2.72;
    sigmas[4] = 0.11;
    meanInf[4] = 2.6;
    meanSup[4] = 2.9;

    means[5] = 2.75;
    meanInf[5] = 2.7;
    meanSup[5] = 2.9;
    sigmas[5] = 0.08;
    

    // Mass range and rebin factors for each bin
    double massMin[Nbins] = {1.6, 1.9, 1.6, 1.8, 2.25, 2.2};
    double massMax[Nbins] = {3.0, 3.0, 3.1, 3.2, 3.1, 3.1};
    int rebin_factor[Nbins] = {2, 2, 1, 1, 2, 3};

    // Data fit parameters (inclusive in run number)
    double LeftLowLim_incl[Nbins] = {1.05, 1, 1, 1.1, 1.5, 1.95};
    double LeftUpLim_incl[Nbins] = {1.8, 1.9, 2, 2, 2.2, 2.5};
    double RightLowLim_incl[Nbins] = {2.97, 2.95, 2.95, 3, 3.05, 3.2};
    double RightUpLim_incl[Nbins] = {4.5, 4.5, 5, 4.5, 5, 4.8};

    // Gaussian parameters for background fit in data
    double gauss_mu_init_incl[Nbins] = {3.0, 3.3, 3, 3.3, 3.14, 3.3};
    double gauss_mu_low_incl[Nbins] = {2.95, 2.8, 2.8, 2.8, 3, 3};
    double gauss_mu_up_incl[Nbins] = {3.4, 3.2, 3.4, 3.4, 3.3, 3.3};
    double gauss_sigma_init_incl[Nbins] = {0.2, 0.2, 0.1, 0.1, 0.13, 0.15};
    double gauss_sigma_low_incl[Nbins] = {0.2, 0.15, 0.1, 0.1, 0.1, 0.1};
    double gauss_sigma_up_incl[Nbins] = {0.4, 0.3, 0.4, 0.4, 0.2, 0.2};

    //fraction of Gaussian in background fit
    double frac_gauss_init[Nbins] = {0.25, 0.25, 0.3, 0.3, 0.4, 0.3};

    // Crystal Ball parameters for data fit (personalized limits)
    double Mucb_init_incl[Nbins] = {2.44, 2.52, 2.54, 2.617, 2.7, 2.82};
    double Mucb_lowlim_incl[Nbins] = {2.38, 2.45, 2.4, 2.5, 2.65, 2.6};
    double Mucb_uplim_incl[Nbins] = {2.5, 2.55, 2.6, 2.7, 2.75, 3};
    double Sigmacb_init_incl[Nbins] = {0.1601, 0.133, 0.1399, 0.136, 0.1317, 0.12};
    double Sigmacb_lowlim_incl[Nbins] = {0.09, 0.12, 0.1, 0.11, 0.1, 0.09};
    double Sigmacb_uplim_incl[Nbins] = {0.25, 0.25, 0.2, 0.2, 0.18, 0.16};
    
    // Arrays to store fitted values for scale and smearing calculations
    double mu_data_fitted[Nbins], sigma_data_fitted[Nbins];
    double inc_mu_data_fitted[Nbins], inc_sigma_data_fitted[Nbins];
    
    // Bin centers and half-widths for Raw SC pT (similar to FitData.C but for Raw SC bins)
    double bincenters[Nbins] = {5.75, 8.25, 10, 12.5, 17, 30};  // Approximate bin centers for Raw SC pT
    double binhalfwidths[Nbins] = {1.75, 0.75, 1, 1.5, 3, 10};  // Approximate bin half-widths
    
    // Create TGraphErrors for scale and smearing vs pT
    TGraphErrors *graph_scale = new TGraphErrors(Nbins);
    TGraphErrors *graph_smearing = new TGraphErrors(Nbins);
    graph_scale->SetName("Graph_Scale_RawSC");
    graph_smearing->SetName("Graph_Smearing_RawSC");
    

    //===========================================
    // STEP 1: FIT MC PROJECTIONS WITH DOUBLE CRYSTAL BALL
    //===========================================
    std::cout << "=== STEP 1: FITTING MC PROJECTIONS ===" << std::endl;
    
    for (int i = 0; i < Nbins; i++) {
        std::cout << "Fitting MC projection for bin " << i+1 << std::endl;
        
        // Project the Y (mass) for pt bin i using bin i+1
        TH1D* hist_mc = h2d_MC->ProjectionY(Form("proj_mc_bin_%d", i+1), i+1, i+1);
        if (!hist_mc || hist_mc->GetEntries() == 0) {
            std::cerr << "Empty MC projection for pt bin " << i+1 << std::endl;
            continue;
        }
        hist_mc->Rebin(rebin_factor[i]);

        // Get mean and RMS for initialization
        means[i] = hist_mc->GetMean();
        sigmas[i] = hist_mc->GetRMS() / 2.0;

        // Define mass variable
        RooRealVar mass("mass", "m(e^{+}e^{-})", massMin[i], massMax[i]); 

        // Convert histogram to RooDataHist
        RooDataHist data_mc("data_mc", "MC Dataset", mass, hist_mc);

        // Define Crystal Ball parameters
        RooRealVar cb_mean(Form("cb_mean_%d", i+1), "Mean of CB", means[i], meanInf[i], meanSup[i]);
        RooRealVar cb_sigma(Form("cb_sigma_%d", i+1), "Sigma of CB", sigmas[i], 0.1 * sigmas[i], 2 * sigmas[i]);
        RooRealVar cb_alphaL(Form("cb_alphaL_%d", i+1), "AlphaL of CB", alphaLs[i], 0, 15.0);
        RooRealVar cb_nL(Form("cb_nL_%d", i+1), "nL of CB", nLs[i], 0, 500.0);
        RooRealVar cb_alphaR(Form("cb_alphaR_%d", i+1), "AlphaR of CB", alphaRs[i], 0, 15.0);
        RooRealVar cb_nR(Form("cb_nR_%d", i+1), "nR of CB", nRs[i], 0, 500.0);

        // Create asymmetric Crystal Ball 
        RooCrystalBall cb(Form("cb_%d", i+1), "Asymmetric Crystal Ball", mass, cb_mean, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR);

        // Perform fit
        RooFitResult *fit_result_simplex = cb.fitTo(data_mc, RooFit::Minimizer("Minuit2", "Simplex"), RooFit::MaxCalls(1000000), RooFit::Save(), RooFit::SumW2Error(true), RooFit::PrintLevel(-1));
        RooFitResult *fit_result = cb.fitTo(data_mc, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::MaxCalls(10000000), RooFit::Save(), RooFit::SumW2Error(true), RooFit::PrintLevel(-1));

        // Plot MC fit
        TCanvas *canvas_mc = new TCanvas(Form("canvas_mc_%d", i+1), Form("MC Fit for proj_bin_%d", i+1), 900, 700);
        RooPlot *frame_mc = mass.frame();
        frame_mc->SetTitle("");
        frame_mc->GetXaxis()->SetRangeUser(0, 5);
        data_mc.plotOn(frame_mc);
        cb.plotOn(frame_mc, RooFit::LineColor(bluCMS), RooFit::LineWidth(5));
        frame_mc->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");

        double chi2_mc = frame_mc->chiSquare();
        TPaveText *paveText_mc = new TPaveText(0.7, 0.75, 0.88, 0.88, "NDC");
        paveText_mc->AddText(Form("#chi^{2} = %.2f", chi2_mc));
        paveText_mc->AddText(Form("#mu = %.4f +/- %.4f", cb_mean.getVal(), cb_mean.getError()));
        paveText_mc->AddText(Form("#sigma = %.4f +/- %.4f", cb_sigma.getVal(), cb_sigma.getError()));
        paveText_mc->SetFillColor(0);
        paveText_mc->SetBorderSize(0);
        paveText_mc->SetShadowColor(0);
        frame_mc->addObject(paveText_mc);

        frame_mc->Draw();
        WriteSimulation();

        // Save MC fit plot
        canvas_mc->SaveAs(Form("PlotConID2022/RawScaleandSmearing/Pt_bin%d/MCFit_RawSC_Pt%d.png", i+1, i+1));
        delete canvas_mc;

        // Store MC fit results for data fit initialization
        mu_mc[i] = cb_mean.getVal();
        inc_mu_mc[i] = cb_mean.getError();
        val_sigma_mc[i] = cb_sigma.getVal();
        inc_sigma_mc[i] = cb_sigma.getError();
        val_nL[i] = cb_nL.getVal();
        inc_nL[i] = cb_nL.getError();
        val_alphL[i] = cb_alphaL.getVal();
        inc_alphL[i] = cb_alphaL.getError();
        val_nR[i] = cb_nR.getVal();
        inc_nR[i] = cb_nR.getError();
        val_alphR[i] = cb_alphaR.getVal();
        inc_alphR[i] = cb_alphaR.getError();

        std::cout << "MC Bin " << i+1 << ": mu = " << mu_mc[i] << " +/- " << inc_mu_mc[i] 
                  << ", sigma = " << val_sigma_mc[i] << " +/- " << inc_sigma_mc[i] << std::endl;
    }

    //===========================================
    // STEP 2: FIT DATA PROJECTIONS
    //===========================================
    std::cout << "\n=== STEP 2: FITTING DATA PROJECTIONS ===" << std::endl;

    for (int i = 0; i < Nbins; i++) {
        std::cout << "Fitting Data projection for bin " << i+1 << std::endl;

        // Project the Y (mass) for pt bin i using bin i+1
        TH1D* hist_data = h2d_Data->ProjectionY(Form("proj_data_bin_%d", i+1), i+1, i+1);
        if (!hist_data || hist_data->GetEntries() == 0) {
            std::cerr << "Empty Data projection for pt bin " << i+1 << std::endl;
            continue;
        }
       

        // Define mass variable for data fit
        RooRealVar mass_data("mass_data", "m(e^{+}e^{-})", 0, 6); 

        // Convert histogram to RooDataHist
        RooDataHist data("data", "Data Dataset", mass_data, hist_data);

        //--- BACKGROUND FIT FIRST ---
        // Define polynomial background parameters
        RooRealVar A(Form("A_%d", i+1), "4th deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
        RooRealVar B(Form("B_%d", i+1), "3rd deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
        RooRealVar C(Form("C_%d", i+1), "2nd deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
        RooRealVar D(Form("D_%d", i+1), "1st deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
        RooRealVar E(Form("E_%d", i+1), "0 deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
        RooPolynomial poly("poly", "Polynomial of 4th degree", mass_data, RooArgList(A, B, C, D, E));

        // Define Gaussian for background
        RooRealVar gauss_mu(Form("gauss_mu_%d", i+1), "Gaussian mean", gauss_mu_init_incl[i], gauss_mu_low_incl[i], gauss_mu_up_incl[i]);
        RooRealVar gauss_sigma(Form("gauss_sigma_%d", i+1), "Gaussian sigma", gauss_sigma_init_incl[i], gauss_sigma_low_incl[i], gauss_sigma_up_incl[i]);
        RooGaussian gauss("gauss", "Gaussian component", mass_data, gauss_mu, gauss_sigma);

        // Background fractions
        RooRealVar frac_gauss(Form("frac_gauss_%d", i+1), "Fraction of Gaussian", frac_gauss_init[i], 0.0, 1.0);
        RooAddPdf background("background", "Background model", RooArgList(poly, gauss), RooArgList(frac_gauss));

        // Fit background in sidebands
        mass_data.setRange("range1", LeftLowLim_incl[i], LeftUpLim_incl[i]);
        mass_data.setRange("range2", RightLowLim_incl[i], RightUpLim_incl[i]);

        RooFitResult *bg_fit_result = background.fitTo(data, RooFit::Range("range1,range2"), RooFit::Save(), RooFit::PrintLevel(-1));

        // Extract fitted background parameter values and errors
        double A_val = A.getVal(), A_err = A.getError();
        double B_val = B.getVal(), B_err = B.getError();
        double C_val = C.getVal(), C_err = C.getError();
        double D_val = D.getVal(), D_err = D.getError();
        double E_val = E.getVal(), E_err = E.getError();
        double gauss_mu_val = gauss_mu.getVal(), gauss_mu_err = gauss_mu.getError();
        double gauss_sigma_val = gauss_sigma.getVal(), gauss_sigma_err = gauss_sigma.getError();
        double frac_gauss_val = frac_gauss.getVal(), frac_gauss_err = frac_gauss.getError();

        // Plot background fit
        PlotBackgroundFit(mass_data, data, background, i);

        //--- SIGNAL + BACKGROUND FIT ---
        // Initialize Crystal Ball parameters with personalized values and limits
        RooRealVar mu_cb(Form("mu_cb_%d", i+1), "CB mean", Mucb_init_incl[i], Mucb_lowlim_incl[i], Mucb_uplim_incl[i]);
        RooRealVar sigma_cb(Form("sigma_cb_%d", i+1), "CB sigma", Sigmacb_init_incl[i], Sigmacb_lowlim_incl[i], Sigmacb_uplim_incl[i]);
        RooRealVar alphaL_cb(Form("alphaL_cb_%d", i+1), "CB alphaL", val_alphL[i], 0, 15.0);
        RooRealVar nL_cb(Form("nL_cb_%d", i+1), "CB nL", val_nL[i], 0, 500.0);
        RooRealVar alphaR_cb(Form("alphaR_cb_%d", i+1), "CB alphaR", val_alphR[i], 0, 15.0);
        RooRealVar nR_cb(Form("nR_cb_%d", i+1), "CB nR", val_nR[i], 0, 500.0);

        // Create Crystal Ball for signal
        RooCrystalBall crystal("crystal", "Crystal Ball signal", mass_data, mu_cb, sigma_cb, alphaL_cb, nL_cb, alphaR_cb, nR_cb);

        // Constrain background parameters within Nsigma from fitted values
        gauss_mu.setVal(gauss_mu_val);
        gauss_mu.setConstant(true);  // Keep Gaussian mean fixed
        gauss_sigma.setVal(gauss_sigma_val);
        gauss_sigma.setRange(gauss_sigma_val - Nsigma*gauss_sigma_err, gauss_sigma_val + Nsigma*gauss_sigma_err);
        A.setVal(A_val);
        A.setRange(A_val - Nsigma*A_err, A_val + Nsigma*A_err);
        B.setVal(B_val);
        B.setRange(B_val - Nsigma*B_err, B_val + Nsigma*B_err);
        C.setVal(C_val);
        C.setRange(C_val - Nsigma*C_err, C_val + Nsigma*C_err);
        D.setVal(D_val);
        D.setRange(D_val - Nsigma*D_err, D_val + Nsigma*D_err);
        E.setVal(E_val);
        E.setRange(E_val - Nsigma*E_err, E_val + Nsigma*E_err);
        frac_gauss.setVal(frac_gauss_val);
        frac_gauss.setRange(frac_gauss_val - Nsigma*frac_gauss_err, frac_gauss_val + Nsigma*frac_gauss_err);

        // Combined model fractions
        RooRealVar nsig(Form("nsig_%d", i+1), "Number of signal events", hist_data->Integral() * 0.1, 0, hist_data->Integral());
        RooRealVar nbkg(Form("nbkg_%d", i+1), "Number of background events", hist_data->Integral() * 0.9, 0, hist_data->Integral());

        // Combined model
        RooAddPdf model("model", "Signal + Background", RooArgList(crystal, background), RooArgList(nsig, nbkg));

        // Set fit range for signal + background fit
        mass_data.setRange("fit_range", LeftLowLim_incl[i], RightUpLim_incl[i]);

        // Fit full model
        RooFitResult *full_fit_result = model.fitTo(data, RooFit::Range("fit_range"), RooFit::Save(), RooFit::PrintLevel(-1));
        int fitstatus = full_fit_result->status();

        // Plot combined fit
        PlotDataFit(mass_data, data, model, background, crystal, LeftLowLim_incl[i], RightUpLim_incl[i], i, fitstatus);

        // Store fitted values for scale and smearing calculations
        mu_data_fitted[i] = mu_cb.getVal();
        sigma_data_fitted[i] = sigma_cb.getVal();
        inc_mu_data_fitted[i] = mu_cb.getError();
        inc_sigma_data_fitted[i] = sigma_cb.getError();

        std::cout << "Data Bin " << i+1 << ": mu = " << mu_cb.getVal() << " +/- " << mu_cb.getError() 
                  << ", sigma = " << sigma_cb.getVal() << " +/- " << sigma_cb.getError() << std::endl;
    }

    //===========================================
    // STEP 3: CALCULATE AND PLOT SCALE AND SMEARING VS PT
    //===========================================
    std::cout << "\n=== STEP 3: PLOTTING SCALE AND SMEARING VS PT ===" << std::endl;
    
    // Fill the graphs with scale and smearing values
    for (int i = 0; i < Nbins; i++) {
        // Calculate scale: 1 - mu_data/mu_mc
        double scale_value = 1.0 - mu_data_fitted[i]/mu_mc[i];
        double scale_error = (mu_data_fitted[i]/mu_mc[i]) * sqrt(
            (inc_mu_data_fitted[i]/mu_data_fitted[i])*(inc_mu_data_fitted[i]/mu_data_fitted[i]) + 
            (inc_mu_mc[i]/mu_mc[i])*(inc_mu_mc[i]/mu_mc[i])
        );
        
        // Calculate smearing: sigma_data/sigma_mc
        double smearing_value = sigma_data_fitted[i]/val_sigma_mc[i];
        double smearing_error = (sigma_data_fitted[i]/val_sigma_mc[i]) * sqrt(
            (inc_sigma_data_fitted[i]/sigma_data_fitted[i])*(inc_sigma_data_fitted[i]/sigma_data_fitted[i]) + 
            (inc_sigma_mc[i]/val_sigma_mc[i])*(inc_sigma_mc[i]/val_sigma_mc[i])
        );
        
        graph_scale->SetPoint(i, bincenters[i], scale_value);
        graph_scale->SetPointError(i, binhalfwidths[i], scale_error);
        
        graph_smearing->SetPoint(i, bincenters[i], smearing_value);
        graph_smearing->SetPointError(i, binhalfwidths[i], smearing_error);
        
        std::cout << "Bin " << i+1 << ": Scale = " << scale_value << " ± " << scale_error 
                  << ", Smearing = " << smearing_value << " ± " << smearing_error << std::endl;
    }
    
    // Set graph properties
    graph_scale->SetTitle("");
    graph_smearing->SetTitle("");
    
    // Create and plot scale and smearing vs pT
    TCanvas *c_scale_smear = new TCanvas("c_scale_smear", "Raw SC Scale and Smearing vs p_{T}", 1600, 600);
    c_scale_smear->Divide(2,1);
    
    // Scale plot
    c_scale_smear->cd(1);
    gPad->SetLeftMargin(0.15);
    graph_scale->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_scale->GetYaxis()->SetTitle("1 - m^{(J/#psi)}_{data}/m^{(J/#psi)}_{MC}");
    graph_scale->SetMarkerStyle(21);
    graph_scale->SetLineColor(rossoCMS);
    graph_scale->SetMarkerColor(rossoCMS);
    graph_scale->Draw("APE");
    
    // Smearing plot
    c_scale_smear->cd(2);
    gPad->SetLeftMargin(0.15);
    graph_smearing->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_smearing->GetYaxis()->SetTitle("#sigma_{data}/#sigma_{MC}");
    graph_smearing->SetMarkerStyle(21);
    graph_smearing->SetLineColor(violaCMS);
    graph_smearing->SetMarkerColor(violaCMS);
    graph_smearing->Draw("APE");
    
    // Add CMS labels to both pads
    c_scale_smear->cd(1);
    WritePreliminary();
    
    c_scale_smear->cd(2);
    WritePreliminary();
    
    // Save the canvas
    gSystem->Exec("mkdir -p PlotConID2022/RawScaleandSmearing");
    c_scale_smear->SaveAs("PlotConID2022/RawScaleandSmearing/ScaleAndSmearing_RawSC_vsPt.png");
    delete c_scale_smear;

    std::cout << "\n=== FIT COMPLETED FOR ALL BINS ===" << std::endl;
    std::cout << "Results saved in PlotConID2022/RawScaleandSmearing/Pt_bin<N>/" << std::endl;
    std::cout << "Scale and Smearing plot saved as PlotConID2022/RawScaleandSmearing/ScaleAndSmearing_RawSC_vsPt.png" << std::endl;

    fileMC->Close();
    fileData->Close();
}
