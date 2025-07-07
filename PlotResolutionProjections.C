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

#define NbinsPt 6

// Function to fit histogram with Crystal Ball and return the fit result
RooFitResult* FitWithCB(TH1D* hist, RooRealVar& mass, 
                     double mean_init, double sigma_init, 
                     double alphaL_init, double nL_init, 
                     double alphaR_init, double nR_init,
                     double mean_min, double mean_max,
                     int color, RooPlot* frame, const char* name) {
    
    // Convert histogram to RooDataHist
    RooDataHist data(Form("data_%s", name), "Dataset from histogram", mass, hist);
    
    // Define Crystal Ball parameters
    RooRealVar cb_mean(Form("cb_mean_%s", name), "Mean of CB", mean_init, mean_min, mean_max);
    RooRealVar cb_sigma(Form("cb_sigma_%s", name), "Sigma of CB", sigma_init, 0.001, 0.2);
    RooRealVar cb_alphaL(Form("cb_alphaL_%s", name), "AlphaL of CB", alphaL_init, 0, 150);
    RooRealVar cb_nL(Form("cb_nL_%s", name), "nL of CB", nL_init, 0, 200.0);
    RooRealVar cb_alphaR(Form("cb_alphaR_%s", name), "AlphaR of CB", alphaR_init, 0, 150);
    RooRealVar cb_nR(Form("cb_nR_%s", name), "nR of CB", nR_init, 0, 200.0);
    
    // Create Crystal Ball
    RooCrystalBall cb(Form("cb_%s", name), "Asymmetric Crystal Ball", mass, cb_mean, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR);
    
    // First fit with Simplex for stability
    RooFitResult *fit_result_simplex = cb.fitTo(data, 
                                              RooFit::Minimizer("Minuit2", "Simplex"), 
                                              RooFit::MaxCalls(1000000), 
                                              RooFit::Save(), 
                                              RooFit::SumW2Error(true), 
                                              RooFit::PrintLevel(-1));
    
    // Then fit with Migrad for precision
    RooFitResult *fit_result = cb.fitTo(data, 
                                      RooFit::Minimizer("Minuit2", "Migrad"), 
                                      RooFit::MaxCalls(10000000), 
                                      RooFit::Save(), 
                                      RooFit::SumW2Error(true), 
                                      RooFit::PrintLevel(-1));
    
    // Plot data with error bars but no markers
    data.plotOn(frame, 
               RooFit::Name(Form("data_%s", name)), 
               RooFit::MarkerSize(0),
               RooFit::LineColor(color),
               RooFit::DrawOption("E1"));
    
    // Plot the fit
    cb.plotOn(frame, 
             RooFit::Name(Form("fit_%s", name)), 
             RooFit::LineColor(color), 
             RooFit::LineWidth(2));
    
    return fit_result;
}

// Function to add CMS and Simulation labels
void AddCMSLabels() {
    TLatex latex;
    latex.SetTextSize(0.05);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextAngle(0);
    //latex.DrawLatexNDC(0.12, 0.93, "#it{Simulation}");
    
    TLatex cmsLabel;
    cmsLabel.SetNDC();
    cmsLabel.SetTextFont(62);
    cmsLabel.SetTextSize(0.05);
    cmsLabel.SetTextAlign(11);
    cmsLabel.DrawLatex(0.15, 0.92, "CMS");
    
    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextFont(52);
    prelimLabel.SetTextSize(0.05);
    prelimLabel.SetTextAlign(11);
    prelimLabel.DrawLatex(0.24, 0.92, "Preliminary");
}

void FitResolutionHistograms(const char* filename = "outputHistograms_MC.root") {
    // Define colors
    int redColor = kRed;
    int blueColor = kBlue;
    int greenColor = kGreen+2;
    
    // Open file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    // Get histograms
    TH2D *h_rawResolution_vs_pt = (TH2D*)file->Get("h_rawResolution_vs_pt");
    TH2D *h_correctedSCResolution_vs_pt = (TH2D*)file->Get("h_correctedSCResolution_vs_pt");
    TH2D *h_ecaltrkregResolution_vs_pt = (TH2D*)file->Get("h_ecaltrkregResolution_vs_pt");
    
    if (!h_rawResolution_vs_pt || !h_correctedSCResolution_vs_pt || !h_ecaltrkregResolution_vs_pt) {
        std::cerr << "Error: One or more histograms not found" << std::endl;
        file->Close();
        return;
    }

    // Personalized fitting ranges
    double raw_lowlim[NbinsPt] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double raw_uplim[NbinsPt] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    double corrected_lowlim[NbinsPt] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double corrected_uplim[NbinsPt] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};   
    double ecaltrkreg_lowlim[NbinsPt] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double ecaltrkreg_uplim[NbinsPt] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    
    // Parameters for each histogram type and pt bin
    // These are initial values that can be modified for each bin
    double means[NbinsPt] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double sigmas[NbinsPt] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double alphaLs[NbinsPt] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    double nLs[NbinsPt] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    double alphaRs[NbinsPt] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
    double nRs[NbinsPt] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
    double meanMin[NbinsPt] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
    double meanMax[NbinsPt] = {1.4, 1.4, 1.4, 1.2, 1.2, 1.2};
    
    // Arrays to store values for TGraphErrors
    double ptCenters[NbinsPt];
    double ptErrors[NbinsPt];
    double rawMeans[NbinsPt];
    double rawMeanErrors[NbinsPt];
    double corrMeans[NbinsPt];
    double corrMeanErrors[NbinsPt];
    double ecalMeans[NbinsPt];
    double ecalMeanErrors[NbinsPt];
    double rawSigmas[NbinsPt];
    double rawSigmaErrors[NbinsPt];
    double corrSigmas[NbinsPt];
    double corrSigmaErrors[NbinsPt];
    double ecalSigmas[NbinsPt];
    double ecalSigmaErrors[NbinsPt];
    
    // Create canvas for all pt bins
    TCanvas *canvas = new TCanvas("canvas", "Resolution Comparison", 1800, 900);
    canvas->Divide(3, 2); // 3x2 grid for 6 pt bins
    
    // Loop over pt bins
    for (int i = 0; i < NbinsPt; i++) {
        // Get pt bin edges for labels
        double ptLow = h_rawResolution_vs_pt->GetXaxis()->GetBinLowEdge(i+1);
        double ptHigh = h_rawResolution_vs_pt->GetXaxis()->GetBinUpEdge(i+1);
        
        // Calculate bin center and half-width for TGraphErrors
        ptCenters[i] = (ptLow + ptHigh) / 2.0;
        ptErrors[i] = (ptHigh - ptLow) / 2.0;
        
        // Create projections for this pt bin
        TH1D *h_raw_proj = h_rawResolution_vs_pt->ProjectionY(Form("h_raw_proj_bin%d", i+1), i+1, i+1);
        TH1D *h_corr_proj = h_correctedSCResolution_vs_pt->ProjectionY(Form("h_corr_proj_bin%d", i+1), i+1, i+1);
        TH1D *h_ecal_proj = h_ecaltrkregResolution_vs_pt->ProjectionY(Form("h_ecal_proj_bin%d", i+1), i+1, i+1);
        
        // Create RooRealVar for fitting
        RooRealVar massraw("massraw", "E_{gen}/E_{reco}", raw_lowlim[i], raw_uplim[i]);
        RooRealVar masscorr("masscorr", "E_{gen}/E_{reco}", corrected_lowlim[i], corrected_uplim[i]);
        RooRealVar massecaltrk("massecaltrk", "E_{gen}/E_{reco}", ecaltrkreg_lowlim[i], ecaltrkreg_uplim[i]);
        
        // Create plot frame
        RooPlot *frameraw = massraw.frame();
        frameraw->SetTitle("");
        RooPlot *framecorr = masscorr.frame();
        framecorr->SetTitle("");
        RooPlot *frameecaltrk = massecaltrk.frame();
        frameecaltrk->SetTitle("");
        
        // Fit histograms and add to frame
        RooFitResult *raw_fit = FitWithCB(h_raw_proj, massraw, means[i], sigmas[i], alphaLs[i], nLs[i], alphaRs[i], nRs[i], meanMin[i], meanMax[i], redColor, frameraw, Form("raw_bin%d", i+1));
        RooFitResult *corr_fit = FitWithCB(h_corr_proj, masscorr, means[i], sigmas[i], alphaLs[i], nLs[i], alphaRs[i], nRs[i], meanMin[i], meanMax[i], blueColor, framecorr, Form("corr_bin%d", i+1));
        RooFitResult *ecal_fit = FitWithCB(h_ecal_proj, massecaltrk, means[i], sigmas[i], alphaLs[i], nLs[i], alphaRs[i], nRs[i], meanMin[i], meanMax[i], greenColor, frameecaltrk, Form("ecal_bin%d", i+1));
        
        // Set titles and range
        frameraw->GetXaxis()->SetTitle("E_{gen}/E_{reco}");
        frameraw->GetYaxis()->SetTitle("Entries");
        frameraw->GetXaxis()->SetRangeUser(0.5, 3.0);
        framecorr->GetXaxis()->SetTitle("E_{gen}/E_{reco}");
        framecorr->GetYaxis()->SetTitle("Entries");
        framecorr->GetXaxis()->SetRangeUser(0.5, 3.0);
        frameecaltrk->GetXaxis()->SetTitle("E_{gen}/E_{reco}");
        frameecaltrk->GetYaxis()->SetTitle("Entries");
        frameecaltrk->GetXaxis()->SetRangeUser(0.5, 3.0);
        
        // Draw frame
        frameraw->Draw();
        framecorr->Draw("SAME");
        frameecaltrk->Draw("SAME");
        
        // Add pt range label
        TPaveText *ptLabel = new TPaveText(0.6, 0.7, 0.85, 0.8, "NDC");
        ptLabel->SetBorderSize(0);
        ptLabel->SetFillStyle(0);
        ptLabel->AddText(Form("%.1f < p_{T} < %.1f GeV", ptLow, ptHigh));
        ptLabel->Draw();
        
        // Add legend
        TLegend *legend = new TLegend(0.6, 0.55, 0.85, 0.7);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(frameraw->findObject(Form("data_raw_bin%d", i+1)), "Raw SC", "l");
        legend->AddEntry(framecorr->findObject(Form("data_corr_bin%d", i+1)), "Corrected SC", "l");
        legend->AddEntry(frameecaltrk->findObject(Form("data_ecal_bin%d", i+1)), "ECAL+Track", "l");
        legend->Draw();
        
        // Add CMS labels
        AddCMSLabels();
        
        // Create individual canvas for this pt bin for higher resolution
        TCanvas *individualCanvas = new TCanvas(Form("canvas_bin%d", i+1), 
                                            Form("Resolution Comparison %.1f < p_{T} < %.1f GeV", ptLow, ptHigh), 
                                            800, 600);
        individualCanvas->SetLeftMargin(0.15);
        
        // Clone frame for individual canvas
        RooPlot *frameClone = (RooPlot*)frameraw->Clone(Form("frameraw_bin%d", i+1));
        RooPlot *frameCloneCorr = (RooPlot*)framecorr->Clone(Form("framecorr_bin%d", i+1));
        RooPlot *frameCloneEcal = (RooPlot*)frameecaltrk->Clone(Form("frameecal_bin%d", i+1));
        frameCloneEcal->Draw();
        frameClone->Draw("SAME");
        frameCloneCorr->Draw("SAME");
        ptLabel->Draw();
        legend->Draw();
        AddCMSLabels();
        
        // Save individual canvas
        individualCanvas->SaveAs(Form("PlotConID2022/ResolutionfromEGMmcprompt/Resolution_PtBin%d_%.1f_%.1f_GeV.png", i+1, ptLow, ptHigh));
        delete individualCanvas;
        
        // NEW: Create canvas with 3 separate subcanvases for this pt bin
        TCanvas *separateCanvas = new TCanvas(Form("separate_canvas_bin%d", i+1),
                                           Form("Resolution Components %.1f < p_{T} < %.1f GeV", ptLow, ptHigh),
                                           3000, 400);
        separateCanvas->Divide(3, 1); // Divide into 3 pads horizontally
        
        // First pad: Raw resolution
        separateCanvas->cd(1);
        gPad->SetLeftMargin(0.15);
        frameraw->SetTitle("Raw SC");
        frameraw->Draw();
        TPaveText *ptLabel1 = (TPaveText*)ptLabel->Clone("ptLabel1");
        ptLabel1->Draw();
        // CMS labels removed as requested
        
        // Add textbox with fit parameters for Raw
        TPaveText *paramBox1 = new TPaveText(0.65, 0.35, 0.89, 0.70, "NDC");
        paramBox1->SetBorderSize(1);
        paramBox1->SetFillColor(0);
        paramBox1->SetTextAlign(12);
        paramBox1->SetTextSize(0.035);
        
        // Extract fit parameters and errors for Raw
        RooRealVar *cb_mean_raw = (RooRealVar*)raw_fit->floatParsFinal().find(Form("cb_mean_raw_bin%d", i+1));
        RooRealVar *cb_sigma_raw = (RooRealVar*)raw_fit->floatParsFinal().find(Form("cb_sigma_raw_bin%d", i+1));
        RooRealVar *cb_alphaL_raw = (RooRealVar*)raw_fit->floatParsFinal().find(Form("cb_alphaL_raw_bin%d", i+1));
        RooRealVar *cb_alphaR_raw = (RooRealVar*)raw_fit->floatParsFinal().find(Form("cb_alphaR_raw_bin%d", i+1));
        RooRealVar *cb_nL_raw = (RooRealVar*)raw_fit->floatParsFinal().find(Form("cb_nL_raw_bin%d", i+1));
        RooRealVar *cb_nR_raw = (RooRealVar*)raw_fit->floatParsFinal().find(Form("cb_nR_raw_bin%d", i+1));
        
        paramBox1->AddText(Form("Mean: %.3f #pm %.3f", cb_mean_raw->getVal(), cb_mean_raw->getError()));
        paramBox1->AddText(Form("#sigma: %.4f #pm %.4f", cb_sigma_raw->getVal(), cb_sigma_raw->getError()));
        paramBox1->AddText(Form("#alpha_{L}: %.3f #pm %.3f", cb_alphaL_raw->getVal(), cb_alphaL_raw->getError()));
        paramBox1->AddText(Form("n_{L}: %.2f #pm %.2f", cb_nL_raw->getVal(), cb_nL_raw->getError()));
        paramBox1->AddText(Form("#alpha_{R}: %.3f #pm %.3f", cb_alphaR_raw->getVal(), cb_alphaR_raw->getError()));
        paramBox1->AddText(Form("n_{R}: %.2f #pm %.2f", cb_nR_raw->getVal(), cb_nR_raw->getError()));
        paramBox1->Draw();
        
        // Second pad: Corrected resolution
        separateCanvas->cd(2);
        gPad->SetLeftMargin(0.15);
        framecorr->SetTitle("Corrected SC");
        framecorr->Draw();
        TPaveText *ptLabel2 = (TPaveText*)ptLabel->Clone("ptLabel2");
        ptLabel2->Draw();
        // CMS labels removed as requested
        
        // Add textbox with fit parameters for Corrected
        TPaveText *paramBox2 = new TPaveText(0.65, 0.35, 0.89, 0.70, "NDC");
        paramBox2->SetBorderSize(1);
        paramBox2->SetFillColor(0);
        paramBox2->SetTextAlign(12);
        paramBox2->SetTextSize(0.035);
        
        // Extract fit parameters and errors for Corrected
        RooRealVar *cb_mean_corr = (RooRealVar*)corr_fit->floatParsFinal().find(Form("cb_mean_corr_bin%d", i+1));
        RooRealVar *cb_sigma_corr = (RooRealVar*)corr_fit->floatParsFinal().find(Form("cb_sigma_corr_bin%d", i+1));
        RooRealVar *cb_alphaL_corr = (RooRealVar*)corr_fit->floatParsFinal().find(Form("cb_alphaL_corr_bin%d", i+1));
        RooRealVar *cb_alphaR_corr = (RooRealVar*)corr_fit->floatParsFinal().find(Form("cb_alphaR_corr_bin%d", i+1));
        RooRealVar *cb_nL_corr = (RooRealVar*)corr_fit->floatParsFinal().find(Form("cb_nL_corr_bin%d", i+1));
        RooRealVar *cb_nR_corr = (RooRealVar*)corr_fit->floatParsFinal().find(Form("cb_nR_corr_bin%d", i+1));
        
        paramBox2->AddText(Form("Mean: %.3f #pm %.3f", cb_mean_corr->getVal(), cb_mean_corr->getError()));
        paramBox2->AddText(Form("#sigma: %.4f #pm %.4f", cb_sigma_corr->getVal(), cb_sigma_corr->getError()));
        paramBox2->AddText(Form("#alpha_{L}: %.3f #pm %.3f", cb_alphaL_corr->getVal(), cb_alphaL_corr->getError()));
        paramBox2->AddText(Form("n_{L}: %.2f #pm %.2f", cb_nL_corr->getVal(), cb_nL_corr->getError()));
        paramBox2->AddText(Form("#alpha_{R}: %.3f #pm %.3f", cb_alphaR_corr->getVal(), cb_alphaR_corr->getError()));
        paramBox2->AddText(Form("n_{R}: %.2f #pm %.2f", cb_nR_corr->getVal(), cb_nR_corr->getError()));
        paramBox2->Draw();
        
        // Third pad: ECAL+Track resolution
        separateCanvas->cd(3);
        gPad->SetLeftMargin(0.15);
        frameecaltrk->SetTitle("ECAL+Track");
        frameecaltrk->Draw();
        TPaveText *ptLabel3 = (TPaveText*)ptLabel->Clone("ptLabel3");
        ptLabel3->Draw();
        // CMS labels removed as requested
        
        // Add textbox with fit parameters for ECAL+Track
        TPaveText *paramBox3 = new TPaveText(0.65, 0.35, 0.89, 0.70, "NDC");
        paramBox3->SetBorderSize(1);
        paramBox3->SetFillColor(0);
        paramBox3->SetTextAlign(12);
        paramBox3->SetTextSize(0.035);
        
        // Extract fit parameters and errors for ECAL+Track
        RooRealVar *cb_mean_ecal = (RooRealVar*)ecal_fit->floatParsFinal().find(Form("cb_mean_ecal_bin%d", i+1));
        RooRealVar *cb_sigma_ecal = (RooRealVar*)ecal_fit->floatParsFinal().find(Form("cb_sigma_ecal_bin%d", i+1));
        RooRealVar *cb_alphaL_ecal = (RooRealVar*)ecal_fit->floatParsFinal().find(Form("cb_alphaL_ecal_bin%d", i+1));
        RooRealVar *cb_alphaR_ecal = (RooRealVar*)ecal_fit->floatParsFinal().find(Form("cb_alphaR_ecal_bin%d", i+1));
        RooRealVar *cb_nL_ecal = (RooRealVar*)ecal_fit->floatParsFinal().find(Form("cb_nL_ecal_bin%d", i+1));
        RooRealVar *cb_nR_ecal = (RooRealVar*)ecal_fit->floatParsFinal().find(Form("cb_nR_ecal_bin%d", i+1));
        
        paramBox3->AddText(Form("Mean: %.3f #pm %.3f", cb_mean_ecal->getVal(), cb_mean_ecal->getError()));
        paramBox3->AddText(Form("#sigma: %.4f #pm %.4f", cb_sigma_ecal->getVal(), cb_sigma_ecal->getError()));
        paramBox3->AddText(Form("#alpha_{L}: %.3f #pm %.3f", cb_alphaL_ecal->getVal(), cb_alphaL_ecal->getError()));
        paramBox3->AddText(Form("n_{L}: %.2f #pm %.2f", cb_nL_ecal->getVal(), cb_nL_ecal->getError()));
        paramBox3->AddText(Form("#alpha_{R}: %.3f #pm %.3f", cb_alphaR_ecal->getVal(), cb_alphaR_ecal->getError()));
        paramBox3->AddText(Form("n_{R}: %.2f #pm %.2f", cb_nR_ecal->getVal(), cb_nR_ecal->getError()));
        paramBox3->Draw();
        
        // Save the separate canvas
        separateCanvas->SaveAs(Form("PlotConID2022/ResolutionfromEGMmcprompt/ResolutionComponents_PtBin%d_%.1f_%.1f_GeV.png", i+1, ptLow, ptHigh));
        //separateCanvas->SaveAs(Form("PlotConID2022/ResolutionfromEGMmcprompt/ResolutionComponents_PtBin%d_%.1f_%.1f_GeV.pdf", i+1, ptLow, ptHigh));
        delete separateCanvas;

        // Store means and sigmas with their errors
        rawMeans[i] = cb_mean_raw->getVal();
        rawMeanErrors[i] = cb_mean_raw->getError();
        corrMeans[i] = cb_mean_corr->getVal();
        corrMeanErrors[i] = cb_mean_corr->getError();
        ecalMeans[i] = cb_mean_ecal->getVal();
        ecalMeanErrors[i] = cb_mean_ecal->getError();
        
        rawSigmas[i] = cb_sigma_raw->getVal();
        rawSigmaErrors[i] = cb_sigma_raw->getError();
        corrSigmas[i] = cb_sigma_corr->getVal();
        corrSigmaErrors[i] = cb_sigma_corr->getError();
        ecalSigmas[i] = cb_sigma_ecal->getVal();
        ecalSigmaErrors[i] = cb_sigma_ecal->getError();
    }
    
    // Save combined canvas
    canvas->SaveAs("PlotConID2022/ResolutionfromEGMmcprompt/Resolution_AllPtBins.png");
    
    // Create TGraphErrors for means vs pt
    TGraphErrors *gr_raw_mean = new TGraphErrors(NbinsPt, ptCenters, rawMeans, ptErrors, rawMeanErrors);
    gr_raw_mean->SetName("gr_raw_mean");
    gr_raw_mean->SetTitle("Raw SC");
    gr_raw_mean->SetMarkerStyle(20);
    gr_raw_mean->SetMarkerColor(redColor);
    gr_raw_mean->SetLineColor(redColor);
    gr_raw_mean->SetMarkerSize(1.2);
    
    TGraphErrors *gr_corr_mean = new TGraphErrors(NbinsPt, ptCenters, corrMeans, ptErrors, corrMeanErrors);
    gr_corr_mean->SetName("gr_corr_mean");
    gr_corr_mean->SetTitle("Corrected SC");
    gr_corr_mean->SetMarkerStyle(21);
    gr_corr_mean->SetMarkerColor(blueColor);
    gr_corr_mean->SetLineColor(blueColor);
    gr_corr_mean->SetMarkerSize(1.2);
    
    TGraphErrors *gr_ecal_mean = new TGraphErrors(NbinsPt, ptCenters, ecalMeans, ptErrors, ecalMeanErrors);
    gr_ecal_mean->SetName("gr_ecal_mean");
    gr_ecal_mean->SetTitle("ECAL+Track");
    gr_ecal_mean->SetMarkerStyle(22);
    gr_ecal_mean->SetMarkerColor(greenColor);
    gr_ecal_mean->SetLineColor(greenColor);
    gr_ecal_mean->SetMarkerSize(1.2);
    
    // Create TGraphErrors for sigmas vs pt
    TGraphErrors *gr_raw_sigma = new TGraphErrors(NbinsPt, ptCenters, rawSigmas, ptErrors, rawSigmaErrors);
    gr_raw_sigma->SetName("gr_raw_sigma");
    gr_raw_sigma->SetTitle("Raw SC");
    gr_raw_sigma->SetMarkerStyle(20);
    gr_raw_sigma->SetMarkerColor(redColor);
    gr_raw_sigma->SetLineColor(redColor);
    gr_raw_sigma->SetMarkerSize(1.2);
    
    TGraphErrors *gr_corr_sigma = new TGraphErrors(NbinsPt, ptCenters, corrSigmas, ptErrors, corrSigmaErrors);
    gr_corr_sigma->SetName("gr_corr_sigma");
    gr_corr_sigma->SetTitle("Corrected SC");
    gr_corr_sigma->SetMarkerStyle(21);
    gr_corr_sigma->SetMarkerColor(blueColor);
    gr_corr_sigma->SetLineColor(blueColor);
    gr_corr_sigma->SetMarkerSize(1.2);
    
    TGraphErrors *gr_ecal_sigma = new TGraphErrors(NbinsPt, ptCenters, ecalSigmas, ptErrors, ecalSigmaErrors);
    gr_ecal_sigma->SetName("gr_ecal_sigma");
    gr_ecal_sigma->SetTitle("ECAL+Track");
    gr_ecal_sigma->SetMarkerStyle(22);
    gr_ecal_sigma->SetMarkerColor(greenColor);
    gr_ecal_sigma->SetLineColor(greenColor);
    gr_ecal_sigma->SetMarkerSize(1.2);
    
    // Create canvas for Mean vs pT
    TCanvas *canvasMean = new TCanvas("canvasMean", "Mean vs p_{T}", 800, 600);
    canvasMean->SetLeftMargin(0.15);
    canvasMean->SetRightMargin(0.05);
    
    // Create a multi-graph to hold all three graphs
    TMultiGraph *mg_mean = new TMultiGraph();
    mg_mean->SetTitle("Mean vs p_{T}");
    mg_mean->Add(gr_raw_mean, "PL");
    mg_mean->Add(gr_corr_mean, "PL");
    mg_mean->Add(gr_ecal_mean, "PL");
    mg_mean->Draw("A");
    
    // Set axis titles and ranges
    mg_mean->GetXaxis()->SetTitle("p_{T} [GeV]");
    mg_mean->GetYaxis()->SetTitle("Mean of E_{gen}/E_{reco}");
    mg_mean->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust as needed
    
    // Add legend
    TLegend *legendMean = new TLegend(0.6, 0.7, 0.89, 0.89);
    legendMean->SetBorderSize(0);
    legendMean->SetFillStyle(0);
    legendMean->AddEntry(gr_raw_mean, "Raw SC", "P");
    legendMean->AddEntry(gr_corr_mean, "Corrected SC", "P");
    legendMean->AddEntry(gr_ecal_mean, "ECAL+Track", "P");
    legendMean->Draw();
    
    // Add reference line at 1.0
    TLine *lineMean = new TLine(mg_mean->GetXaxis()->GetXmin(), 1.0, 
                               mg_mean->GetXaxis()->GetXmax(), 1.0);
    lineMean->SetLineColor(kGray+2);
    lineMean->SetLineStyle(2);
    lineMean->Draw();
    
    // Add CMS labels
    AddCMSLabels();
    
    // Create canvas for Sigma vs pT
    TCanvas *canvasSigma = new TCanvas("canvasSigma", "Sigma vs p_{T}", 800, 600);
    canvasSigma->SetLeftMargin(0.15);
    canvasSigma->SetRightMargin(0.05);
    
    // Create a multi-graph to hold all three graphs
    TMultiGraph *mg_sigma = new TMultiGraph();
    mg_sigma->SetTitle("");
    mg_sigma->Add(gr_raw_sigma, "P");
    mg_sigma->Add(gr_corr_sigma, "P");
    mg_sigma->Add(gr_ecal_sigma, "P");
    mg_sigma->Draw("A");
    
    // Set axis titles and ranges
    mg_sigma->GetXaxis()->SetTitle("p_{T} [GeV]");
    mg_sigma->GetYaxis()->SetTitle("#sigma of E_{gen}/E_{reco}");
    mg_sigma->GetYaxis()->SetRangeUser(0.0, 0.15); // Adjust as needed
    
    // Add legend
    TLegend *legendSigma = new TLegend(0.6, 0.7, 0.89, 0.89);
    legendSigma->SetBorderSize(0);
    legendSigma->SetFillStyle(0);
    legendSigma->AddEntry(gr_raw_sigma, "Raw SC", "P");
    legendSigma->AddEntry(gr_corr_sigma, "Corrected SC", "P");
    legendSigma->AddEntry(gr_ecal_sigma, "ECAL+Track", "P");
    legendSigma->Draw();
    
    // Add CMS labels
    AddCMSLabels();
    
    // Save canvases
    canvasMean->SaveAs("PlotConID2022/ResolutionfromEGMmcprompt/MeanVsPt_Comparison.png");
    canvasSigma->SaveAs("PlotConID2022/ResolutionfromEGMmcprompt/SigmaVsPt_Comparison.png");
    
    // Save graphs to a root file for future use
    TFile *outFile = new TFile("PlotConID2022/ResolutionfromEGMmcprompt/ResolutionGraphs.root", "RECREATE");
    gr_raw_mean->Write();
    gr_corr_mean->Write();
    gr_ecal_mean->Write();
    gr_raw_sigma->Write();
    gr_corr_sigma->Write();
    gr_ecal_sigma->Write();
    outFile->Close();
    
    // Clean up
    delete canvasMean;
    delete canvasSigma;
    delete outFile;
    file->Close();
}