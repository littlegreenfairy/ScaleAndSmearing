#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooCrystalBall.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <iostream>
#include <TLatex.h>
#include <TLegend.h>
#include <TColor.h>
#include <TSystem.h>
#include <TPaveText.h>

#define nbinsPt 6
#define nBinsRun 9
#define Nsigma 1.0

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

void PlotBackgroundFit(RooRealVar& mass, RooDataHist& data, RooAddPdf& background, int i, int j) {
    // Create a canvas to draw the background fit
    TCanvas* cBackgroundFit = new TCanvas(Form("cBackgroundFit_%d_%d", i+1, j+1), "Background Fit", 800, 600);
    gPad->SetLeftMargin(0.13);

    // Create a frame for the mass variable to hold the plot
    RooPlot* frame = mass.frame();
    
    // Plot the data points on the frame
    data.plotOn(frame, RooFit::Name("Data"), RooFit::MarkerStyle(20), RooFit::MarkerSize(0.7));

    // Plot the background fit function
    background.plotOn(frame, RooFit::Name("BackgroundFit"), RooFit::LineColor(kMagenta), RooFit::LineStyle(kDashed));

    // Add labels and titles
    frame->SetTitle(Form("Background Fit for Full Electron Pt Bin %d and Run Bin %d", i+1, j+1));
    frame->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle("Events");

    // Draw the frame on the canvas
    frame->Draw();

    // Add legend for clarity
    TLegend* legend = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend->AddEntry(frame->findObject("Data"), "Data", "p");
    legend->AddEntry(frame->findObject("BackgroundFit"), "Background Fit", "l");
    legend->Draw();

    WritePreliminary();

    // Save the canvas as a file
    cBackgroundFit->SaveAs(Form("PlotConID2022/ScaleonFullEle/Pt_bin%d/BackgroundFit_FullEle_Pt%d_Run%d.png", i+1, i+1, j+1));
    delete cBackgroundFit;
}

void PlotDataFit(RooRealVar& mass, RooDataHist& data, RooAddPdf& model, RooAddPdf& background, RooCrystalBall& crystal, double leftlim, double rightlim, int i, int j, int fitstatus) {

    // Create canvas for plotting
    TCanvas *c = new TCanvas(Form("c_bin_%d_%d", i+1, j+1), Form("Bin %d_%d", i+1, j+1), 950, 700);
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
    c->SaveAs(Form("PlotConID2022/ScaleonFullEle/Pt_bin%d/DataFit_FullEle_Pt%d_Run%d.png", i+1, i+1, j+1));
    delete c;
}

void ScaleonFullElectron() {
    // Create directories for output plots
    for(int i = 0; i < nbinsPt; i++) {
        gSystem->Exec(Form("mkdir -p PlotConID2022/ScaleonFullEle/Pt_bin%d", i+1));
    }

    // Open the data and MC files
    TFile *dataFile = TFile::Open("outputHistograms_DATA_partF.root", "READ");
    TFile *mcFile = TFile::Open("outputHistograms_MC.root", "READ");
    
    if (!dataFile || dataFile->IsZombie()) {
        std::cerr << "Error: Cannot open outputHistograms_DATA_partF.root" << std::endl;
        return;
    }
    
    if (!mcFile || mcFile->IsZombie()) {
        std::cerr << "Error: Cannot open outputHistograms_MC.root" << std::endl;
        return;
    }
    
    // Get the histograms
    TH3D* h_fullEle_invMass_pt_runN = (TH3D*)dataFile->Get("h_fullEle_invMass_pt_runN");
    TH2D* h_invmass_vs_pt_fullele = (TH2D*)mcFile->Get("h_invmass_vs_pt_fullele");
    
    if (!h_fullEle_invMass_pt_runN) {
        std::cerr << "Error: h_fullEle_invMass_pt_runN not found in data file!" << std::endl;
        dataFile->Close();
        mcFile->Close();
        return;
    }
    
    if (!h_invmass_vs_pt_fullele) {
        std::cerr << "Error: h_invmass_vs_pt_fullele not found in MC file!" << std::endl;
        dataFile->Close();
        mcFile->Close();
        return;
    }

    // Arrays for MC fit parameters (customizable)
    double means[nbinsPt], sigmas[nbinsPt], alphaLs[nbinsPt], nLs[nbinsPt], alphaRs[nbinsPt], nRs[nbinsPt];
    double meanInf[nbinsPt], meanSup[nbinsPt];
    double mu_mc[nbinsPt], inc_mu_mc[nbinsPt], val_sigma_mc[nbinsPt], inc_sigma_mc[nbinsPt];
    double val_nL[nbinsPt], inc_nL[nbinsPt], val_alphL[nbinsPt], inc_alphL[nbinsPt];
    double val_nR[nbinsPt], inc_nR[nbinsPt], val_alphR[nbinsPt], inc_alphR[nbinsPt];

    // Initialize MC fit parameters (customizable)
    for (int i = 0; i < nbinsPt; i++) {
        means[i] = 3.0;
        sigmas[i] = 0.14;
        alphaLs[i] = 0.92;
        nLs[i] = 16.2;
        alphaRs[i] = 1.67;
        nRs[i] = 4.59;
        meanInf[i] = means[i] - 2 * sigmas[i];
        meanSup[i] = means[i] + 2 * sigmas[i];
    }

    // Mass range and rebin factors for each bin (customizable)
    double massMin[nbinsPt] = {2.2, 2.2, 2.3, 2.1, 2.65, 2.65};
    double massMax[nbinsPt] = {3.9, 3.9, 3.8, 3.8, 3.5, 3.5};
    int rebin_factor[nbinsPt] = {1, 1, 1, 1, 1, 2};

    // Data fit parameters - sideband ranges (customizable for each pt and run bin)
    double LeftLowLim[nbinsPt][nBinsRun] = {
        {1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15},
        {1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3},
        {1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6},
        {1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6},
        {1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65},
        {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}
    };
    
    double LeftUpLim[nbinsPt][nBinsRun] = {
        {2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4},
        {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5},
        {2.45, 2.45, 2.45, 2.45, 2.45, 2.45, 2.45, 2.45, 2.45},
        {2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4},
        {2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65},
        {2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7}
    };
    
    double RightLowLim[nbinsPt][nBinsRun] = {
        {3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 3.65},
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55},
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55},
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55},
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55},
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55}
    };
    
    double RightUpLim[nbinsPt][nBinsRun] = {
        {5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2},
        {5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2},
        {4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8},
        {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0},
        {5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2},
        {5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 5.2}
    };

    // Crystal Ball parameters for data fit (customizable for each pt and run bin)
    double Mucb_ini[nbinsPt][nBinsRun] = {
        {3.0449, 3.0449, 3.0449, 3.0449, 3.0449, 3.0449, 3.0449, 3.0449, 3.0449},
        {3.0421, 3.0421, 3.0421, 3.0421, 3.0421, 3.0421, 3.0421, 3.0421, 3.0421},
        {3.0343, 3.0343, 3.0343, 3.0343, 3.0343, 3.0343, 3.0343, 3.0343, 3.0343},
        {3.0714, 3.0714, 3.0714, 3.0714, 3.0714, 3.0714, 3.0714, 3.0714, 3.0714},
        {3.0839, 3.0839, 3.0839, 3.0839, 3.0839, 3.0839, 3.0839, 3.0839, 3.0839},
        {3.1217, 3.1217, 3.1217, 3.1217, 3.1217, 3.1217, 3.1217, 3.1217, 3.1217}
    };
    
    double Mucb_lowlim[nbinsPt][nBinsRun] = {
        {2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9},
        {2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9},
        {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5},
        {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5},
        {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5},
        {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0}
    };
    
    double Mucb_uplim[nbinsPt][nBinsRun] = {
        {3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2},
        {3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2},
        {3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2},
        {3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2},
        {3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2},
        {3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2}
    };
    
    double Sigmacb_ini[nbinsPt][nBinsRun] = {
        {0.1613, 0.1613, 0.1613, 0.1613, 0.1613, 0.1613, 0.1613, 0.1613, 0.1613},
        {0.149, 0.149, 0.149, 0.149, 0.149, 0.149, 0.149, 0.149, 0.149},
        {0.1440, 0.1440, 0.1440, 0.1440, 0.1440, 0.1440, 0.1440, 0.1440, 0.1440},
        {0.1203, 0.1203, 0.1203, 0.1203, 0.1203, 0.1203, 0.1203, 0.1203, 0.1203},
        {0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11},
        {0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09}
    };
    
    double Sigmacb_uplim[nbinsPt][nBinsRun] = {
        {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
        {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
        {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
        {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
        {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2},
        {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}
    };
    
    double Sigmacb_lowlim[nbinsPt][nBinsRun] = {
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}
    };

    // Gaussian parameters for background fit in data (customizable for each pt and run bin)
    double gauss_mu_init[nbinsPt][nBinsRun] = {
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55},
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55},
        {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5},
        {3.62, 3.62, 3.62, 3.62, 3.62, 3.62, 3.62, 3.62, 3.62},
        {3.6365, 3.6365, 3.6365, 3.6365, 3.6365, 3.6365, 3.6365, 3.6365, 3.6365},
        {3.6195, 3.6195, 3.6195, 3.6195, 3.6195, 3.6195, 3.6195, 3.6195, 3.6195}
    };
    
    double gauss_mu_low[nbinsPt][nBinsRun] = {
        {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5},
        {3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4},
        {3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4},
        {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5},
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55},
        {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5}
    };
    
    double gauss_mu_up[nbinsPt][nBinsRun] = {
        {3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7}
    };
    
    double gauss_sigma_init[nbinsPt][nBinsRun] = {
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15}
    };
    
    double gauss_sigma_low[nbinsPt][nBinsRun] = {
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}
    };
    
    double gauss_sigma_up[nbinsPt][nBinsRun] = {
        {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4},
        {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4},
        {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4},
        {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4},
        {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4},
        {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3}
    };

    // Inclusive parameters for pT-inclusive fits (1D arrays)
    double LeftLowLim_incl[nbinsPt] = {1.15, 1.3, 1.6, 1.6, 1.65, 2};
    double LeftUpLim_incl[nbinsPt] = {2.4, 2.5, 2.45, 2.4, 2.65, 2.7};
    double RightLowLim_incl[nbinsPt] = {3.65, 3.55, 3.55, 3.55, 3.55, 3.55};
    double RightUpLim_incl[nbinsPt] = {5.2, 5.2, 4.8, 5, 5.2, 5.2};
    double Mucb_ini_incl[nbinsPt] = {3.0449, 3.0421, 3.0343, 3.0714, 3.0839, 3.1217};
    double Mucb_lowlim_incl[nbinsPt] = {2.9, 2.9, 2.5, 2.5, 2.5, 3};
    double Mucb_uplim_incl[nbinsPt] = {3.2, 3.2, 3.2, 3.2, 3.2, 3.2};
    double Sigmacb_ini_incl[nbinsPt] = {0.1613, 0.149, 0.1440, 0.1203, 0.11, 0.09};
    double Sigmacb_uplim_incl[nbinsPt] = {0.3, 0.3, 0.3, 0.3, 0.2, 0.2};
    double Sigmacb_lowlim_incl[nbinsPt] = {0.05, 0.05, 0.05, 0.05, 0, 0.05};
    double gauss_mu_init_incl[nbinsPt] = {3.55, 3.55, 3.5, 3.62, 3.6365, 3.6195};
    double gauss_mu_low_incl[nbinsPt] = {3.5, 3.4, 3.4, 3.5, 3.55, 3.5};
    double gauss_mu_up_incl[nbinsPt] = {3.7, 3.7, 3.7, 3.7, 3.7, 3.7};
    double gauss_sigma_init_incl[nbinsPt] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.15};
    double gauss_sigma_low_incl[nbinsPt] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.1};
    double gauss_sigma_up_incl[nbinsPt] = {0.4, 0.4, 0.4, 0.4, 0.4, 0.3};

    //===========================================
    // STEP 1: FIT MC PROJECTIONS WITH DOUBLE CRYSTAL BALL
    //===========================================
    std::cout << "=== STEP 1: FITTING MC PROJECTIONS ===" << std::endl;
    
    for (int i = 0; i < nbinsPt; i++) {
        // Project the Y (mass) for pt bin i using bin i+1
        TH1D* proj_mc = h_invmass_vs_pt_fullele->ProjectionY(Form("fullEle_MC_proj_bin_%d", i+1), i+1, i+1);
        
        if (!proj_mc || proj_mc->GetEntries() == 0) {
            std::cerr << "Empty MC projection for pt bin " << i+1 << std::endl;
            continue;
        }
        proj_mc->Rebin(rebin_factor[i]);

        // Get mean and RMS for initialization
        means[i] = proj_mc->GetMean();
        sigmas[i] = proj_mc->GetRMS() / 2.0;

        // Define mass variable
        RooRealVar mass("mass", "m(e^{+}e^{-})", massMin[i], massMax[i]); 

        // Convert histogram to RooDataHist
        RooDataHist data_mc("data_mc", "MC Dataset", mass, proj_mc);

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
        canvas_mc->SaveAs(Form("PlotConID2022/ScaleonFullEle/Pt_bin%d/MCFit_FullEle_Pt%d.png", i+1, i+1));
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
    }

    //===========================================
    // STEP 2: FIT DATA PROJECTIONS
    //===========================================
    std::cout << "\n=== STEP 2: FITTING DATA PROJECTIONS ===" << std::endl;

    // Create 2D histograms for scale and smearing corrections
    double Ptbins[] = {4, 7, 9, 11, 14, 20, 40}; 
    double runBins[] = {356309, 356900, 357538, 357732, 360000, 360400, 361000, 361600, 362200, 362760};
    TH2D *h_scale_fullEle = new TH2D("h_scale_fullEle", "Scale between data and MC (Full Electron); p_{T} [GeV]; Run number", nbinsPt, Ptbins, nBinsRun, runBins);
    TH2D *h_smearing_fullEle = new TH2D("h_smearing_fullEle", "Smearing corrections (Full Electron); p_{T} [GeV]; Run number", nbinsPt, Ptbins, nBinsRun, runBins);

    for(int i = 0; i < nbinsPt; i++) { // Loop over pT bins
        for(int j = 0; j < nBinsRun; j++) { // Loop over run number bins
            
            // Project onto invariant mass axis (Z) for specific pt and run bins
            TH1D* proj_invmass = h_fullEle_invMass_pt_runN->ProjectionZ(Form("fullEle_proj_bins_%d_%d", i+1, j+1), i+1, i+1, j+1, j+1);
            
            if (!proj_invmass || proj_invmass->GetEntries() == 0) {
                std::cerr << "Empty Data projection for pt bin " << i+1 << ", run bin " << j+1 << std::endl;
                continue;
            }

            // Define mass variable for data fit
            RooRealVar mass_data("mass_data", "m(e^{+}e^{-})", 0, 6); 

            // Convert histogram to RooDataHist
            RooDataHist data("data", "Data Dataset", mass_data, proj_invmass);

            //--- BACKGROUND FIT FIRST ---
            // Define polynomial background parameters
            RooRealVar A(Form("A_%d_%d", i+1, j+1), "4th deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar B(Form("B_%d_%d", i+1, j+1), "3rd deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar C(Form("C_%d_%d", i+1, j+1), "2nd deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar D(Form("D_%d_%d", i+1, j+1), "1st deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar E(Form("E_%d_%d", i+1, j+1), "0 deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
            RooPolynomial poly("poly", "Polynomial of 4th degree", mass_data, RooArgList(A, B, C, D, E));

            // Define Gaussian for background
            RooRealVar gauss_mu(Form("gauss_mu_%d_%d", i+1, j+1), "Gaussian mean", gauss_mu_init[i][j], gauss_mu_low[i][j], gauss_mu_up[i][j]);
            RooRealVar gauss_sigma(Form("gauss_sigma_%d_%d", i+1, j+1), "Gaussian sigma", gauss_sigma_init[i][j], gauss_sigma_low[i][j], gauss_sigma_up[i][j]);
            RooGaussian gauss("gauss", "Gaussian component", mass_data, gauss_mu, gauss_sigma);

            // Background fractions
            RooRealVar frac_gauss(Form("frac_gauss_%d_%d", i+1, j+1), "Fraction of Gaussian", 0.3, 0.0, 1.0);
            RooAddPdf background("background", "Background model", RooArgList(poly, gauss), RooArgList(frac_gauss));

            // Fit background in sidebands
            mass_data.setRange("range1", LeftLowLim[i][j], LeftUpLim[i][j]);
            mass_data.setRange("range2", RightLowLim[i][j], RightUpLim[i][j]);

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
            PlotBackgroundFit(mass_data, data, background, i, j);

            //--- SIGNAL + BACKGROUND FIT ---
            // Initialize Crystal Ball parameters with values from MC fit and 2D arrays
            RooRealVar mu_cb(Form("mu_cb_%d_%d", i+1, j+1), "CB mean", Mucb_ini[i][j], Mucb_lowlim[i][j], Mucb_uplim[i][j]);
            RooRealVar sigma_cb(Form("sigma_cb_%d_%d", i+1, j+1), "CB sigma", Sigmacb_ini[i][j], Sigmacb_lowlim[i][j], Sigmacb_uplim[i][j]);
            RooRealVar alphaL_cb(Form("alphaL_cb_%d_%d", i+1, j+1), "CB alphaL", val_alphL[i], 0, 15.0);
            RooRealVar nL_cb(Form("nL_cb_%d_%d", i+1, j+1), "CB nL", val_nL[i], 0, 500.0);
            RooRealVar alphaR_cb(Form("alphaR_cb_%d_%d", i+1, j+1), "CB alphaR", val_alphR[i], 0, 15.0);
            RooRealVar nR_cb(Form("nR_cb_%d_%d", i+1, j+1), "CB nR", val_nR[i], 0, 500.0);

            // Create Crystal Ball for signal
            RooCrystalBall crystal("crystal", "Crystal Ball signal", mass_data, mu_cb, sigma_cb, alphaL_cb, nL_cb, alphaR_cb, nR_cb);

            // Constrain Crystal Ball tail parameters within Nsigma from MC fit
            alphaL_cb.setRange(val_alphL[i] - Nsigma*inc_alphL[i], val_alphL[i] + Nsigma*inc_alphL[i]);
            nL_cb.setRange(val_nL[i] - Nsigma*inc_nL[i], val_nL[i] + Nsigma*inc_nL[i]);
            alphaR_cb.setRange(val_alphR[i] - Nsigma*inc_alphR[i], val_alphR[i] + Nsigma*inc_alphR[i]);
            nR_cb.setRange(val_nR[i] - Nsigma*inc_nR[i], val_nR[i] + Nsigma*inc_nR[i]);

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
            RooRealVar nsig(Form("nsig_%d_%d", i+1, j+1), "Number of signal events", proj_invmass->Integral() * 0.1, 0, proj_invmass->Integral());
            RooRealVar nbkg(Form("nbkg_%d_%d", i+1, j+1), "Number of background events", proj_invmass->Integral() * 0.9, 0, proj_invmass->Integral());

            // Combined model
            RooAddPdf model("model", "Signal + Background", RooArgList(crystal, background), RooArgList(nsig, nbkg));

            // Set fit range for signal + background fit
            mass_data.setRange("fit_range", LeftLowLim[i][j], RightUpLim[i][j]);

            // Fit full model
            RooFitResult *full_fit_result = model.fitTo(data, RooFit::Range("fit_range"), RooFit::Save(), RooFit::PrintLevel(-1));
            int fitstatus = full_fit_result->status();

            // Plot combined fit
            PlotDataFit(mass_data, data, model, background, crystal, LeftLowLim[i][j], RightUpLim[i][j], i, j, fitstatus);

            //===========================================
            // CALCULATE SCALE AND SMEARING
            //===========================================
            // Extract fitted parameters from data
            double mu_data = mu_cb.getVal();
            double inc_mu_data = mu_cb.getError();
            double sigma_data = sigma_cb.getVal();
            double inc_sigma_data = sigma_cb.getError();

            // Calculate scale correction following FitData.C formula
            double scale = 1 - (mu_data / mu_mc[i]);
            double inc_scale = (mu_data / mu_mc[i]) * sqrt((inc_mu_data / mu_data)*(inc_mu_data / mu_data) + (inc_mu_mc[i] / mu_mc[i])*(inc_mu_mc[i] / mu_mc[i]));

            // Calculate smearing correction following FitData.C formula
            double smearing = sigma_data / val_sigma_mc[i];
            double inc_smearing = smearing * sqrt((inc_sigma_mc[i] / val_sigma_mc[i])*(inc_sigma_mc[i] / val_sigma_mc[i]) + (inc_sigma_data / sigma_data)*(inc_sigma_data / sigma_data));

            // Fill 2D histograms
            h_scale_fullEle->SetBinContent(i+1, j+1, scale);
            h_scale_fullEle->SetBinError(i+1, j+1, inc_scale);
            h_smearing_fullEle->SetBinContent(i+1, j+1, smearing);
            h_smearing_fullEle->SetBinError(i+1, j+1, inc_smearing);

            std::cout << "Pt bin " << i+1 << ", Run bin " << j+1 << ": Scale = " << scale << " ± " << inc_scale 
                      << ", Smearing = " << smearing << " ± " << inc_smearing << std::endl;
        }
    }

    //===========================================
    // STEP 3: PLOT SCALE AND SMEARING
    //===========================================
    std::cout << "\n=== STEP 3: PLOTTING SCALE AND SMEARING ===" << std::endl;

    // Create canvas for scale and smearing plots
    TCanvas *c_scale_smear_fullEle = new TCanvas("c_scale_smear_fullEle", "Scale and Smearing vs p_{T} and Run (Full Electron)", 1600, 600);
    c_scale_smear_fullEle->Divide(2,1);

    // Plot scale
    c_scale_smear_fullEle->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_scale_fullEle->SetTitle("");
    h_scale_fullEle->GetXaxis()->SetTitle("p_{T} [GeV]");
    h_scale_fullEle->GetYaxis()->SetTitle("Run Number");
    h_scale_fullEle->GetZaxis()->SetTitle("1 - m^{(J/#psi)}_{data}/m^{(J/#psi)}_{MC}");
    h_scale_fullEle->SetMarkerSize(1.5);
    h_scale_fullEle->SetMarkerColor(rossoCMS);
    h_scale_fullEle->Draw("COLZ");

    // Add CMS labels to first pad
    TLatex cmsLabel;
    cmsLabel.SetNDC();
    cmsLabel.SetTextFont(62);
    cmsLabel.SetTextSize(0.05);
    cmsLabel.SetTextAlign(11);
    cmsLabel.DrawLatex(0.18, 0.92, "CMS");

    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextFont(52);
    prelimLabel.SetTextSize(0.05);
    prelimLabel.SetTextAlign(11);
    prelimLabel.DrawLatex(0.27, 0.92, "Preliminary");

    TLatex lumiLabel;
    lumiLabel.SetNDC();
    lumiLabel.SetTextFont(42);
    lumiLabel.SetTextSize(0.045);
    lumiLabel.SetTextAlign(31);
    lumiLabel.DrawLatex(0.85, 0.92, "38.01 fb^{-1} (2022, 13.7 TeV)");

    // Plot smearing
    c_scale_smear_fullEle->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_smearing_fullEle->SetTitle("");
    h_smearing_fullEle->GetXaxis()->SetTitle("p_{T} [GeV]");
    h_smearing_fullEle->GetYaxis()->SetTitle("Run Number");
    h_smearing_fullEle->GetZaxis()->SetTitle("#sigma_{data}/#sigma_{MC}");
    h_smearing_fullEle->SetMarkerSize(1.5);
    h_smearing_fullEle->SetMarkerColor(violaCMS);
    h_smearing_fullEle->Draw("COLZ TEXT");

    // Add CMS labels to second pad
    cmsLabel.DrawLatex(0.18, 0.92, "CMS");
    prelimLabel.DrawLatex(0.27, 0.92, "Preliminary");
    lumiLabel.DrawLatex(0.85, 0.92, "38.01 fb^{-1} (2022, 13.7 TeV)");

    // Save the canvas
    c_scale_smear_fullEle->SaveAs("PlotConID2022/ScaleonFullEle/ScaleAndSmearing_FullElectron_vsPt_runN.png");
    delete c_scale_smear_fullEle;

    // Save histograms to file
    TFile *file_corrections_fullEle = new TFile("scale_corrections_fullEle.root", "RECREATE");
    h_scale_fullEle->Write();
    h_smearing_fullEle->Write();
    file_corrections_fullEle->Close();
    delete file_corrections_fullEle;

    std::cout << "\n=== ANALYSIS COMPLETED SUCCESSFULLY ===" << std::endl;
    std::cout << "Results saved in:" << std::endl;
    std::cout << "- Individual fit plots: PlotConID2022/ScaleonFullEle/Pt_bin%d/" << std::endl;
    std::cout << "- Scale and Smearing plot: PlotConID2022/ScaleonFullEle/ScaleAndSmearing_FullElectron_vsPt_runN.png" << std::endl;
    std::cout << "- Scale and Smearing histograms: scale_corrections_fullEle.root" << std::endl;

    // Close files
    dataFile->Close();
    mcFile->Close();
    delete dataFile;
    delete mcFile;
}
