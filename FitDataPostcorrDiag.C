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

#define NbinsPt 6
#define NbinsRun 9
#define Nsigma 1

// Define CMS colors
    Color_t bluCMS = TColor::GetColor("#5790fc");
    Color_t gialloCMS = TColor::GetColor("#ffcc00");
    Color_t violaCMS = TColor::GetColor("#660099");
    Color_t rossoCMS = TColor::GetColor("#cc0000");
    Color_t arancioneCMS = TColor::GetColor("#ca5200");

void PlotBackgroundFit(RooRealVar& mass, RooDataHist& data, RooAddPdf& background, int i, int j) {
    // Create a canvas to draw the background fit
    TCanvas* cBackgroundFit = new TCanvas(Form("cBackgroundFit_%d_%d", i+1, j+1), "Background Fit", 800, 600);

    // Create a frame for the mass variable to hold the plot
    RooPlot* frame = mass.frame();
    
    // Plot the data points on the frame
    data.plotOn(frame, RooFit::Name("Data"), RooFit::MarkerStyle(20), RooFit::MarkerSize(0.7));

    // Plot the background fit function
    background.plotOn(frame, RooFit::Name("BackgroundFit"), RooFit::LineColor(kMagenta), RooFit::LineStyle(kDashed));

    // Add labels and titles
    frame->SetTitle(Form("Background Fit for Pt Bin %d and Run Bin %d", i+1, j+1));
    frame->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle("Events");

    // Draw the frame on the canvas
    frame->Draw();

    // Add legend for clarity
    TLegend* legend = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend->AddEntry(frame->findObject("Data"), "Data", "p");
    legend->AddEntry(frame->findObject("BackgroundFit"), "Background Fit", "l");
    legend->Draw();

    // Estrai i parametri di interesse dal modello per la stampa
    RooRealVar* mu = (RooRealVar*)background.getVariables()->find(Form("gauss_mu_%d", i+1));
    RooRealVar* sigma = (RooRealVar*)background.getVariables()->find(Form("gauss_sigma_%d", i+1));

    //background.getVariables()->Print("v");

    // Optional: add text annotations for fit parameters
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.15, 0.85, Form("Fit for Pt bin %d, Run bin %d", i+1, j+1));
    latex.DrawLatex(0.15, 0.85, Form("Fit for Pt bin %d", i+1));
    latex.DrawLatex(0.15, 0.8, Form("#mu = %.4f +/- %.4f, #sigma = %.4f +/- %.4f", mu->getVal(), mu->getError(), sigma->getVal(), sigma->getError()));

    // Save the canvas as a file (optional)
    if(j == 1000) {
        // Create a canvas for inclusive run fit
        cBackgroundFit->SetTitle(Form("Background Fit for Pt Bin %d (All Runs)", i+1));
        // Save with special name for inclusive run
        cBackgroundFit->SaveAs(Form("PlotConID2022/FitCorrectedDataDiag_id2022/Pt_bin%d/BackgroundFit_Pt%d_AllRuns.png", i+1, i+1));
    } else if(j == 999) {
        // Keep existing special case for j=999
        cBackgroundFit->SaveAs(Form("PlotConID2022/FitCorrectedDataDiag_id2022/Pt_bin%d/BackgroundFit_Pt%d_Run%d.png", i+1, i+1, j+1));
    } else {
        // Regular case - save with run bin number
        cBackgroundFit->SaveAs(Form("PlotConID2022/FitCorrectedDataDiag_id2022/Pt_bin%d/BackgroundFit_Pt%d_Run%d.png", i+1, i+1, j+1));
    }
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
    legend->AddEntry(frame->findObject("Data"), "Data", "p");
    legend->AddEntry(frame->findObject("Background"), "Background Fit", "l");
    legend->AddEntry(frame->findObject("CrystalBall"), "Signal (Crystal Ball)", "l");
    legend->AddEntry(frame->findObject("CombinedFit"), "Combined Fit", "l");
    legend->Draw();

    // Estrai i parametri di interesse dal modello per la stampa
    RooRealVar* mu = (RooRealVar*)model.getVariables()->find(Form("mu_cb_%d", i+1));
    RooRealVar* sigma = (RooRealVar*)model.getVariables()->find(Form("sigma_cb_%d", i+1));
    double chi2 = frame->chiSquare();

    // Optional: Add annotations for fit information
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    //latex.DrawLatex(0.15, 0.85, Form("Fit for Pt bin %d, Run bin %d", i+1, j+1));
    /*if(fitstatus == 0){
        latex.DrawLatex(0.15, 0.8, "Fit converged");
    }else{
        latex.DrawLatex(0.15, 0.8, Form("Fit did not converge. Status %d", fitstatus));
    }*/
    latex.DrawLatex(0.15, 0.75, Form("#mu = %.4f #pm %.4f GeV", mu->getVal(), mu->getError()));
    latex.DrawLatex(0.15, 0.7, Form("#sigma = %.4f #pm %.4f GeV", sigma->getVal(), sigma->getError())); 
    //latex.DrawLatex(0.15, 0.66, Form("#chi^{2} = %.2f", chi2));

    // Save plot as image (optional)
    if(j == 1000) {
        // Set title for inclusive run
        frame->SetTitle(Form("Data Fit for Pt Bin %d (All Runs)", i+1));
        // Save with special name for inclusive run
        c->SaveAs(Form("PlotConID2022/FitCorrectedDataDiag_id2022/Pt_bin%d/DataFit_Pt%d_AllRuns.png", i+1, i+1));
    } else if(j == 999) {
        // Keep existing special case
        c->SaveAs(Form("PlotConID2022/FitCorrectedDataDiag_id2022/Pt_bin%d/DataFit_Pt%d_Run%d.png", i+1, i+1, j+1));
    } else {
        // Regular case
        c->SaveAs(Form("PlotConID2022/FitCorrectedDataDiag_id2022/Pt_bin%d/DataFit_Pt%d_Run%d.png", i+1, i+1, j+1));
    }
    delete c;
}


void FitDataPostcorrDiag(){
///////////////
    // Apri il file ROOT contenente gli istogrammi
    TFile *file = TFile::Open("outputHistograms_DATA_partF.root");

    // Definisci i limiti personalizzati e fattore di rebinning
    double LeftLowLim[NbinsPt][NbinsRun] = {
        {1.6, 1.4, 1.4, 1.6, 1.7, 1.6, 1.5, 1.6, 1.5}, 
        {1.6, 1.4, 1.4, 1.6, 1.7, 1.4, 1.4, 1.6, 1.7}, 
        {1.6, 1.3, 1.4, 1.6, 1.7, 1.6, 1.4, 1.6, 1.5}, 
        {1.6, 1.4, 1.4, 1.5, 1.7, 1.6, 1.5, 1.6, 1.4},
        {1.5, 1.5, 1.5, 1.5, 1.5, 1.6, 1.6, 1.6, 1.5}, 
        {1.8, 1.8, 1.8, 1.8, 1.8, 2, 1.5, 2, 2.45}
    };  // Limiti sinistri personalizzati
    
    double LeftUpLim[NbinsPt][NbinsRun] = {
        {2.4, 2.5, 2.2, 2.2, 2.5, 2.4, 2.3, 2.4, 2.4}, 
        {2.4, 2.5, 2.2, 2.2, 2.5, 2.5, 2.3, 2.4, 2.5}, 
        {2.4, 2.4, 2.2, 2.2, 2.5, 2.4, 2.45, 2.4, 2.45}, 
        {2.4, 2.5, 2.2, 2.45, 2.5, 2.5, 2.4, 2.4, 2.4},
        {2.5, 2.5, 2.4, 2.5, 2.6, 2.4, 2.6, 2.4, 2.5}, 
        {2.5, 2.5, 2.5, 2.5, 2.5, 2.4, 2.3, 2.4, 2.8}
    };
    
    double RightLowLim[NbinsPt][NbinsRun] = {
        {3.6, 3.4, 3.6, 3.6, 3.6, 3.5, 3.6, 3.5, 3.5}, 
        {3.6, 3.4, 3.6, 3.6, 3.6, 3.5, 3.6, 3.5, 3.6}, 
        {3.6, 3.4, 3.5, 3.5, 3.5, 3.5, 3.6, 3.5, 3.6}, 
        {3.6, 3.4, 3.5, 3.6, 3.6, 3.5, 3.6, 3.5, 3.5}, 
        {3.5, 3.55, 3.55, 3.55, 3.55, 3.55, 3.5, 3.5, 3.5}, 
        {3.6, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6}
    }; // Limiti destri personalizzati
    
    double RightUpLim[NbinsPt][NbinsRun] = {
        {5, 5, 5, 5, 5, 5, 5, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5.2, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5, 5, 5}, 
        {5, 5.2, 5, 5, 5, 5, 5, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5, 5, 5}
    };

    //limiti inclusivi in run number (aggiungere)
    double LeftLowLim_incl[NbinsPt] = {1.3, 1.3, 1.6, 1.6, 1.8, 2.1};
    double LeftUpLim_incl[NbinsPt] = {2.45, 2.45, 2.45, 2.4, 2.6, 2.7};
    double RightLowLim_incl[NbinsPt] = {3.65, 3.5, 3.55, 3.55, 3.5, 3.55};
    double RightUpLim_incl[NbinsPt] = {5.3, 5.2, 4.8, 5, 5.2, 5.2};
    

    double BackgroundYlims[NbinsPt][NbinsRun] = {
        {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000}, 
        {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000}, 
        {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000},
        {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000}, 
        {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000}, 
        {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000}
    };

    // Array per parametri crystal ball
    ///////////////// Parametri Crystal Ball
    double Mucb_ini[NbinsPt][NbinsRun] = {
        {3.0442, 3.0337, 3.0422, 3.0489, 3.0, 3.0400, 3.0450, 3.0500, 3.0550},
        {3.0199, 3.0337, 3.0422, 2.9975, 3.0782, 3.0550, 3.0600, 3.0650, 3.0700},
        {3.0442, 3.0337, 3.0350, 3.0489, 3.0782, 3.0600, 3.0650, 3.0700, 3.0750},
        {3.0442, 3.0337, 3.0422, 3.0489, 3.0782, 3.0650, 3.0700, 3.0750, 3.0800},
        {3.0442, 3.0337, 3.0422, 3.0489, 3.1098, 3.0700, 3.0750, 3.0800, 3.0850},
        {3.0442, 3.0337, 3.0422, 3.0489, 3.0782, 3.0750, 3.0800, 3.0850, 3.0900}
    };
    
    double Mucb_lowlim[NbinsPt][NbinsRun] = {
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.9, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.9, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8}
    };
    
    double Mucb_uplim[NbinsPt][NbinsRun] = {
        {3.2, 3.5, 3.5, 3.5, 3.2, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.5, 3.2, 3.5, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.2, 3.5, 3.5, 3.3, 3.4, 3.3, 3.2},
        {3.1, 3.5, 3.5, 3.5, 3.5, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.5, 3.5, 3.5, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.5, 3.5, 3.5, 3.3, 3.4, 3.3, 3.2}
    };
    
    double Sigmacb_ini[NbinsPt][NbinsRun] = {
        {0.1637, 0.1503, 0.1351, 0.1271, 0.15, 0.14, 0.13, 0.14, 0.15},
        {0.1333, 0.1413, 0.1351, 0.1271, 0.1115, 0.1342, 0.12, 0.13, 0.1452},
        {0.1637, 0.1503, 0.1337, 0.1271, 0.1115, 0.12, 0.11, 0.12, 0.13},
        {0.1637, 0.1503, 0.1351, 0.1271, 0.1115, 0.11, 0.10, 0.11, 0.12},
        {0.1637, 0.1503, 0.1351, 0.1271, 0.1115, 0.11, 0.10, 0.11, 0.12},
        {0.1637, 0.1503, 0.1351, 0.1271, 0.1115, 0.10, 0.09, 0.10, 0.11}
    };
    
    double Sigmacb_uplim[NbinsPt][NbinsRun] = {
        {0.3, 0.5, 0.5, 0.5, 0.5, 0.4, 0.4, 0.4, 0.4},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.2, 0.3, 0.3, 0.3},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.3, 0.3, 0.3, 0.3},
        {0.2, 0.5, 0.5, 0.5, 0.15, 0.3, 0.3, 0.3, 0.3},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.3, 0.3, 0.3, 0.3},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.3, 0.3, 0.3, 0.3}
    };
    
    double Sigmacb_lowlim[NbinsPt][NbinsRun] = {
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}
    };

    //parametri inclusivi in pt
    double Mucb_ini_incl[NbinsPt] = {3.0470, 3.0421, 3.0343, 3.0714, 3.0867, 3.1217};
    double Mucb_lowlim_incl[NbinsPt] = {2.9, 2.9, 2.5, 2.5, 2.8, 3};
    double Mucb_uplim_incl[NbinsPt] = {3.2, 3.2, 3.2, 3.2, 3.2, 3.2};
    double Sigmacb_ini_incl[NbinsPt] = {0.1671, 0.149, 0.1440, 0.1203, 0.1267, 0.09};
    double Sigmacb_uplim_incl[NbinsPt] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.2};
    double Sigmacb_lowlim_incl[NbinsPt] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05};


    // Aggiungi gli array per i parametri della Gaussiana
    double gauss_mu_init[NbinsPt][NbinsRun] = {
        {3.6, 3.4, 3.7, 3.4, 3.5, 3.55, 3.6, 3.55, 3.5},
        {3.6, 3.4, 3.7, 3.4, 3.4, 3.5, 3.55, 3.5, 3.45},
        {3.6, 3.4, 3.7, 3.4, 3.4, 3.5, 3.55, 3.5, 3.45},
        {3.6, 3.4, 3.7, 3.4, 3.4, 3.5, 3.55, 3.5, 3.45},
        {3.6, 3.65, 3.7, 3.4, 3.6, 3.5, 3.55, 3.5, 3.45},
        {3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 3.6, 3.65, 3.65}
    };  // Valori iniziali per mu della Gaussiana
    
    double gauss_mu_low[NbinsPt][NbinsRun] = {
        {3.5, 3.2, 3.6, 3.2, 3.4, 3.4, 3.5, 3.4, 3.4},
        {3.5, 3.2, 3.6, 3.2, 3.2, 3.4, 3.45, 3.4, 3.3},
        {3.5, 3.2, 3.6, 3.2, 3.2, 3.4, 3.45, 3.4, 3.3},
        {3.5, 3.2, 3.6, 3.2, 3.2, 3.4, 3.45, 3.4, 3.3},
        {3.5, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45},
        {3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.5, 3.55, 3.55}
    };   // Limiti inferiori per mu
    
    double gauss_mu_up[NbinsPt][NbinsRun] = {
        {3.7, 3.75, 3.8, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.75, 3.8, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.75, 3.8, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.75, 3.8, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7},
        {3.7, 3.8, 3.8, 3.7, 3.8, 3.7, 3.7, 3.7, 3.7},
        {3.75, 3.75, 3.8, 3.75, 3.75, 3.75, 3.7, 3.75, 3.75}
    };    // Limiti superiori per mu
    
    double gauss_sigma_init[NbinsPt][NbinsRun] = {
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1457},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.1, 0.0959, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
        {0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15}
    }; // Valori iniziali per sigma della Gaussiana
    
    double gauss_sigma_low[NbinsPt][NbinsRun] = {
        {0.05, 0.05, 0.02, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.02, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.02, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.05, 0.05, 0.02, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.1, 0.05, 0.02, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}
    }; // Limiti inferiori per sigma
    
    double gauss_sigma_up[NbinsPt][NbinsRun] = {
        {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2},
        {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2},
        {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2},
        {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2},
        {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2},
        {0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18}
    }; // Limiti superiori per sigma

    //parametri gaussiana inclusivi in pT
    double gauss_mu_init_incl[NbinsPt] = {3.55, 3.55, 3.5, 3.62, 3.65, 3.6195};
    double gauss_mu_low_incl[NbinsPt] = {3.5, 3.4, 3.4, 3.5, 3.62, 3.5};
    double gauss_mu_up_incl[NbinsPt] = {3.7, 3.7, 3.7, 3.7, 3.8, 3.7};

    double bincenters[NbinsPt] = {5.5, 8, 10, 12.5, 17, 30};
    double binhalfwidths[NbinsPt] = {1.5, 1, 1, 1.5, 3, 10};

    // Leggi i parametri della Crystal Ball dal file
    TFile *file_param = TFile::Open("fit_results.root", "READ");
    TTree *tree = (TTree*)file_param->Get("fitResults");
    if (!tree) {
        std::cerr << "Impossibile trovare il TTree fitResults." << std::endl;
        return;
    }

    // Variabili per contenere i valori letti dal TTree
    int bin;
    double mu_mc, inc_mu_mc, sigma_ini, inc_sigma, nL_ini, inc_nL, alphaL_ini, inc_alphaL, nR_ini, inc_nR, alphaR_ini, inc_alphaR;
    // Collega le variabili ai rami del TTree
    tree->SetBranchAddress("Bin", &bin);
    tree->SetBranchAddress("Mean", &mu_mc);
    tree->SetBranchAddress("MeanError", &inc_mu_mc);
    tree->SetBranchAddress("Sigma", &sigma_ini);
    tree->SetBranchAddress("SigmaError", &inc_sigma);
    tree->SetBranchAddress("nL", &nL_ini);
    tree->SetBranchAddress("nLError", &inc_nL);
    tree->SetBranchAddress("alphaL", &alphaL_ini);
    tree->SetBranchAddress("alphaLError", &inc_alphaL);
    tree->SetBranchAddress("nR", &nR_ini);
    tree->SetBranchAddress("nRError", &inc_nR);
    tree->SetBranchAddress("alphaR", &alphaR_ini);
    tree->SetBranchAddress("alphaRError", &inc_alphaR);

    //Istogrammi 2D binnati su Pt e run number che conterranno i valori delle correzioni di scala
    double Ptbins[] = {4, 7, 9, 11, 14, 20, 40}; 
    double runBins[] = {356309, 356900, 357538, 357732, 360000, 360400, 361000, 361600, 362200, 362760}; //9 bins

    //Estraggo gli istogrammi e proietto su ciascun bin di Pt
    
    // Extract the 3D histograms from the input file
    TH3D *h_invMass_ECAL_corrected = (TH3D*)file->Get("h_invMass_ECAL_corrected");
    TH3D *h_invMass_ECAL_check = (TH3D*)file->Get("h_invMass_ECAL_check");
    //TH3D *h_invMass_ECAL_corr_offdiag = (TH3D*)file->Get("h_invMass_ECAL_corr_offdiag");
    TH3D *h_invMass_ECAL_check_offdiag = (TH3D*)file->Get("h_invMass_ECAL_check_offdiag");
    TH3D *h_invMass_ECAL_corr_diag = (TH3D*)file->Get("h_invMass_ECAL_corr_diag");
    TH3D *h_invMass_ECAL_check_diag = (TH3D*)file->Get("h_invMass_ECAL_check_diag");

    // Check if histograms were successfully retrieved
    if (!h_invMass_ECAL_corrected || !h_invMass_ECAL_check || 
        //!h_invMass_ECAL_corr_offdiag || 
        !h_invMass_ECAL_check_offdiag ||
        !h_invMass_ECAL_corr_diag || !h_invMass_ECAL_check_diag) {
        std::cerr << "Error: One or more required histograms not found in the file." << std::endl;
        return;
    }

    TH1D *h_smearing_postcorr = new TH1D("h_smearing_postcorr", "Smearing after correction", NbinsPt, Ptbins);

    // Create directory for output plots
    gSystem->Exec("mkdir -p PlotConID2022/FitCorrectedDataDiag_id2022");
    for (int i = 0; i < NbinsPt; i++) {
        gSystem->Exec(Form("mkdir -p PlotConID2022/FitCorrectedDataDiag_id2022/Pt_bin%d", i+1));
    }
    
    //....oooOO0OOooo........oooOO0OOooo.... Loop su tutti gli istogrammi ....oooOO0OOooo........oooOO0OOooo....

    for(int i=0; i < NbinsPt; i++){
        // Set the pt bin range (i+1 because bin counting starts at 1)
        int binLow = i+1;
        int binHigh = i+1;
        
        // Create proper histogram names with pt bin information
        TString checkName = Form("h2D_runN_invM_check_ptbin%d", i+1);
        TString checkDiagName = Form("h2D_runN_invM_check_diag_ptbin%d", i+1);
        TString checkOffdiagName = Form("h2D_runN_invM_check_offdiag_ptbin%d", i+1);
        TString corrDiagName = Form("h2D_runN_invM_corr_diag_ptbin%d", i+1);
        TString corrName = Form("h2D_runN_invM_corr_ptbin%d", i+1);
        TString corrOffdiagName = Form("h2D_runN_invM_corr_offdiag_ptbin%d", i+1);

        //Create a histogram for the scale 
            // Definisci gli istogrammi 2D per le correzioni di scala
        TH1D *h_scale_postcorr = new TH1D(Form("h_scale_postcorr%d", i+1), "Scale between data and MC post corrections -all categories; Pt ; Run number", NbinsRun, runBins);
        
        // Get the projections for this pt bin
        h_invMass_ECAL_check->GetXaxis()->SetRange(binLow, binHigh);
        TH2D *h2D_runN_invM_check = (TH2D*)h_invMass_ECAL_check->Project3D("yz");
        h2D_runN_invM_check->SetName(checkName);
        
        h_invMass_ECAL_check_diag->GetXaxis()->SetRange(binLow, binHigh);
        TH2D *h2D_runN_invM_check_diag = (TH2D*)h_invMass_ECAL_check_diag->Project3D("yz");
        h2D_runN_invM_check_diag->SetName(checkDiagName);
        
        h_invMass_ECAL_check_offdiag->GetXaxis()->SetRange(binLow, binHigh);
        TH2D *h2D_runN_invM_check_offdiag = (TH2D*)h_invMass_ECAL_check_offdiag->Project3D("yz");
        h2D_runN_invM_check_offdiag->SetName(checkOffdiagName);
        
        h_invMass_ECAL_corr_diag->GetXaxis()->SetRange(binLow, binHigh);
        TH2D *h2D_runN_invM_corr_diag = (TH2D*)h_invMass_ECAL_corr_diag->Project3D("yz");
        h2D_runN_invM_corr_diag->SetName(corrDiagName);
        
        h_invMass_ECAL_corrected->GetXaxis()->SetRange(binLow, binHigh);
        TH2D *h2D_runN_invM_corr = (TH2D*)h_invMass_ECAL_corrected->Project3D("yz");
        h2D_runN_invM_corr->SetName(corrName);
        
        //h_invMass_ECAL_corr_offdiag->GetXaxis()->SetRange(binLow, binHigh);
        //TH2D *h2D_runN_invM_corr_offdiag = (TH2D*)h_invMass_ECAL_corr_offdiag->Project3D("yz");
        //h2D_runN_invM_corr_offdiag->SetName(corrOffdiagName);

        for(int j=0; j < NbinsRun; j++){

            // Proietto l'istogramma 2D su un istogramma 1D per il run number
            TH1D *hist_corr = (TH1D*)h2D_runN_invM_corr->ProjectionX(Form("h1D_runN_invM_corr_ptbin%d_runbin%d", i+1, j+1), j+1, j+1);
            TH1D *hist_corr_diag = (TH1D*)h2D_runN_invM_corr_diag->ProjectionX(Form("h1D_runN_invM_corr_diag_ptbin%d_runbin%d", i+1, j+1), j+1, j+1);
            //TH1D *hist_corr_offdiag = (TH1D*)h2D_runN_invM_corr_offdiag->ProjectionY(Form("h1D_runN_invM_corr_offdiag_ptbin%d_runbin%d", i+1, j+1), j+1, j+1);

            
            //hist->Rebin(rebin_factor[i][j]);

            //....oooOO0OOooo.... fit fondo ....oooOO0OOooo....
            // Definisci la variabile di massa invariante
            RooRealVar mass("mass", "m(e^{+}e^{-})", 0, 6); 

            // Converto l'istogramma in un RooDataHist
            RooDataHist data("data", "Dataset from histogram", mass, hist_corr_diag);

            // Definisci i parametri del polinomio di secondo grado
            RooRealVar A(Form("A_%d_%d", i+1, j+1), "4th deg coeff", 0,  -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar B(Form("B_%d_%d", i+1, j+1), "3rd deg coeff", 0,  -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar C(Form("C_%d_%d", i+1, j+1), "2nd deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar D(Form("D_%d_%d", i+1, j+1), "1st deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar E(Form("E_%d_%d", i+1, j+1), "0 deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
            RooPolynomial poly("poly", "Polynomial of 4th degree", mass, RooArgList(A, B, C, D, E));

            // Definisci i parametri della Gaussiana
            RooRealVar gauss_mu(Form("gauss_mu_%d", i+1), "Gaussian mean", gauss_mu_init[i][j], gauss_mu_low[i][j], gauss_mu_up[i][j]);
            RooRealVar gauss_sigma(Form("gauss_sigma_%d", i+1), "Gaussian sigma", gauss_sigma_init[i][j], gauss_sigma_low[i][j], gauss_sigma_up[i][j]);
            RooGaussian gauss("gauss", "Gaussian component", mass, gauss_mu, gauss_sigma);

            // Definisci il modello di fondo combinando il polinomio con la Gaussiana
            RooRealVar frac_gauss("frac_gauss", "fraction of Gaussian", 0.3, 0.0, 1.0);
            RooAddPdf background("background", "Background Model", RooArgList(poly, gauss), RooArgList(frac_gauss));

            //Definisci il range per il fit del fondo
            mass.setRange("range1", LeftLowLim[i][j], LeftUpLim[i][j]);
            mass.setRange("range2", RightLowLim[i][j], RightUpLim[i][j]);

            //////////////////// Esegui fit sul solo fondo //////////////////////
            RooFitResult *fit_result = background.fitTo(data, RooFit::Range("range1,range2"),RooFit::MaxCalls(10000000), RooFit::Save(), RooFit::SumW2Error(true));

            // Ottieni i valori dei parametri e i loro errori dal primo fit
            double A_val = A.getVal(), A_err = A.getError();
            double B_val = B.getVal(), B_err = B.getError();
            double C_val = C.getVal(), C_err = C.getError();
            double D_val = D.getVal(), D_err = D.getError();
            double E_val = E.getVal(), E_err = E.getError();
            double gauss_mu_val = gauss_mu.getVal(), gauss_mu_err = gauss_mu.getError();
            double gauss_sigma_val = gauss_sigma.getVal(), gauss_sigma_err = gauss_sigma.getError();

            //////////////////////////////////////////////////////// Plot solo fondo
            PlotBackgroundFit(mass, data, background, i, j);

            //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
            //                                          FIT SEGNALE + FONDO
            //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

            mass.setRange("range_full", LeftLowLim[i][j], RightUpLim[i][j]);

            // Leggi i parametri della Crystal Ball dal TTree per il bin corrente
            tree->GetEntry(i);

            RooRealVar mu_cb(Form("mu_cb_%d", i+1), "mu_cb", Mucb_ini[i][j], Mucb_lowlim[i][j], Mucb_uplim[i][j]);
            RooRealVar sigma_cb(Form("sigma_cb_%d", i+1), "sigma_cb", Sigmacb_ini[i][j], Sigmacb_lowlim[i][j], Sigmacb_uplim[i][j]);
            RooRealVar alphaL_cb(Form("alphaL_cb_%d", i+1), "alphaL_cb", alphaL_ini, alphaL_ini - Nsigma*inc_alphaL, alphaL_ini + Nsigma*inc_alphaL);
            RooRealVar nL_cb(Form("nL_cb_%d", i+1), "nL_cb", nL_ini, nL_ini - Nsigma*inc_nL, nL_ini + Nsigma*inc_nL);
            RooRealVar alphaR_cb(Form("alphaR_cb_%d", i+1), "alphaR_cb", alphaR_ini, alphaR_ini - Nsigma*inc_alphaR, alphaR_ini + Nsigma*inc_alphaR);
            RooRealVar nR_cb(Form("nR_cb_%d", i+1), "nR_cb", nR_ini, nR_ini - Nsigma*inc_nR, nR_ini + Nsigma*inc_nR);

            if (i == 4) {
                // Use wider ranges for this specific pt bin (2*Nsigma instead of Nsigma)
                alphaL_cb.setRange(alphaL_ini - 2*Nsigma*inc_alphaL, alphaL_ini + 2*Nsigma*inc_alphaL);
                nL_cb.setRange(nL_ini - 2*Nsigma*inc_nL, nL_ini + 2*Nsigma*inc_nL);
                alphaR_cb.setRange(alphaR_ini - 3*Nsigma*inc_alphaR, alphaR_ini + 3*Nsigma*inc_alphaR);
                nR_cb.setRange(nR_ini - 3*Nsigma*inc_nR, nR_ini + 3*Nsigma*inc_nR);
            }


            // Definisci la Crystal Ball
            RooCrystalBall crystal("crystal", "crystal ball", mass, mu_cb, sigma_cb, alphaL_cb, nL_cb, alphaR_cb, nR_cb);

            // Aggiorna i parametri della Gaussiana con i vincoli per il secondo fit
            gauss_mu.setVal(gauss_mu_val);
            gauss_sigma.setVal(gauss_sigma_val);
            //gauss_mu.setRange(gauss_mu_val - Nsigma*gauss_mu_err, gauss_mu_val + Nsigma*gauss_mu_err);
            gauss_mu.setConstant(true);
            gauss_sigma.setRange(gauss_sigma_val - Nsigma*gauss_sigma_err, gauss_sigma_val + Nsigma*gauss_sigma_err);
            A.setVal(A_val);
            A.setRange(A_val - Nsigma * A_err, A_val + Nsigma * A_err);
            B.setVal(B_val);
            B.setRange(B_val - Nsigma * B_err, B_val + Nsigma * B_err);
            C.setVal(C_val);
            C.setRange(C_val - Nsigma * C_err, C_val + Nsigma * C_err);
            D.setVal(D_val);
            D.setRange(D_val - Nsigma * D_err, D_val + Nsigma * D_err);
            E.setVal(E_val);
            E.setRange(E_val - Nsigma * E_err, E_val + Nsigma * E_err);

            RooRealVar frac("frac", "fraction of background", 0.5, 0.0, 1.0);
            RooAddPdf model("model", "signal + background", RooArgList(crystal, background), RooArgList(frac));

            // Esegui il fit segnale + fondo
            fit_result = model.fitTo(data, RooFit::Range("range_full"),RooFit::MaxCalls(10000000), RooFit::Save(), RooFit::SumW2Error(true));
            //Plot del fit completo
            PlotDataFit(mass, data, model, background, crystal, LeftLowLim[i][j], RightUpLim[i][j], i, j, fit_result->status());

            //Calcolo la scala 
            double mu_data = mu_cb.getVal();
            double inc_mu_data = mu_cb.getError();
            double scale_postcorr = 1 - (mu_data / mu_mc);
            double inc_scale_postcorr = (mu_data / mu_mc) * sqrt((inc_mu_data / mu_data)*(inc_mu_data / mu_data) + (inc_mu_mc / mu_mc)*(inc_mu_mc / mu_mc));
            //Scrivo i parametri di interesse negli istogrammi 2D
            h_scale_postcorr->SetBinContent(j+1, scale_postcorr);
            h_scale_postcorr->SetBinError(j+1, inc_scale_postcorr);



        } //FINE LOOP SU RUN NUMBER

        //Confronto la scala prima e dopo la correzione in ogni bin di Pt

        //Estraggo la scala prima della correzione
        // Open the file containing scale corrections
        TFile *scaleFile = TFile::Open("scale_corrections.root");
        if (!scaleFile || scaleFile->IsZombie()) {
            std::cerr << "Error: Could not open scale_corrections.root" << std::endl;
            return;
        }

        // Extract the h_scale histogram
        TH2D *h_scale = (TH2D*)scaleFile->Get("h_scale");
        if (!h_scale) {
            std::cerr << "Error: Could not find h_scale histogram in scale_corrections.root" << std::endl;
            scaleFile->Close();
            return;
        }

        // Create a TCanvas to compare scales before and after correction
        TCanvas *cCompare = new TCanvas(Form("cCompare_pt%d", i+1), Form("Scale Comparison for Pt Bin %d", i+1), 900, 600);
        gStyle->SetOptStat(0); // Turn off statistics box

        // Setup the histogram for before correction
        TH1D *h_scale_before = (TH1D*)h_scale->ProjectionY(Form("h_scale_before_pt%d", i+1), i+1, i+1);
        h_scale_before->SetTitle(Form("Scale Comparison - Pt Bin %d (%0.1f < p_{T} < %0.1f GeV)", 
                   i+1, bincenters[i]-binhalfwidths[i], bincenters[i]+binhalfwidths[i]));
        h_scale_before->GetXaxis()->SetTitle("Run Number");
        h_scale_before->GetYaxis()->SetTitle("Scale Factor");
        h_scale_before->SetMarkerStyle(20);
        h_scale_before->SetMarkerColor(kBlue);
        h_scale_before->SetLineColor(kBlue);

        // Setup the histogram for after correction
        h_scale_postcorr->SetMarkerStyle(21);
        h_scale_postcorr->SetMarkerColor(kRed);
        h_scale_postcorr->SetLineColor(kRed);

        // Set y-axis range with maximum at 0.02
        double minY = h_scale_before->GetMinimum();
        //minY = std::min(minY, h_scale_postcorr->GetMinimum());
        h_scale_before->GetYaxis()->SetRangeUser(-0.04, 0.02);

        // Draw both histograms
        h_scale_before->Draw("PE");
        h_scale_postcorr->Draw("PE SAME");

        // Create and draw the legend
        TLegend *legend = new TLegend(0.65, 0.15, 0.89, 0.29);
        legend->AddEntry(h_scale_before, "Before Correction", "lp");
        legend->AddEntry(h_scale_postcorr, "After Correction", "lp");
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->Draw();
        
        // Save the canvas
        cCompare->SaveAs(Form("PlotConID2022/FitCorrectedDataDiag_id2022/Pt_bin%d/Scale_Comparison_Pt%d.png", i+1, i+1));
        delete cCompare;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//                                      Add inclusive run fit for each pt bin
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Project h2D_runN_invM_corr_diag to get inclusive run histogram
TH1D *hist_incl_run = (TH1D*)h2D_runN_invM_corr_diag->ProjectionX(Form("h1D_invM_corr_diag_ptbin%d_inclRun", i+1));
if (!hist_incl_run) {
    std::cerr << "Error: Could not create inclusive run projection for pt bin " << i+1 << std::endl;
    continue;
}

// Now fit this histogram with the same model used for individual run bins
RooRealVar mass("mass", "m(e^{+}e^{-})", 0, 6); 

// Convert histogram to RooDataHist
RooDataHist data_incl_run("data_incl_run", "Dataset from inclusive run histogram", mass, hist_incl_run);

// Define background model parameters
RooRealVar A(Form("A_incl_%d", i+1), "4th deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
RooRealVar B(Form("B_incl_%d", i+1), "3rd deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
RooRealVar C(Form("C_incl_%d", i+1), "2nd deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
RooRealVar D(Form("D_incl_%d", i+1), "1st deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
RooRealVar E(Form("E_incl_%d", i+1), "0 deg coeff", 0, -RooNumber::infinity(), RooNumber::infinity());
RooPolynomial poly("poly", "Polynomial of 4th degree", mass, RooArgList(A, B, C, D, E));

// Define Gaussian parameters using inclusive values
RooRealVar gauss_mu(Form("gauss_mu_%d", i+1), "Gaussian mean", gauss_mu_init_incl[i], gauss_mu_low_incl[i], gauss_mu_up_incl[i]);
RooRealVar gauss_sigma(Form("gauss_sigma_%d", i+1), "Gaussian sigma", gauss_sigma_init[i][0], gauss_sigma_low[i][0], gauss_sigma_up[i][0]);
RooGaussian gauss("gauss", "Gaussian component", mass, gauss_mu, gauss_sigma);

// Combine into background model
RooRealVar frac_gauss("frac_gauss", "fraction of Gaussian", 0.3, 0.0, 1.0);
RooAddPdf background("background", "Background Model", RooArgList(poly, gauss), RooArgList(frac_gauss));

// Define fit ranges for background
mass.setRange("range1", LeftLowLim_incl[i], LeftUpLim_incl[i]);
mass.setRange("range2", RightLowLim_incl[i], RightUpLim_incl[i]);

// Perform background-only fit
RooFitResult *fit_result_incl = background.fitTo(data_incl_run, 
                                               RooFit::Range("range1,range2"),
                                               RooFit::MaxCalls(10000000), 
                                               RooFit::Save(), 
                                               RooFit::SumW2Error(true));

// Extract parameter values from background fit
double A_val = A.getVal(), A_err = A.getError();
double B_val = B.getVal(), B_err = B.getError();
double C_val = C.getVal(), C_err = C.getError();
double D_val = D.getVal(), D_err = D.getError();
double E_val = E.getVal(), E_err = E.getError();
double gauss_mu_val = gauss_mu.getVal(), gauss_mu_err = gauss_mu.getError();
double gauss_sigma_val = gauss_sigma.getVal(), gauss_sigma_err = gauss_sigma.getError();

// Plot background fit
PlotBackgroundFit(mass, data_incl_run, background, i, 1000); // Using 1000 to indicate inclusive run

// Set up for signal + background fit
mass.setRange("range_full", LeftLowLim_incl[i], RightUpLim_incl[i]);

// Get Crystal Ball parameters from tree
tree->GetEntry(i);

// Define Crystal Ball parameters using inclusive values
RooRealVar mu_cb(Form("mu_cb_%d", i+1), "mu_cb", Mucb_ini_incl[i], Mucb_lowlim_incl[i], Mucb_uplim_incl[i]);
RooRealVar sigma_cb(Form("sigma_cb_%d", i+1), "sigma_cb", Sigmacb_ini_incl[i], Sigmacb_lowlim_incl[i], Sigmacb_uplim_incl[i]);
RooRealVar alphaL_cb(Form("alphaL_cb_%d", i+1), "alphaL_cb", alphaL_ini, alphaL_ini - Nsigma*inc_alphaL, alphaL_ini + Nsigma*inc_alphaL);
RooRealVar nL_cb(Form("nL_cb_%d", i+1), "nL_cb", nL_ini, nL_ini - Nsigma*inc_nL, nL_ini + Nsigma*inc_nL);
RooRealVar alphaR_cb(Form("alphaR_cb_%d", i+1), "alphaR_cb", alphaR_ini, alphaR_ini - Nsigma*inc_alphaR, alphaR_ini + Nsigma*inc_alphaR);
RooRealVar nR_cb(Form("nR_cb_%d", i+1), "nR_cb", nR_ini, nR_ini - Nsigma*inc_nR, nR_ini + Nsigma*inc_nR);



// Definisci la Crystal Ball
RooCrystalBall crystal("crystal", "crystal ball", mass, mu_cb, sigma_cb, alphaL_cb, nL_cb, alphaR_cb, nR_cb);

// Aggiorna i parametri della Gaussiana con i vincoli per il secondo fit
gauss_mu.setVal(gauss_mu_val);
gauss_sigma.setVal(gauss_sigma_val);
gauss_mu.setConstant(true);
gauss_sigma.setRange(gauss_sigma_val - Nsigma*gauss_sigma_err, gauss_sigma_val + Nsigma*gauss_sigma_err);
A.setVal(A_val);
A.setRange(A_val - Nsigma * A_err, A_val + Nsigma * A_err);
B.setVal(B_val);
B.setRange(B_val - Nsigma * B_err, B_val + Nsigma * B_err);
C.setVal(C_val);
C.setRange(C_val - Nsigma * C_err, C_val + Nsigma * C_err);
D.setVal(D_val);
D.setRange(D_val - Nsigma * D_err, D_val + Nsigma * D_err);
E.setVal(E_val);
E.setRange(E_val - Nsigma * E_err, E_val + Nsigma * E_err);

// Create the combined model
RooRealVar frac("frac", "fraction of background", 0.5, 0.0, 1.0);
RooAddPdf model("model", "signal + background", RooArgList(crystal, background), RooArgList(frac));

// Perform signal + background fit
fit_result_incl = model.fitTo(data_incl_run, 
                            RooFit::Range("range_full"),
                            RooFit::MaxCalls(10000000), 
                            RooFit::Save(), 
                            RooFit::SumW2Error(true));

// Plot the fit result
PlotDataFit(mass, data_incl_run, model, background, crystal, 
           LeftLowLim_incl[i], RightUpLim_incl[i], 
           i, 1000, fit_result_incl->status());

// Calculate scale correction for inclusive fit
double mu_data_incl = mu_cb.getVal();
double inc_mu_data_incl = mu_cb.getError();
double scale_postcorr_incl = 1 - (mu_data_incl / mu_mc);
double inc_scale_postcorr_incl = (mu_data_incl / mu_mc) * 
                                sqrt((inc_mu_data_incl / mu_data_incl)*(inc_mu_data_incl / mu_data_incl) + 
                                    (inc_mu_mc / mu_mc)*(inc_mu_mc / mu_mc));

//calculate smearing 
double smearing_postcorr = sigma_cb.getVal() / sigma_ini;
double inc_smearing_postcorr = smearing_postcorr * sqrt((sigma_cb.getError() / sigma_cb.getVal())*(sigma_cb.getError() / sigma_cb.getVal()) + 
                                    (inc_sigma / sigma_ini)*(inc_sigma / sigma_ini));
h_smearing_postcorr->SetBinContent(i+1, smearing_postcorr);
h_smearing_postcorr->SetBinError(i+1, inc_smearing_postcorr);

// You could save this inclusive result to a different histogram if needed

    }    //....oooOO0OOooo........oooOO0OOooo.... FINE LOOP SU PT BIN ....oooOO0OOooo........oooOO0OOooo....

    // After your existing code, add this section to extract h_smearing_inclusiveRun

    // Open the file containing the scale and smearing corrections
    TFile *smearingFile = TFile::Open("scale_corrections.root");
    if (!smearingFile || smearingFile->IsZombie()) {
        std::cerr << "Error: Could not open scale_corrections.root for smearing extraction" << std::endl;
        return;
    }

    // Extract the h_smearing_inclusiveRun histogram
    TH1D *h_smearing_inclusiveRun = (TH1D*)smearingFile->Get("h_smearing_inclusiveRun");
    if (!h_smearing_inclusiveRun) {
        std::cerr << "Error: Could not find h_smearing_inclusiveRun histogram in scale_corrections.root" << std::endl;
        smearingFile->Close();
        return;
    }

    // Create a canvas to display the smearing histogram
    TCanvas *cSmearing = new TCanvas("cSmearing", "", 900, 600);
    gStyle->SetOptStat(0); // Turn off statistics box

    // Set up the appearance for original smearing histogram
    h_smearing_inclusiveRun->SetTitle("");
    h_smearing_inclusiveRun->GetXaxis()->SetTitle("p_{T} [GeV]");
    h_smearing_inclusiveRun->GetYaxis()->SetTitle("#sigma_{data}/#sigma_{MC}");
    h_smearing_inclusiveRun->SetMarkerStyle(20);
    h_smearing_inclusiveRun->SetMarkerColor(violaCMS);
    h_smearing_inclusiveRun->SetLineColor(violaCMS);

    // Explicitly set the x-axis range to extend to 42 GeV
    h_smearing_inclusiveRun->GetXaxis()->SetRangeUser(0, 42);
    h_smearing_inclusiveRun->GetYaxis()->SetRangeUser(0.8, 2.6);

    // Set up the appearance for post-correction smearing histogram
    h_smearing_postcorr->SetMarkerStyle(21);
    h_smearing_postcorr->SetLineColor(arancioneCMS); // Changed to gialloCMS as requested earlier
    h_smearing_postcorr->SetMarkerColor(arancioneCMS);

    // Draw the histograms
    h_smearing_inclusiveRun->Draw("PE");
    h_smearing_postcorr->Draw("PE SAME");

    // Add a horizontal line at y=1 for reference - use the extended x-range
    TLine *smearingLine = new TLine(4, 1, 42, 1);
    smearingLine->SetLineColor(bluCMS);
    smearingLine->SetLineWidth(2);
    smearingLine->SetLineStyle(kDashed);
    smearingLine->Draw();

    // Create and draw the legend
    TLegend *legend = new TLegend(0.65, 0.75, 0.89, 0.89);
    legend->AddEntry(h_smearing_inclusiveRun, "Before Correction", "lp");
    legend->AddEntry(h_smearing_postcorr, "After Correction", "lp");
    legend->AddEntry(smearingLine, "Perfect agreement", "l");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();

    // Add CMS label in bold
    TLatex cmsLabel;
    cmsLabel.SetNDC();
    cmsLabel.SetTextFont(62);  // Bold font
    cmsLabel.SetTextSize(0.05);
    cmsLabel.SetTextAlign(11);
    cmsLabel.DrawLatex(0.10, 0.92, "CMS");

    // Add Preliminary in italics (not bold)
    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextFont(52);  // Italic font
    prelimLabel.SetTextSize(0.05);
    prelimLabel.SetTextAlign(11);
    prelimLabel.DrawLatex(0.18, 0.92, "Preliminary");

    // Draw lumi/energy information
    TLatex lumiLabel;
    lumiLabel.SetNDC();
    lumiLabel.SetTextFont(42);  // Regular font
    lumiLabel.SetTextSize(0.045);
    lumiLabel.SetTextAlign(31);
    lumiLabel.DrawLatex(0.90, 0.92, "38.01 fb^{-1} (2022, 13.7 TeV)");

    // Save the canvas
    cSmearing->SaveAs("PlotConID2022/FitCorrectedDataDiag_id2022/Smearing_Comparison_vs_Pt.png");
    delete cSmearing;

    // Close the file
    smearingFile->Close();
}