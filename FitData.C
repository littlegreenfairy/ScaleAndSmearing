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
    Color_t gialloCMS = TColor::GetColor("#ffcc00");
    Color_t violaCMS = TColor::GetColor("#660099");
    Color_t rossoCMS = TColor::GetColor("#cc0000");

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
    cBackgroundFit->SaveAs(Form("PlotConID2022/FitData_id2022/Pt_bin%d/BackgroundFit_Pt%d_Run%d.png", i+1, i+1, j+1));
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
    latex.DrawLatex(0.15, 0.66, Form("#chi^{2} = %.2f", chi2));

    // Save plot as image (optional)
    c->SaveAs(Form("PlotConID2022/FitData_id2022/Pt_bin%d/DataFit_Pt%d_Run%d.png", i+1, i+1, j+1));
    delete c;
}


void FitData(){
///////////////
    // Apri il file ROOT contenente gli istogrammi
    TFile *file = TFile::Open("outputHistograms_DATA_partF.root");

    // Definisci i limiti personalizzati e fattore di rebinning
    double LeftLowLim[NbinsPt][NbinsRun] = {
        {1.6, 1.4, 1.4, 1.6, 1.7, 1.6, 1.5, 1.6, 1.5}, 
        {1.6, 1.4, 1.4, 1.6, 1.7, 1.4, 1.5, 1.6, 1.7}, 
        {1.6, 1.4, 1.4, 1.6, 1.7, 1.6, 1.5, 1.6, 1.7}, 
        {1.6, 1.4, 1.4, 1.6, 1.7, 1.6, 1.5, 1.6, 1.7},
        {1.6, 1.9, 1.9, 1.6, 1.7, 1.6, 1.6, 1.6, 1.7}, 
        {2.45, 2.45, 2.45, 2, 2.45, 2, 1.5, 2, 2.45}
    };  // Limiti sinistri personalizzati
    
    double LeftUpLim[NbinsPt][NbinsRun] = {
        {2.4, 2.5, 2.2, 2.2, 2.5, 2.4, 2.3, 2.4, 2.4}, 
        {2.4, 2.5, 2.2, 2.2, 2.5, 2.5, 2.3, 2.4, 2.5}, 
        {2.4, 2.5, 2.2, 2.2, 2.5, 2.4, 2.3, 2.4, 2.5}, 
        {2.4, 2.5, 2.2, 2.2, 2.5, 2.4, 2.3, 2.4, 2.5},
        {2.4, 2.5, 2.2, 2.2, 2.6, 2.4, 2.6, 2.4, 2.5}, 
        {2.8, 2.8, 2.8, 2.2, 2.8, 2.4, 2.3, 2.4, 2.8}
    };
    
    double RightLowLim[NbinsPt][NbinsRun] = {
        {3.6, 3.4, 3.6, 3.6, 3.6, 3.5, 3.6, 3.5, 3.5}, 
        {3.6, 3.4, 3.6, 3.6, 3.6, 3.5, 3.6, 3.5, 3.6}, 
        {3.6, 3.4, 3.6, 3.6, 3.6, 3.5, 3.6, 3.5, 3.6}, 
        {3.6, 3.4, 3.6, 3.6, 3.6, 3.5, 3.6, 3.5, 3.6}, 
        {3.6, 3.4, 3.6, 3.6, 3.6, 3.5, 3.6, 3.5, 3.6}, 
        {3.6, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6}
    }; // Limiti destri personalizzati
    
    double RightUpLim[NbinsPt][NbinsRun] = {
        {5, 5, 5, 5, 5, 5, 5, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5.2, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5, 5, 5}, 
        {5, 5, 5, 5, 5, 5, 5, 5, 5}
    };

    //limiti inclusivi in run number (aggiungere)
    double LeftLowLim_incl[NbinsPt] = {1.15, 1.3, 1.6, 1.6, 1.65, 2};
    double LeftUpLim_incl[NbinsPt] = {2.4, 2.5, 2.45, 2.4, 2.65, 2.7};
    double RightLowLim_incl[NbinsPt] = {3.65, 3.55, 3.55, 3.55, 3.55, 3.55};
    double RightUpLim_incl[NbinsPt] = {5.2, 5.2, 4.8, 5, 5.2, 5.2};
    

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
        {3.0442, 3.0337, 3.0422, 2.9975, 3.0782, 3.0550, 3.0600, 3.0650, 3.0700},
        {3.0442, 3.0337, 3.0943, 3.0489, 3.0782, 3.0600, 3.0650, 3.0700, 3.0750},
        {3.0442, 3.0337, 3.0422, 3.0489, 3.0782, 3.0650, 3.0700, 3.0750, 3.0800},
        {3.0442, 3.0337, 3.0422, 3.0489, 3.1098, 3.0700, 3.0750, 3.0800, 3.0850},
        {3.0442, 3.0337, 3.0422, 3.0489, 3.0782, 3.0750, 3.0800, 3.0850, 3.0900}
    };
    
    double Mucb_lowlim[NbinsPt][NbinsRun] = {
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.9, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8},
        {2.8, 2.7, 2.5, 2.5, 2.8, 2.7, 2.6, 2.7, 2.8}
    };
    
    double Mucb_uplim[NbinsPt][NbinsRun] = {
        {3.2, 3.5, 3.5, 3.5, 3.2, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.5, 3.2, 3.5, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.2, 3.5, 3.5, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.5, 3.5, 3.5, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.5, 3.5, 3.5, 3.3, 3.4, 3.3, 3.2},
        {3.2, 3.5, 3.5, 3.5, 3.5, 3.3, 3.4, 3.3, 3.2}
    };
    
    double Sigmacb_ini[NbinsPt][NbinsRun] = {
        {0.1637, 0.1503, 0.1351, 0.1271, 0.15, 0.14, 0.13, 0.14, 0.15},
        {0.1637, 0.1503, 0.1351, 0.1271, 0.1115, 0.1342, 0.12, 0.13, 0.1452},
        {0.1637, 0.1503, 0.1355, 0.1271, 0.1115, 0.12, 0.11, 0.12, 0.13},
        {0.1637, 0.1503, 0.1351, 0.1271, 0.1115, 0.11, 0.10, 0.11, 0.12},
        {0.1637, 0.1503, 0.1351, 0.1271, 0.1115, 0.11, 0.10, 0.11, 0.12},
        {0.1637, 0.1503, 0.1351, 0.1271, 0.1115, 0.10, 0.09, 0.10, 0.11}
    };
    
    double Sigmacb_uplim[NbinsPt][NbinsRun] = {
        {0.3, 0.5, 0.5, 0.5, 0.5, 0.4, 0.4, 0.4, 0.4},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.2, 0.3, 0.3, 0.3},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.3, 0.3, 0.3, 0.3},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.3, 0.3, 0.3, 0.3},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.3, 0.3, 0.3, 0.3},
        {0.3, 0.5, 0.5, 0.5, 0.15, 0.3, 0.3, 0.3, 0.3}
    };
    
    double Sigmacb_lowlim[NbinsPt][NbinsRun] = {
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
        {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}
    };

    //parametri inclusivi in pt
    double Mucb_ini_incl[NbinsPt] = {3.0449, 3.0421, 3.0343, 3.0714, 3.0839, 3.1217};
    double Mucb_lowlim_incl[NbinsPt] = {2.9, 2.9, 2.5, 2.5, 2.5, 3};
    double Mucb_uplim_incl[NbinsPt] = {3.2, 3.2, 3.2, 3.2, 3.2, 3.2};
    double Sigmacb_ini_incl[NbinsPt] = {0.1613, 0.149, 0.1440, 0.1203, 0.11, 0.09};
    double Sigmacb_uplim_incl[NbinsPt] = {0.3, 0.3, 0.3, 0.3, 0.2, 0.2};
    double Sigmacb_lowlim_incl[NbinsPt] = {0.05, 0.05, 0.05, 0.05, 0, 0.05};


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
        {3.5, 3.6, 3.6, 3.2, 3.5, 3.4, 3.45, 3.4, 3.3},
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
        {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
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
    double gauss_mu_init_incl[NbinsPt] = {3.55, 3.55, 3.5, 3.62, 3.6365, 3.6195};
    double gauss_mu_low_incl[NbinsPt] = {3.5, 3.4, 3.4, 3.5, 3.55, 3.5};
    double gauss_mu_up_incl[NbinsPt] = {3.7, 3.7, 3.7, 3.7, 3.7, 3.7};

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
    TH2D *h_scale = new TH2D("h_scale", "Scale between data and MC ; Pt ; Run number", NbinsPt, Ptbins, NbinsRun, runBins);
    TH2D *h_corr_1ele = new TH2D("h_corr_1ele", "Scale correction for single electron energy; Pt ; Run number", NbinsPt, Ptbins, NbinsRun, runBins);
    TH1D *h_corr_1ele_inclusiveRun = new TH1D("h_corr_1ele_inclusiveRun", "Scale corrections for single electron energy; Pt; correction", NbinsPt, Ptbins);
    TH1D *h_scale_inclusiveRun = new TH1D("h_scale_inclusiveRun", "Scale between data and MC; Pt; correction", NbinsPt, Ptbins);
    TH1D *h_smearing_inclusiveRun = new TH1D("h_smearing_inclusiveRun", "Smearing corrections for single electron energy; Pt; correction", NbinsPt, Ptbins);
    TGraphErrors *graph_scale = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_smearing = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_meandata = new TGraphErrors(NbinsPt);
    graph_meandata->SetName("Graph_MeanData");

    //File .root per salvare gli istogrammi con le correzioni
    TFile *file_corrections = new TFile("scale_corrections.root", "RECREATE");
    //....oooOO0OOooo........oooOO0OOooo.... Loop su tutti gli istogrammi ....oooOO0OOooo........oooOO0OOooo....

    double compatibility_psi2s[NbinsPt];

    for(int i=0; i < NbinsPt; i++){

        for(int j=0; j < NbinsRun; j++){

            // Ottieni l'istogramma
            TH1D *hist = (TH1D*)file->Get(Form("proj_bins_%d_%d", i+1, j+1));
            if (!hist) {
                std::cerr << "Istogramma non trovato: proj_bins_" << i + 1 << "_" << j + 1 << std::endl;
                continue;
            }
            //hist->Rebin(rebin_factor[i][j]);

            //....oooOO0OOooo.... fit fondo ....oooOO0OOooo....
            // Definisci la variabile di massa invariante
            RooRealVar mass("mass", "m(e^{+}e^{-})", 0, 6); 

            // Converto l'istogramma in un RooDataHist
            RooDataHist data("data", "Dataset from histogram", mass, hist);

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

            //Calcolo la scala e la correzione per il singolo elettrone
            double mu_data = mu_cb.getVal();
            double inc_mu_data = mu_cb.getError();
            double scale = 1 - (mu_data / mu_mc);
            double inc_scale = (mu_data / mu_mc) * sqrt((inc_mu_data / mu_data)*(inc_mu_data / mu_data) + (inc_mu_mc / mu_mc)*(inc_mu_mc / mu_mc));
            double corr_1ele = mu_mc / mu_data;
            double inc_corr_1ele = corr_1ele * sqrt((inc_mu_data / mu_data)*(inc_mu_data / mu_data) + (inc_mu_mc / mu_mc)*(inc_mu_mc / mu_mc));
            //Scrivo i parametri di interesse negli istogrammi 2D
            h_scale->SetBinContent(i+1, j+1, scale);
            h_scale->SetBinError(i+1, j+1, inc_scale);
            h_corr_1ele->SetBinContent(i+1, j+1, corr_1ele);
            h_corr_1ele->SetBinError(i+1, j+1, inc_corr_1ele);

        } //FINE LOOP SU RUN NUMBER


        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
        //                                      Aggiungo un fit inclusivo in Run Number
        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
        TH1D *hist_inclusiveRun = (TH1D*)file->Get(Form("proj_bin_%d", i+1));
            if (!hist_inclusiveRun) {
                std::cerr << "Istogramma non trovato: proj_bin_" << i + 1 << std::endl;
                continue;
            }
        //Ora fitto questo istogramma con lo stesso modello di prima ma senza il binning su j
        RooRealVar mass("mass", "m(e^{+}e^{-})", 0, 6); 

        // Converto l'istogramma in un RooDataHist
        RooDataHist data_inclusiveRun("data_inclusiveRun", "Dataset from histogram", mass, hist_inclusiveRun);

        // Definisci i parametri del polinomio di secondo grado
            RooRealVar A(Form("A_%d", i+1), "4th deg coeff", 0,  -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar B(Form("B_%d", i+1), "3rd deg coeff", 0,  -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar C(Form("C_%d", i+1), "2nd deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar D(Form("D_%d", i+1), "1st deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
            RooRealVar E(Form("E_%d", i+1), "0 deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
            RooPolynomial poly("poly", "Polynomial of 4th degree", mass, RooArgList(A, B, C, D, E));

            // Definisci i parametri della Gaussiana
            RooRealVar gauss_mu(Form("gauss_mu_%d", i+1), "Gaussian mean", gauss_mu_init_incl[i], gauss_mu_low_incl[i], gauss_mu_up_incl[i]);
            RooRealVar gauss_sigma(Form("gauss_sigma_%d", i+1), "Gaussian sigma", gauss_sigma_init[i][0], gauss_sigma_low[i][0], gauss_sigma_up[i][0]);
            RooGaussian gauss("gauss", "Gaussian component", mass, gauss_mu, gauss_sigma);

            // Definisci il modello di fondo combinando il polinomio con la Gaussiana
            RooRealVar frac_gauss("frac_gauss", "fraction of Gaussian", 0.3, 0.0, 1.0);
            RooAddPdf background("background", "Background Model", RooArgList(poly, gauss), RooArgList(frac_gauss));

            //Definisci il range per il fit del fondo
            mass.setRange("range1", LeftLowLim_incl[i], LeftUpLim_incl[i]);
            mass.setRange("range2", RightLowLim_incl[i], RightUpLim_incl[i]);

            //////////////////// Esegui fit sul solo fondo //////////////////////
            RooFitResult *fit_result_incl = background.fitTo(data_inclusiveRun, RooFit::Range("range1,range2"),RooFit::MaxCalls(10000000), RooFit::Save(), RooFit::SumW2Error(true));

            // Ottieni i valori dei parametri e i loro errori dal primo fit
            double A_val = A.getVal(), A_err = A.getError();
            double B_val = B.getVal(), B_err = B.getError();
            double C_val = C.getVal(), C_err = C.getError();
            double D_val = D.getVal(), D_err = D.getError();
            double E_val = E.getVal(), E_err = E.getError();
            double gauss_mu_val = gauss_mu.getVal(), gauss_mu_err = gauss_mu.getError();
            double gauss_sigma_val = gauss_sigma.getVal(), gauss_sigma_err = gauss_sigma.getError();
            //////////////////////////////////////////////////////// Plot solo fondo
            PlotBackgroundFit(mass, data_inclusiveRun, background, i, 999); //999 indica che è inclusivo in pt (volevo riciclare la funzione)

            //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
            //                                          FIT SEGNALE + FONDO
            //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

            mass.setRange("range_full", LeftLowLim[i][0], RightUpLim[i][0]);

            // Leggi i parametri della Crystal Ball dal TTree per il bin corrente
            tree->GetEntry(i);

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
            fit_result_incl = model.fitTo(data_inclusiveRun, RooFit::Range("range_full"),RooFit::MaxCalls(10000000), RooFit::Save(), RooFit::SumW2Error(true));
            //Plot del fit completo
            PlotDataFit(mass, data_inclusiveRun, model, background, crystal, LeftLowLim[i][0], RightUpLim[i][0], i, 999, fit_result_incl->status());
        
        //calcolo la correzione inclusiva in pt
            double mu_data = mu_cb.getVal();
            double inc_mu_data = mu_cb.getError();
            double scale = 1 - (mu_data / mu_mc);
            double inc_scale = (mu_data / mu_mc) * sqrt((inc_mu_data / mu_data)*(inc_mu_data / mu_data) + (inc_mu_mc / mu_mc)*(inc_mu_mc / mu_mc));
            double corr_1ele = mu_mc / mu_data;
            double inc_corr_1ele = corr_1ele * sqrt((inc_mu_data / mu_data)*(inc_mu_data / mu_data) + (inc_mu_mc / mu_mc)*(inc_mu_mc / mu_mc));
            
            //Scrivo i parametri di interesse negli istogrammi 2D
            h_scale_inclusiveRun->SetBinContent(i+1, scale);
            h_scale_inclusiveRun->SetBinError(i+1, inc_scale);
            h_smearing_inclusiveRun->SetBinContent(i+1, sigma_cb.getVal()/sigma_ini);
            h_smearing_inclusiveRun->SetBinError(i+1, (sigma_cb.getVal()/sigma_ini)*sqrt((inc_sigma/sigma_ini)*(inc_sigma/sigma_ini) + (sigma_cb.getError()/sigma_cb.getVal())*(sigma_cb.getError()/sigma_cb.getVal())));
            h_corr_1ele_inclusiveRun->SetBinContent(i+1, corr_1ele);
            h_corr_1ele_inclusiveRun->SetBinError(i+1, inc_corr_1ele);

            //verifico se la differenza tra il centro della gaussiana e il centro della crystal ball è compatibile con 0.5
            compatibility_psi2s[i] = (gauss_mu_val - mu_data - 0.5) / sqrt(gauss_mu_err*gauss_mu_err + inc_mu_data*inc_mu_data);

            graph_meandata->SetPoint(i, bincenters[i], mu_data);
            graph_meandata->SetPointError(i, binhalfwidths[i], inc_mu_data);
            //calcolo scala e smearing
            graph_scale->SetPoint(i, bincenters[i], 1 - mu_data/mu_mc);
            graph_scale->SetPointError(i, binhalfwidths[i], (mu_data/mu_mc)*sqrt((inc_mu_data/mu_data)*(inc_mu_data/mu_data) + (inc_mu_mc/mu_mc)*(inc_mu_mc/mu_mc)));

            graph_smearing->SetPoint(i, bincenters[i], sigma_cb.getVal()/sigma_ini);
            graph_smearing->SetPointError(i, binhalfwidths[i], (sigma_cb.getVal()/sigma_ini)*sqrt((sigma_cb.getError()/sigma_cb.getVal())*(sigma_cb.getError()/sigma_cb.getVal()) + (inc_sigma/sigma_ini)*(inc_sigma/sigma_ini)));
    }

    //plot correzioni vs Pt
    TCanvas *c_corr = new TCanvas("c_corr", "Single electron corrections vs p_{T}", 800, 800);
    c_corr->SetLeftMargin(0.15);
    h_corr_1ele_inclusiveRun->SetLineColor(violaCMS);
    h_corr_1ele_inclusiveRun->GetXaxis()->SetTitle("p_{T} [GeV]");
    h_corr_1ele_inclusiveRun->GetYaxis()->SetTitle("m^{(J/#psi)}_{MC} / m^{(J/#psi)}_{data}");
    h_corr_1ele_inclusiveRun->SetMarkerStyle(21);
    h_corr_1ele_inclusiveRun->SetMarkerColor(violaCMS);
    h_corr_1ele_inclusiveRun->SetStats(kFALSE);
    h_corr_1ele_inclusiveRun->Draw("E1");

    c_corr->SaveAs("PlotConID2022/Single_ele_corrections_vspT.png");
    delete c_corr;

    graph_scale->SetTitle("");
    graph_smearing->SetTitle("");
    //plot scala e smearing vs p_t
    TCanvas *c_scale_smear = new TCanvas("c_scale_smear", "Scale and spearing vs p_{T}", 1600, 600);
    c_scale_smear->Divide(2,1);
    c_scale_smear->cd(1);
    gPad->SetLeftMargin(0.15);
    graph_scale->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_scale->GetYaxis()->SetTitle("1 - m^{(J/#psi)}_{data}/m^{(J/#psi)}_{MC}");
    graph_scale->SetMarkerStyle(21);
    graph_scale->SetLineColor(rossoCMS);
    graph_scale->SetMarkerColor(rossoCMS);
    graph_scale->Draw("APE");
    c_scale_smear->cd(2);
    
    graph_smearing->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_smearing->GetYaxis()->SetTitle("#sigma_{data}/#sigma_{MC}");
    graph_smearing->SetMarkerStyle(21);
    graph_smearing->SetLineColor(violaCMS);
    graph_smearing->SetMarkerColor(violaCMS);
    graph_smearing->Draw("APE");

    c_scale_smear->cd(1);
    TLatex cmsLabel;
    cmsLabel.SetNDC();                          // use normalized coordinates
    cmsLabel.SetTextFont(62);                   // Helvetica Bold
    cmsLabel.SetTextSize(0.05);                 // adjust to match CMS style
    cmsLabel.SetTextAlign(11);                  // align top-left corner of text
    cmsLabel.DrawLatex(0.15, 0.92, "CMS");

    // Add Preliminary in italics (not bold)
    TLatex prelimLabel;
    prelimLabel.SetNDC();
    prelimLabel.SetTextFont(52);   // Italic font
    prelimLabel.SetTextSize(0.05);
    prelimLabel.SetTextAlign(11);
    prelimLabel.DrawLatex(0.24, 0.92, "Preliminary");

    // Draw lumi/energy on the first pad
    TLatex lumiLabel;
    lumiLabel.SetNDC();
    lumiLabel.SetTextFont(42);                  // regular font for the detail
    lumiLabel.SetTextSize(0.045);
    lumiLabel.SetTextAlign(31);                 // align top-right corner of text
    lumiLabel.DrawLatex(0.88, 0.92, "38.01 fb^{-1} (2022, 13.7 TeV)");

    // (Optional) replicate on the second pad if you want them on both panels:
    c_scale_smear->cd(2);
    cmsLabel.DrawLatex(0.10, 0.92, "CMS");
    prelimLabel.DrawLatex(0.19, 0.92, "Preliminary");
    lumiLabel.DrawLatex(0.9, 0.92, "38.01 fb^{-1} (2022, 13.7 TeV)");

    c_scale_smear->SaveAs("PlotConID2022/ScaleAndSmearing_vsPt.png");
    delete c_scale_smear;

//stampo compatibilità con psi2s
for(int i=0; i< NbinsPt; i++){
    std::cout << "compatibilità tra Deltamu e 0.5: Bin " << i+1 << "  valore:" << compatibility_psi2s[i] << std::endl;
}

file_corrections->cd();
h_scale->Write();
h_corr_1ele->Write();
h_scale_inclusiveRun->Write();
h_smearing_inclusiveRun->Write();
h_corr_1ele_inclusiveRun->Write();
graph_meandata->Write();

file_corrections->Close();
delete file_corrections;
}