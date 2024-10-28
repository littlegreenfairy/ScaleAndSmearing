#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooCrystalBall.h>
#include <RooFitResult.h>
#include <iostream>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TF1.h>
#include <RooPolynomial.h>
#include <TColor.h>

#define NbinsPt 6
#define NbinsRun 5

void Jpsi_binPt_and_runN() {
    // Apri il file ROOT
    TFile *file_data = TFile::Open("outputHistograms_DATA_partF.root");
    if (!file_data || file_data->IsZombie()) {
        std::cout << "Errore: Impossibile aprire il file ROOT!" << std::endl;
        return;
    }

    // Array di canvas, uno per ogni bin di pT
    TCanvas *canvas[NbinsPt];
    for (int ptBin = 1; ptBin <= NbinsPt; ++ptBin) {
        // Crea un canvas per ogni bin di pT
        canvas[ptBin-1] = new TCanvas(Form("canvas_ptBin_%d", ptBin), Form("Bin pT %d", ptBin), 2500, 1000);
        canvas[ptBin-1]->Divide(3,2);  // Layout 2x3 per mostrare i 6 istogrammi
        
        for (int runBin = 1; runBin <= NbinsRun; ++runBin) {
            // Costruisci il nome dell'istogramma
            TString histName = Form("proj_bins_%d_%d", ptBin, runBin);
            TH1D *hist = (TH1D*)file_data->Get(histName);
            
            if (hist) {
                // Posiziona l'istogramma nel pad corretto
                canvas[ptBin-1]->cd(runBin);
                hist->SetLineWidth(4);
                hist->SetLineColor(kBlue);
                hist->Draw();
            } else {
                std::cout << "Istogramma " << histName << " non trovato!" << std::endl;
            }
        }
        
        // Visualizza il canvas
        canvas[ptBin-1]->Update();
        canvas[ptBin-1]->SaveAs(Form("Histos_bin_runNumber/canvas_ptBin_%d.png", ptBin));
        delete canvas[ptBin-1];
    }

    // Chiudi il file alla fine
    file_data->Close();
}
