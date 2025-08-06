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
#include <TColor.h>

#define Nbins 6

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....Funzione globale per scritta MC....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void WriteSimulation(){
// Aggiungi una scritta in testo italic al bordo esterno della cornice
    TLatex latex;
    latex.SetTextSize(0.05); // Dimensione del testo
    latex.SetTextFont(42);   // Tipo di font
    latex.SetTextAlign(11);  // Allineamento: 11 = sinistra, alto
    latex.SetTextAngle(0);   // Angolo del testo
    latex.DrawLatexNDC(0.12, 0.93, "#it{Simulation}"); // Posizione e testo

}

void FitProjectionMCAllFullEle() {
    // Create output directories if they don't exist
    gSystem->Exec("mkdir -p PlotConID2022/FitMCAllFullEle_id2022");
    
    ///////////////Definisco i colori
    int bluCMS = TColor::GetColor("#5790FC");
    int rosaCMS = TColor::GetColor("#964A8B");
    int rossoCMS = TColor::GetColor("#A42536");
    int gialloCMS = TColor::GetColor("#F89C20");
    ///////////////
    // Apri il file ROOT contenente gli istogrammi
    TFile *file = TFile::Open("outputHistograms_MC.root");

    // Array per i parametri dei fit
    double means[Nbins], sigmas[Nbins], alphaLs[Nbins], nLs[Nbins], alphaRs[Nbins], nRs[Nbins], meanInf[Nbins], meanSup[Nbins], mean_histo[Nbins];
    double sigma_over_mu[Nbins], inc_sigma_over_mu[Nbins];

     // Initialize MC fit parameters (customizable)
    for (int i = 0; i < Nbins; i++) {
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
    double massMin[Nbins] = {2.6, 2.2, 2.6, 2.6, 2.5, 2.65};
    double massMax[Nbins] = {3.3, 3.9, 3.4, 3.6, 3.7, 3.5};
    int rebin_factor[Nbins] = {1, 1, 1, 1, 1, 1};

    double Bin_centers[Nbins] = {5.75, 8.25, 10.5, 12.5, 17, 30}; //Centri dei bin in P_t
    double Bin_halfwidth[Nbins] = {1.75, 0.75, 1.5, 1.5, 3, 10}; //Larghezze dei bin in P_t

    double mu[Nbins], inc_mu[Nbins], val_sigma[Nbins], inc_sigma[Nbins];
    double val_nL[Nbins], inc_nL[Nbins], val_alphL[Nbins], inc_alphL[Nbins], val_nR[Nbins], inc_nR[Nbins], val_alphR[Nbins], inc_alphR[Nbins]; //array per salvare i parametri

    // Prendi l'istogramma 2D h_invmass_vs_pt_all_fullele
    TH2D* h2d = (TH2D*)file->Get("h_invmass_vs_pt_all_fullele");
    if (!h2d) {
        std::cerr << "Istogramma 2D h_invmass_vs_pt_all_fullele non trovato!" << std::endl;
        return;
    }

    // Loop su tutti i bin di pt
    for (int i = 0; i < Nbins; i++) {
        // Proietta la Y (massa) per il bin di pt i usando il bin i+1
        TH1D* hist = h2d->ProjectionY(Form("proj_bin_%d", i+1), i+1, i+1);
        if (!hist || hist->GetEntries() == 0) {
            std::cerr << "Proiezione Y vuota per il bin pt " << i+1 << std::endl;
            continue;
        }
        hist->Rebin(rebin_factor[i]);

        // Ottieni media e RMS per l'inizializzazione
        means[i] = hist->GetMean();
        sigmas[i] = hist->GetRMS() / 2.0;  // Utilizza RMS/2 come stima iniziale per sigma

        // Definisci la variabile di massa invariante
        RooRealVar mass("mass", "m(e^{+}e^{-})", massMin[i], massMax[i]); 

        // Converto l'istogramma in un RooDataHist
        RooDataHist data("data", "Dataset from histogram", mass, hist);

        // Definisci i parametri della Crystal Ball asimmetrica
        RooRealVar cb_mean(Form("cb_mean_%d", i+1), "Mean of CB", means[i], meanInf[i], meanSup[i]);
        RooRealVar cb_sigma(Form("cb_sigma_%d", i+1), "Sigma of CB", sigmas[i], 0.1 * sigmas[i], 2 * sigmas[i]);
        RooRealVar cb_alphaL(Form("cb_alphaL_%d", i+1), "AlphaL of CB", alphaLs[i], 0.0, 15.0);
        RooRealVar cb_nL(Form("cb_nL_%d", i+1), "nL of CB", nLs[i], 0.0, 500.0);
        RooRealVar cb_alphaR(Form("cb_alphaR_%d", i+1), "AlphaR of CB", alphaRs[i], 0.0, 15.0);
        RooRealVar cb_nR(Form("cb_nR_%d", i+1), "nR of CB", nRs[i], 0.0, 500.0);

        // Crea la Crystal Ball asimmetrica 
        RooCrystalBall cb(Form("cb_%d", i+1), "Asymmetric Crystal Ball", mass, cb_mean, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR);

        // Esegui il fit
        RooFitResult *fit_result_simplex = cb.fitTo(data, RooFit::Minimizer("Minuit2", "Simplex"), RooFit::MaxCalls(1000000), RooFit::Save(),  RooFit::SumW2Error(true), RooFit::PrintLevel(-1));

        RooFitResult *fit_result = cb.fitTo(data, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::MaxCalls(10000000), RooFit::Save(), RooFit::SumW2Error(true), RooFit::PrintLevel(-1));
        ////////////////////////////////////////////////////////
        TCanvas *canvas = new TCanvas(Form("canvas_%d", i+1), Form("Fit for proj_bin_%d", i+1), 900, 700);

        // Plotta l'istogramma e il fit
        RooPlot *frame = mass.frame();
        frame->SetTitle("");
        data.plotOn(frame);
        cb.plotOn(frame, RooFit::LineColor(bluCMS), RooFit::LineWidth(5));
        frame->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");

        double chi2 = frame->chiSquare();
        // Stampare il chi-quadro sul plot
        TPaveText *paveText = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
        paveText->AddText(Form("#chi^{2} = %.2f", chi2));
        paveText->AddText(Form("#mu = %.4f +/- %.4f", cb_mean.getVal(), cb_mean.getError()));
        paveText->AddText(Form("#sigma = %.4f +/- %.4f", cb_sigma.getVal(), cb_sigma.getError()));
        paveText->SetFillColor(0);
        frame->addObject(paveText);

        frame->Draw();
        WriteSimulation();

        // Salva la canvas 
        canvas->SaveAs(Form("PlotConID2022/FitMCAllFullEle_id2022/fit_proj_bin_%d.png", i+1));
        // deallocazioni
        delete canvas;
        delete frame;

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
        //Salvo mu, errore su mu, sigma ed errore su sigma 
        mu[i] = cb_mean.getVal();
        inc_mu[i] = cb_mean.getError();
        val_sigma[i] = cb_sigma.getVal();
        inc_sigma[i] = cb_sigma.getError();
        mean_histo[i] = hist->GetMean();
        //salvo anche i parametri delle code
        val_nL[i] = cb_nL.getVal();
        inc_nL[i] = cb_nL.getError();
        val_alphL[i] = cb_alphaL.getVal();
        inc_alphL[i] = cb_alphaL.getError();

        val_nR[i] = cb_nR.getVal();
        inc_nR[i] = cb_nR.getError();
        val_alphR[i] = cb_alphaR.getVal();
        inc_alphR[i] = cb_alphaR.getError();
        //errore relativo
        sigma_over_mu[i] = cb_sigma.getVal()/cb_mean.getVal();
        inc_sigma_over_mu[i] = cb_sigma.getError()/cb_mean.getVal();
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                          FIT INCLUSIVO IN PT                                           

    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Salvo i parametri per inizializzare il fit dei dati
        // Crea un file ROOT
        TFile *file_param = TFile::Open("fit_results_all_fullele.root", "RECREATE");

        // Crea un albero ROOT per memorizzare i dati
        TTree *tree = new TTree("fitResults", "Fit Results");

        // Definisci variabili per i dati
        int bin;
        double mean, meanError, sigma, sigmaError, nL, nLError, alphaL, alphaLError, nR, nRError, alphaR, alphaRError;

        // Aggiungi le variabili al TTree
        tree->Branch("Bin", &bin);
        tree->Branch("Mean", &mean);
        tree->Branch("MeanError", &meanError);
        tree->Branch("Sigma", &sigma);
        tree->Branch("SigmaError", &sigmaError);
        tree->Branch("nL", &nL);
        tree->Branch("nLError", &nLError);
        tree->Branch("alphaL", &alphaL);
        tree->Branch("alphaLError", &alphaLError);
        tree->Branch("nR", &nR);
        tree->Branch("nRError", &nRError);
        tree->Branch("alphaR", &alphaR);
        tree->Branch("alphaRError", &alphaRError);

        // Riempie l'albero con i dati
        for (int i = 0; i < Nbins; i++) {
         bin = i + 1;
         mean = mu[i];
         meanError = inc_mu[i];
         sigma = val_sigma[i];
         sigmaError = inc_sigma[i];
         nL = val_nL[i];
         nLError = inc_nL[i];
         alphaL = val_alphL[i];
         alphaLError = inc_alphL[i];
         nR = val_nR[i];
         nRError = inc_nR[i];
         alphaR = val_alphR[i];
         alphaRError = inc_alphR[i];
         tree->Fill();
        }

        // Salva e chiudi il file
        file_param->Write();
        file_param->Close();
        ///////////////////////////////////////



    //Plot mu vs Pt

    TGraphErrors *gr = new TGraphErrors(Nbins, Bin_centers, mu, Bin_halfwidth, inc_mu);
    gr->SetTitle("Peak center vs Pt;P_{t} [GeV]; Inv Mass [GeV]");

    //TGraphErrors *gr2 = new TGraphErrors(Nbins, Bin_centers, val_sigma, Bin_halfwidth, inc_sigma);
    TGraphErrors *gr2 = new TGraphErrors(Nbins, Bin_centers, sigma_over_mu, Bin_halfwidth, inc_sigma_over_mu);
    //gr2->SetTitle("Peak width vs Pt;P_{t} [GeV]; #sigma [GeV]");
    gr2->SetTitle("Peak width vs Pt;P_{t} [GeV]; #sigma/#mu");

    // Stile del grafico
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);

    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColor(kMagenta);
    gr2->SetLineColor(kMagenta);

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    gr->Draw("AP");

    c1->cd(2);
    gr2->Draw("AP");

    file->Close();
}

