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

void FitProjectionMC() {
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

    // Inizializza i parametri provvisoriamente
    for (int i = 0; i < Nbins; i++) {
        means[i] = 3;  // Inizializzazione provvisoria, sarà impostato a partire dalla media dell'istogramma
        sigmas[i] = 0.14;  // Stima iniziale per sigma
        alphaLs[i] = 0.92;  // Stima iniziale per alpha1
        nLs[i] = 16.2;  // Stima iniziale per n1
        alphaRs[i] = 1.67;  // Stima iniziale per alpha2
        nRs[i] = 4.59;  // Stima iniziale per n2
        meanInf[i] = means[i] - 2 * sigmas[i];
        meanSup[i] =  means[i] + 2 * sigmas[i];
    }
    //personalizzo per bin 1 //scommentare per limiti massa regressed
    sigmas[0] = 0.1840;
    means[0] = 3;
    //bin2
    means[1] = 3.05;
    sigmas[1] = 0.1702;
    meanInf[1] = 3.01;
    //bin3 
    means[2] = 3.1;
    meanSup[2] = 3.2;
    sigmas[2] = 0.13;
    meanInf[2] = 3.04;
    //bin 4
    //bin 5
    sigmas[4] = 0.07;
    means[4] = 3.13;
    nRs[4] = 15;
    nLs[4] = 4;
    double mu[Nbins], inc_mu[Nbins], val_sigma[Nbins], inc_sigma[Nbins];
    double val_nL[Nbins], inc_nL[Nbins], val_alphL[Nbins], inc_alphL[Nbins], val_nR[Nbins], inc_nR[Nbins], val_alphR[Nbins], inc_alphR[Nbins]; //array per salvare i parametri

    //scommenta qui per limiti su massa regressed
    double massMin[Nbins] = {2.2, 2.25, 2.3, 2.1, 2.65, 2.65};  // Limiti inferiori personalizzati
    double massMax[Nbins] = {3.9, 3.9, 3.8, 3.8, 3.5, 3.5};  // Limiti superiori personalizzati
    double Bin_centers[Nbins] = {5.5, 8, 10.5, 12.5, 17, 30}; //Centri dei bin in P_t
    double Bin_halfwidth[Nbins] = {1.5, 1, 1.5, 1.5, 3, 10}; //Larghezze dei bin in P_t

    


    int rebin_factor[Nbins] = {1, 1, 2, 1, 1, 1};

    // Loop su tutti gli istogrammi
    for (int i = 0; i < Nbins; i++) {
        // Ottieni l'istogramma
        TH1D *hist = (TH1D*)file->Get(Form("proj_bin_%d", i + 1)); //qui per farlo con la massa regressed
        if (!hist) {
            std::cerr << "Istogramma non trovato: proj_bin_" << i + 1 << std::endl;
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
        RooRealVar cb_alphaL(Form("cb_alphaL_%d", i+1), "AlphaL of CB", alphaLs[i], 0.1, 15.0);
        RooRealVar cb_nL(Form("cb_nL_%d", i+1), "nL of CB", nLs[i], 0.1, 200.0);
        RooRealVar cb_alphaR(Form("cb_alphaR_%d", i+1), "AlphaR of CB", alphaRs[i], 0.1, 15.0);
        RooRealVar cb_nR(Form("cb_nR_%d", i+1), "nR of CB", nRs[i], 0.1, 100.0);

        // Crea la Crystal Ball asimmetrica 
        RooCrystalBall cb(Form("cb_%d", i+1), "Asymmetric Crystal Ball", mass, cb_mean, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR);

        // Esegui il fit
        RooFitResult *fit_result_simplex = cb.fitTo(data, RooFit::Minimizer("Minuit2", "Simplex"), RooFit::MaxCalls(1000000), RooFit::Save(),  RooFit::SumW2Error(true), RooFit::PrintLevel(-1));

        RooFitResult *fit_result = cb.fitTo(data, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::MaxCalls(10000000), RooFit::Save(), RooFit::SumW2Error(true), RooFit::PrintLevel(-1));
        ////////////////////////////////////////////////////////
        TCanvas *canvas = new TCanvas(Form("canvas_%d", i+1), Form("Fit for proj_bin_%d", i+1), 900, 700);
        //canvas->Divide(1,2);

        // Plotta l'istogramma e il fit
        RooPlot *frame = mass.frame();
        frame->SetTitle("");
        //frame->SetTitle(hist->GetTitle());
        data.plotOn(frame);
        cb.plotOn(frame, RooFit::LineColor(bluCMS), RooFit::LineWidth(5));
        frame->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");

        double chi2 = frame->chiSquare();
        // Stampare il chi-quadro sul plot
        TPaveText *paveText = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
        paveText->AddText(Form("#chi^{2} = %.2f", chi2));
        paveText->AddText(Form("#mu = %.4f +/- %.4f", cb_mean.getVal(), cb_mean.getError()));
        paveText->AddText(Form("#sigma = %.4f +/- %.4f", cb_sigma.getVal(), cb_sigma.getError()));
        if (fit_result->status() == 0) {
        paveText->AddText("Fit converged");
        } else {
        paveText->AddText("Fit did not converge");
        paveText->AddText(Form("Status: %d", fit_result->status()));
        }
        paveText->SetFillColor(0);
        frame->addObject(paveText);

        /*// Crea un frame per i residui
        
        RooPlot* pullFrame = mass.frame();
        pullFrame->SetTitle("Residui Normalizzati");
        // Calcola i residui normalizzati
        RooHist* pullHist = frame->pullHist();
        pullFrame->addPlotable(pullHist, "P");*/


        // Disegna il frame
        //canvas->cd(1);
        frame->Draw();
        WriteSimulation();

        /*//Disegno i residui
        canvas->cd(2);
        pullFrame->Draw();
        //linea per lo 0
    
        double xMin = pullFrame->GetXaxis()->GetXmin();
        double xMax = pullFrame->GetXaxis()->GetXmax();
        TLine *line = new TLine(xMin, 0, xMax, 0);
        line->SetLineStyle(2);  // Linea tratteggiata
        line->SetLineColor(kRed);  // Colore della linea, opzionale

        // Disegnare la linea sullo stesso canvas
        line->Draw("same");*/

        // Salva la canvas 
        canvas->SaveAs(Form("FitMC/fit_proj_bin_%d.png", i+1));
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

    TH1F *hist_inclusive = (TH1F*)file->Get("h_invmass_ECAL_3reweights");
    // Definisci la variabile di massa invariante
        RooRealVar mass("mass", "m(e^{+}e^{-})", 2.2, 3.9); 

        // Converto l'istogramma in un RooDataHist
        RooDataHist data_incl("data_incl", "Dataset from histogram", mass, hist_inclusive);

        // Definisci i parametri della Crystal Ball asimmetrica
        RooRealVar cb_mean("cb_mean", "Mean of CB", means[2], meanInf[2], meanSup[2]);
        RooRealVar cb_sigma("cb_sigma", "Sigma of CB", sigmas[2], 0.1 * sigmas[2], 2 * sigmas[2]);
        RooRealVar cb_alphaL("cb_alphaL", "AlphaL of CB", alphaLs[2], 0.1, 15.0);
        RooRealVar cb_nL("cb_nL", "nL of CB", nLs[2], 0.1, 200.0);
        RooRealVar cb_alphaR("cb_alphaR", "AlphaR of CB", alphaRs[2], 0.1, 15.0);
        RooRealVar cb_nR("cb_nR", "nR of CB", nRs[2], 0.1, 100.0);

        // Crea la Crystal Ball asimmetrica 
        RooCrystalBall cb("cb", "Asymmetric Crystal Ball", mass, cb_mean, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR);

        // Esegui il fit
        RooFitResult *fit_result_simplex = cb.fitTo(data_incl, RooFit::Minimizer("Minuit2", "Simplex"), RooFit::MaxCalls(1000000), RooFit::Save(),  RooFit::SumW2Error(true), RooFit::PrintLevel(-1));

        RooFitResult *fit_result = cb.fitTo(data_incl, RooFit::Minimizer("Minuit2", "Migrad"), RooFit::MaxCalls(10000000), RooFit::Save(), RooFit::SumW2Error(true), RooFit::PrintLevel(-1));
        ////////////////////////////////////////////////////////
        TCanvas *canvas_incl = new TCanvas("canvas_incl","Fit inclusive in p_{T}", 900, 700);
        RooPlot *frame = mass.frame();
        frame->SetTitle("");
        data_incl.plotOn(frame);
        cb.plotOn(frame, RooFit::LineColor(bluCMS), RooFit::LineWidth(5));
        frame->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");

        double chi2 = frame->chiSquare();
        // Stampare il chi-quadro sul plot
        TPaveText *paveText = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
        //paveText->AddText(Form("#chi^{2} = %.2f", chi2));
        paveText->AddText(Form("#mu = %.4f +/- %.4f", cb_mean.getVal(), cb_mean.getError()));
        paveText->AddText(Form("#sigma = %.4f +/- %.4f", cb_sigma.getVal(), cb_sigma.getError()));
        if (fit_result->status() == 0) {
        paveText->AddText("Fit converged");
        } else {
        paveText->AddText("Fit did not converge");
        paveText->AddText(Form("Status: %d", fit_result->status()));
        }
        paveText->SetFillColor(0);
        frame->addObject(paveText);
        frame->Draw();
        WriteSimulation();
        canvas_incl->SaveAs("FitMC/fit_inclusivePt.png");
        delete canvas_incl;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Salvo i parametri per inizializzare il fit dei dati
        // Crea un file ROOT
        TFile *file_param = TFile::Open("fit_results.root", "RECREATE");

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

