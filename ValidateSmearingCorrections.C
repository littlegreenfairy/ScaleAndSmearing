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

void ValidateSmearingCorrections() {
    ///////////////Definisco i colori
    int bluCMS = TColor::GetColor("#5790FC");
    int rosaCMS = TColor::GetColor("#964A8B");
    int rossoCMS = TColor::GetColor("#A42536");
    int gialloCMS = TColor::GetColor("#F89C20");
    ///////////////
    // Apri il file ROOT contenente gli istogrammi
    TFile *file = TFile::Open("outputHistograms_MC.root");

    // Prendi l'istogramma 2D h_scale_vs_pt
    TH2D* h2d = (TH2D*)file->Get("h_invMass_corrected_smearing");
    if (!h2d) {
        std::cerr << "Istogramma 2D h_invMass_corrected_smearing non trovato!" << std::endl;
        return;
    }

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
    

    double mu[Nbins], inc_mu[Nbins], val_sigma[Nbins], inc_sigma[Nbins];
    double val_nL[Nbins], inc_nL[Nbins], val_alphL[Nbins], inc_alphL[Nbins], val_nR[Nbins], inc_nR[Nbins], val_alphR[Nbins], inc_alphR[Nbins]; //array per salvare i parametri

    
    double massMin[Nbins] = {2.2, 2.2, 2.3, 2.1, 2.4, 2.4};  // Limiti inferiori personalizzati
    double massMax[Nbins] = {3.9, 3.9, 3.8,3.8, 3.8, 3.8};  // Limiti superiori personalizzati
    double Bin_centers[Nbins] = {5.75, 8.25, 10.5, 12.5, 17, 30}; //Centri dei bin in P_t
    double Bin_halfwidth[Nbins] = {1.75, 0.75, 1.5, 1.5, 3, 10}; //Larghezze dei bin in P_t

    


    int rebin_factor[Nbins] = {1, 1, 1, 1, 1, 2}; 

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
        RooRealVar cb_alphaL(Form("cb_alphaL_%d", i+1), "AlphaL of CB", alphaLs[i], 0, 15.0);
        RooRealVar cb_nL(Form("cb_nL_%d", i+1), "nL of CB", nLs[i], 0, 500.0);
        RooRealVar cb_alphaR(Form("cb_alphaR_%d", i+1), "AlphaR of CB", alphaRs[i], 0, 15.0);
        RooRealVar cb_nR(Form("cb_nR_%d", i+1), "nR of CB", nRs[i], 0, 500.0);

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
        frame->GetXaxis()->SetRangeUser(0, 5);
        data.plotOn(frame);
        cb.plotOn(frame, RooFit::LineColor(bluCMS), RooFit::LineWidth(5));
        frame->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");

        double chi2 = frame->chiSquare();
        // Stampare il chi-quadro sul plot
        TPaveText *paveText = new TPaveText(0.7, 0.75, 0.88, 0.88, "NDC");
        paveText->AddText(Form("#chi^{2} = %.2f", chi2));
        paveText->AddText(Form("#mu = %.4f +/- %.4f", cb_mean.getVal(), cb_mean.getError()));
        paveText->AddText(Form("#sigma = %.4f +/- %.4f", cb_sigma.getVal(), cb_sigma.getError()));
        //paveText->AddText(Form("#alpha_{L} = %.3f +/- %.3f", cb_alphaL.getVal(), cb_alphaL.getError()));
        //paveText->AddText(Form("#alpha_{R} = %.3f +/- %.3f", cb_alphaR.getVal(), cb_alphaR.getError()));
        //paveText->AddText(Form("n_{L} = %.3f +/- %.3f", cb_nL.getVal(), cb_nL.getError()));
        //paveText->AddText(Form("n_{R} = %.3f +/- %.3f", cb_nR.getVal(), cb_nR.getError()));
        paveText->SetFillColor(0);
        paveText->SetBorderSize(0);  // Remove border
        paveText->SetShadowColor(0); // Remove shadow
        frame->addObject(paveText);

        frame->Draw();
        WriteSimulation();

        // Salva la canvas 
        canvas->SaveAs(Form("PlotConID2022/FitMC_PostSmearingcorr/fit_proj_bin_%d.png", i+1));
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

    std::cout << "Fit completato per tutti i bin." << std::endl;

    TFile* fsmear = new TFile("smearing_corrections_fullele.root", "UPDATE");
if (fsmear && !fsmear->IsZombie()) {
    TGraphErrors* gr_smearing_vs_pt = (TGraphErrors*)fsmear->Get("gr_smearing_vs_pt");
    TGraphErrors* gr_sigma_vs_pt = (TGraphErrors*)fsmear->Get("gr_sigma_vs_pt");
    if (!gr_smearing_vs_pt || !gr_sigma_vs_pt) {
        std::cerr << "Error: Could not find required graphs in smearing_corrections_fullele.root!" << std::endl;
    } else {
        // Create the ratio graph
        TGraphErrors* gr_smearing_vs_pt_postsmearcorr = new TGraphErrors(Nbins);
        for (int i = 0; i < Nbins; ++i) {
            double pt, sigma_data, err_sigma_data;
            gr_sigma_vs_pt->GetPoint(i, pt, sigma_data);
            err_sigma_data = gr_sigma_vs_pt->GetErrorY(i);

            double sigma_cb = val_sigma[i];
            double err_sigma_cb = inc_sigma[i];

            double ratio = sigma_data / sigma_cb;
            // Error propagation: (Δr/r)^2 = (Δa/a)^2 + (Δb/b)^2
            double err_ratio = ratio * std::sqrt(
                std::pow(err_sigma_data / sigma_data, 2) +
                std::pow(err_sigma_cb / sigma_cb, 2)
            );

            gr_smearing_vs_pt_postsmearcorr->SetPoint(i, pt, ratio);
            gr_smearing_vs_pt_postsmearcorr->SetPointError(i, Bin_halfwidth[i], err_ratio);
        }
        gr_smearing_vs_pt_postsmearcorr->SetName("smearing_vs_pt_postsmearcorr");

        // Save the new graph
        gr_smearing_vs_pt_postsmearcorr->Write("smearing_vs_pt_postsmearcorr");
        // Plot both graphs
        TCanvas* c_smear = new TCanvas("c_smear", "Smearing Correction Comparison", 800, 600);
        gr_smearing_vs_pt->SetLineColor(kBlue);
        gr_smearing_vs_pt->SetMarkerColor(kBlue);
        gr_smearing_vs_pt->SetMarkerStyle(20);
        gr_smearing_vs_pt->SetTitle("Smearing Correction Comparison; p_{T} [GeV]; Smearing");

        gr_smearing_vs_pt->GetYaxis()->SetRangeUser(0.5, 1.5);
        gr_smearing_vs_pt->Draw("AP");
        gr_smearing_vs_pt_postsmearcorr->SetLineColor(kRed);
        gr_smearing_vs_pt_postsmearcorr->SetMarkerColor(kRed);
        gr_smearing_vs_pt_postsmearcorr->SetMarkerStyle(21);
        gr_smearing_vs_pt_postsmearcorr->Draw("P SAME");

        auto leg = new TLegend(0.6, 0.75, 0.88, 0.88);
        leg->AddEntry(gr_smearing_vs_pt, "Post scale corr", "lp");
        leg->AddEntry(gr_smearing_vs_pt_postsmearcorr, "Post smearing corr", "lp");
        leg->Draw();

        c_smear->SaveAs("smearing_comparison.png");

        delete c_smear;
        delete gr_smearing_vs_pt_postsmearcorr;
    }
    fsmear->Close();
    delete fsmear;
} else {
    std::cerr << "Error: Could not open smearing_corrections.root for update!" << std::endl;
}
    

    file->Close();
}

