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
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TF1.h>
#include <TColor.h>
#include <TLatex.h>
#define Nbinsrun 5
#define Nsigma 1

///////////////Definisco i colori
    int bluCMS = TColor::GetColor("#5790fc");
    int rosaCMS = TColor::GetColor("#964a8b");
    int rossoCMS = TColor::GetColor("#e42536");
    int gialloCMS = TColor::GetColor("#f89c20");
    int violaCMS = TColor::GetColor("#7a21dd");
    ///////////////

////....oooOO0OOooo........oooOO0OOooo.... PARAMETRI INIZIALI FIT (GLOBALI) ....oooOO0OOooo........oooOO0OOooo...
    // Definisci i limiti personalizzati e fattore di rebinning
    double LeftLowLim[Nbinsrun] = {1.6,1.4,1.7,1.5,1.7};  // Limiti sinistri personalizzati
    double LeftUpLim[Nbinsrun] = {2.5,2.4,2.3,2.4,2.3};
    double RightLowLim[Nbinsrun] = {3.6, 3.6, 3.6, 3.5, 3.5};
    double RightUpLim[Nbinsrun] = {4.5,5,5,5,4.5};  // Limiti destri personalizzati*/

    double BackgroundYlims[Nbinsrun] = {30000, 30000,30000,30000,30000};

    ///////////////// Parametri Crystal Ball
    double Mucb_ini[Nbinsrun] = {3.0842, 3.0337, 3.0422, 3.0489, 3.0382};
    double Mucb_lowlim[Nbinsrun] = {2.9, 2.9, 2.9, 2.9, 2.9};
    double Mucb_uplim[Nbinsrun] = {3.2, 3.15, 3.15, 3.15, 3.2};
    double Sigmacb_ini[Nbinsrun] = {0.1437, 0.1503, 0.1351, 0.1271, 0.1315};
    double Sigmacb_uplim[Nbinsrun] = {0.3, 0.5, 0.5, 0.5, 0.15};
    double Sigmacb_lowlim[Nbinsrun] = {0, 0.05, 0.05, 0.05, 0.05};

    ///////////////////////
    // Aggiungi gli array per i parametri della Gaussiana
    double gauss_mu_init[Nbinsrun] = {3.6, 3.65, 3.65, 3.65, 3.65};  // Valori iniziali per mu della Gaussiana
    double gauss_mu_low[Nbinsrun] = {3.5, 3.55, 3.5, 3.5, 3.2};   // Limiti inferiori per mu
    double gauss_mu_up[Nbinsrun] = {3.7, 3.8, 3.8, 3.8, 3.8};    // Limiti superiori per mu
    double gauss_sigma_init[Nbinsrun] = {0.1, 0.1, 0.1, 0.1, 0.1}; // Valori iniziali per sigma della Gaussiana
    double gauss_sigma_low[Nbinsrun] = {0.05, 0.05, 0.02, 0.05, 0.05}; // Limiti inferiori per sigma
    double gauss_sigma_up[Nbinsrun] = {0.2, 0.5, 0.2, 0.2, 0.2}; // Limiti superiori per sigma

////....oooOO0OOooo........oooOO0OOooo....oooOO0OOooo........oooOO0OOooo...
//// Centri dei bin e halfwidths
void Fit_binRunN(TH2D* invmass_vsrun, TGraphErrors *graph_mu_vsrun, double* bincenters, double* binhalfwidths); //restituisce l'andamento della media vs bin di run number

//....oooOO0OOooo........oooOO0OOooo.... FUNZIONE PRINCIPALE ....oooOO0OOooo........oooOO0OOooo...

void FitProfilePostcorr(){
    //Apro il file ed estraggo i 3 istogrammi 
    TFile *filehisto = TFile::Open("outputHistograms_DATA_partF.root");
    TH2D *h2D_runN_invM_check = (TH2D*)filehisto->Get("h_invMass_ECAL_check_yz");
    TH2D *h2D_runN_invM_corr = (TH2D*)filehisto->Get("h_invMass_ECAL_corrected_yz");
    TH2D *h2D_runN_invM_corr_offdiag = (TH2D*)filehisto->Get("h_invMass_ECAL_corr_offdiag_yz");
    //Limiti dei bin
    double binedges[6] = {360000, 360800, 361200, 361600, 362070, 362500};
    double bincenters[5];
    double binhalfwidths[5];
    for (int i = 0; i < 5; ++i) {
        bincenters[i] = (binedges[i] + binedges[i + 1]) / 2;
        binhalfwidths[i] = (binedges[i + 1] - binedges[i]) / 2;
    }


    //Dichiaro i TGraphErrors che conterranno gli andamenti vs run number
    std::cout << "print 1" << std::endl;
    TGraphErrors *graph_mucheck = new TGraphErrors(Nbinsrun);
    TGraphErrors *graph_mucorr = new TGraphErrors(Nbinsrun);
    TGraphErrors *graph_mucorr_offdiag = new TGraphErrors(Nbinsrun);
    //chiamo la funzione che fitta i profili su tutti e 3
    Fit_binRunN(h2D_runN_invM_check, graph_mucheck, bincenters, binhalfwidths);
    Fit_binRunN(h2D_runN_invM_corr, graph_mucorr, bincenters, binhalfwidths);
    Fit_binRunN(h2D_runN_invM_corr_offdiag, graph_mucorr_offdiag, bincenters, binhalfwidths);

    //plotto la media in funzione del bin per tutti e 3
    TCanvas *c = new TCanvas("c", "Mean vs run Number", 800, 800);
    graph_mucheck->SetTitle("Mean invariant mass vs Run Number");
    graph_mucheck->GetXaxis()->SetTitle("Run number");
    graph_mucheck->GetYaxis()->SetTitle("m_{J/#psi}");
    graph_mucheck->SetMinimum(3);
    graph_mucheck->SetMaximum(3.1);
    graph_mucheck->SetMarkerStyle(21);
    graph_mucheck->SetLineColor(rossoCMS);
    graph_mucheck->SetMarkerColor(rossoCMS);
    graph_mucheck->Draw("APE");
    graph_mucorr->SetMarkerStyle(21);
    graph_mucorr->SetLineColor(gialloCMS);
    graph_mucorr->SetMarkerColor(gialloCMS);
    graph_mucorr->Draw("PE SAME");
    graph_mucorr_offdiag->SetMarkerStyle(21);
    graph_mucorr_offdiag->SetLineColor(violaCMS);
    graph_mucorr_offdiag->SetMarkerColor(violaCMS);
    graph_mucorr_offdiag->Draw("PE SAME");

    // Add a legend
    TLegend *legend = new TLegend(0.6, 0.8, 0.9, 0.9);  // Position the legend in the top-right corner
    legend->AddEntry(graph_mucheck, "m_{J/#psi} without corrections", "p");
    legend->AddEntry(graph_mucorr, "corrected m_{J/#psi}", "p");
    legend->AddEntry(graph_mucorr_offdiag, "corrected m_{J/#psi} (only off-diag) ", "p");
    legend->Draw();

    c->Update();
    c->SaveAs("MeanvsRunNumber_fromfit.png");
    delete c;

    ////////////// Confronti Istogrammi inclusivi in pT e run number
    TH1F *minv_ECAL = (TH1F*)filehisto->Get("h_invMass_ECAL_ele");
    TH1F *minv_fromscratch = (TH1F*)filehisto->Get("h_invmass_fromscratch");
    TH1F *minv_fromscratch_corr = (TH1F*)filehisto->Get("h_invmass_fromscratch_corr");
    //prendo anche il MC
    TFile *filehisto_mc = TFile::Open("outputHistograms_MC.root");
    TH1F *minv_ECAL_MC = (TH1F*)filehisto_mc->Get("h_invMass_ECAL_ele");

    TCanvas *c2 = new TCanvas("c2", "Invariant mass inclusive in pt", 800, 800);

    minv_ECAL_MC->SetStats(kFALSE);
    minv_ECAL_MC->SetFillColorAlpha(bluCMS, 0.4);
    minv_ECAL_MC->SetLineColor(kBlack);
    minv_ECAL_MC->Draw("HISTO");

    minv_ECAL->SetStats(kFALSE);
    minv_ECAL->SetFillColorAlpha(rossoCMS, 0.4);
    minv_ECAL->SetLineColor(kBlack);
    minv_ECAL->Draw("HISTO SAME");

    minv_fromscratch->Scale(1.0/ minv_fromscratch->Integral());
    minv_fromscratch->SetStats(kFALSE);
    minv_fromscratch->SetFillColorAlpha(gialloCMS, 0.4);
    minv_fromscratch->SetLineColor(kBlack);
    //minv_fromscratch->Draw("HISTO SAME");

    minv_fromscratch_corr->Scale(1.0/ minv_fromscratch_corr->Integral());
    minv_fromscratch_corr->SetStats(kFALSE);
    minv_fromscratch_corr->SetFillColorAlpha(violaCMS, 0.4);
    minv_fromscratch_corr->SetLineColor(kBlack);
    minv_fromscratch_corr->Draw("HISTO SAME");

    c2->SaveAs("Inclusive_invmass_postcorr.png");
    delete c2;

    filehisto->Close();
    filehisto_mc->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
                                    //FUNZIONE CHE FITTA
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Fit_binRunN(TH2D* invmass_vsrun, TGraphErrors *graph_mu_vsrun, double* bincenters, double* binhalfwidths){
    //leggo parametri iniziali dal tree
    TFile *file_param = TFile::Open("fit_results.root", "READ");
    TTree *tree = (TTree*)file_param->Get("fitResults");
    if (!tree) {
        std::cerr << "Impossibile trovare il TTree fitResults." << std::endl;
        return;
    }

    // Variabili per contenere i valori letti dal TTree
    double nL_ini, inc_nL, alphaL_ini, inc_alphaL, nR_ini, inc_nR, alphaR_ini, inc_alphaR;
    // Collega le variabili ai rami del TTree
    tree->SetBranchAddress("nL", &nL_ini);
    tree->SetBranchAddress("nLError", &inc_nL);
    tree->SetBranchAddress("alphaL", &alphaL_ini);
    tree->SetBranchAddress("alphaLError", &inc_alphaL);
    tree->SetBranchAddress("nR", &nR_ini);
    tree->SetBranchAddress("nRError", &inc_nR);
    tree->SetBranchAddress("alphaR", &alphaR_ini);
    tree->SetBranchAddress("alphaRError", &inc_alphaR);
    //Estraggo le projection binwise
    TCanvas* c1 = new TCanvas("c1", "Projections of Invariant Mass in Run Number Bins", 1200, 800);
    c1->Divide(3,2);
    for(int i=1; i<=Nbinsrun; i++){
        TString projName = Form("%s_runbin_%d", invmass_vsrun->GetName(), i);
        TH1D* proj = invmass_vsrun->ProjectionX(projName, i, i);

        ///Fit della projection
        //....oooOO0OOooo........oooOO0OOooo.... FONDO ....oooOO0OOooo........oooOO0OOooo....
        RooRealVar mass("mass", "m(e^{+}e^{-})", 0, 6); 
        // Converto l'istogramma in un RooDataHist
        RooDataHist data("data", "Dataset from histogram", mass, proj);

        // Definisci i parametri del polinomio di secondo grado
        RooRealVar A(Form("A_%d", i), "4th deg coeff", 0,  -RooNumber::infinity(), RooNumber::infinity());
        RooRealVar B(Form("B_%d", i), "3rd deg coeff", 0,  -RooNumber::infinity(), RooNumber::infinity());
        RooRealVar C(Form("C_%d", i), "2nd deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
        RooRealVar D(Form("D_%d", i), "1st deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
        RooRealVar E(Form("E_%d", i), "0 deg coeff",  0,  -RooNumber::infinity(), RooNumber::infinity());
        RooPolynomial poly("poly", "Polynomial of 2nd degree", mass, RooArgList(A, B, C, D, E));

        // Definisci i parametri della Gaussiana
        RooRealVar gauss_mu(Form("gauss_mu_%d", i), "Gaussian mean", gauss_mu_init[i-1], gauss_mu_low[i-1], gauss_mu_up[i-1]);
        RooRealVar gauss_sigma(Form("gauss_sigma_%d", i), "Gaussian sigma", gauss_sigma_init[i-1], gauss_sigma_low[i-1], gauss_sigma_up[i-1]);
        RooGaussian gauss("gauss", "Gaussian component", mass, gauss_mu, gauss_sigma);

        // Definisci il modello di fondo combinando il polinomio con la Gaussiana
        RooRealVar frac_gauss("frac_gauss", "fraction of Gaussian", 0.3, 0.0, 1.0);
        RooAddPdf background("background", "Background Model", RooArgList(poly, gauss), RooArgList(frac_gauss));

        //Definisci il range per il fit del fondo
        mass.setRange("range1", LeftLowLim[i-1], LeftUpLim[i-1]);
        mass.setRange("range2", RightLowLim[i-1], RightUpLim[i-1]);

        // Esegui il fit del fondo
        RooFitResult *fit_result = background.fitTo(data, RooFit::Range("range1,range2"), RooFit::Save(), RooFit::SumW2Error(true));

        // Ottieni i valori dei parametri e i loro errori dal primo fit
        double A_val = A.getVal(), A_err = A.getError();
        double B_val = B.getVal(), B_err = B.getError();
        double C_val = C.getVal(), C_err = C.getError();
        double D_val = D.getVal(), D_err = D.getError();
        double E_val = E.getVal(), E_err = E.getError();
        double gauss_mu_val = gauss_mu.getVal(), gauss_mu_err = gauss_mu.getError();
        double gauss_sigma_val = gauss_sigma.getVal(), gauss_sigma_err = gauss_sigma.getError();

        //plot fondo
        // Crea una canvas per il fit del solo fondo
        TCanvas *canvas_bg = new TCanvas(Form("canvas_bg_%d", i), Form("Fit for proj_bin_%d - Background Only", i), 1000, 1200);
        canvas_bg->Divide(1, 2);

        // Plotta l'istogramma e il fit del solo fondo
        RooPlot *frame_bg = mass.frame(RooFit::Range("range1,range2"));
        frame_bg->SetTitle(Form("Fit Background Only - Bin %d", i));
        data.plotOn(frame_bg);
        background.plotOn(frame_bg, RooFit::LineColor(kMagenta), RooFit::LineWidth(4));
        frame_bg->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");

        // Imposta i limiti dell'asse y manualmente
        frame_bg->SetMinimum(0); // Limite inferiore
        frame_bg->SetMaximum(BackgroundYlims[i-1]); // Limite superiore

        double chi2_bg = frame_bg->chiSquare();
        // Stampare il chi-quadro sul plot
        TPaveText *paveText_bg = new TPaveText(0.8, 0.85, 0.9, 0.9, "NDC");
        paveText_bg->AddText(Form("#chi^{2} = %.2f", chi2_bg));
        if (fit_result->status() == 0) {
        paveText_bg->AddText("Fit converged");
        } else {
        paveText_bg->AddText("Fit did not converge");
        }
        canvas_bg->cd(1);
        frame_bg->Draw();
        paveText_bg->Draw();
        //legendbg->Draw();

        // Calcola e plotta i residui normalizzati
        RooHist *pulls_bg = frame_bg->pullHist();
        RooPlot *frame_pull_bg = mass.frame();
        frame_pull_bg->SetTitle("Normalized residuals");
        frame_pull_bg->addPlotable(pulls_bg, "P");
        frame_pull_bg->GetYaxis()->SetTitle("Pull");
        frame_pull_bg->GetYaxis()->SetTitleOffset(1.4);

        canvas_bg->cd(2);
        frame_pull_bg->Draw();

        // Salva la canvas del solo fondo
        canvas_bg->SaveAs(Form("FitVerifications/Background/FitBackgroundOnly_proj_bin_%d.png", i));
        delete canvas_bg;

        //....oooOO0OOooo........oooOO0OOooo.... FONDO + SEGNALE ....oooOO0OOooo........oooOO0OOooo....
        tree->GetEntry(i-1);
        mass.setRange("range_full", LeftLowLim[i-1], RightUpLim[i-1]);
        // Definisci i parametri della Crystal Ball con i vincoli richiesti
        RooRealVar mu_cb(Form("mu_cb_%d",i), "mu_cb", Mucb_ini[i-1], Mucb_lowlim[i-1], Mucb_uplim[i-1]);
        RooRealVar sigma_cb(Form("sigma_cb_%d", i), "sigma_cb", Sigmacb_ini[i-1], Sigmacb_lowlim[i-1], Sigmacb_uplim[i-1]);
        RooRealVar alphaL_cb(Form("alphaL_cb_%d", i), "alphaL_cb",  alphaL_ini, alphaL_ini - Nsigma*inc_alphaL, alphaL_ini + Nsigma*inc_alphaL);
        RooRealVar nL_cb(Form("nL_cb_%d", i), "nL_cb", nL_ini, nL_ini - Nsigma*inc_nL, nL_ini + Nsigma*inc_nL);
        RooRealVar alphaR_cb(Form("alphaR_cb_%d", i), "alphaR_cb", alphaR_ini, alphaR_ini - Nsigma*inc_alphaR, alphaR_ini + Nsigma*inc_alphaR);
        RooRealVar nR_cb(Form("nR_cb_%d", i), "nR_cb", nR_ini, nR_ini - Nsigma*inc_nR, nR_ini + Nsigma*inc_nR);


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

        // Definisci il modello segnale + fondo
        //RooRealVar nSignal("nSignal", "signal yield", 100, 0, 1000000);
        //RooRealVar nBackground("nBackground", "background yield", 100, 0, 1000000);
        RooRealVar frac("frac", "fraction of background", 0.5, 0.0, 1.0);
        RooAddPdf model("model", "signal + background", RooArgList(crystal, background), RooArgList(frac));

        // Esegui il fit segnale + fondo
        fit_result = model.fitTo(data, RooFit::Range("range_full"), RooFit::Save(), RooFit::SumW2Error(true));

        //....oooOO0OOooo........oooOO0OOooo.... PLOT DEL FIT ....oooOO0OOooo........oooOO0OOooo....
        TCanvas *c = new TCanvas(Form("c_bin_%d", i), Form("Bin %d", i), 950, 700);
         gPad->SetMargin(0.15, 0.1, 0.1, 0.1); 
        RooPlot* frame = mass.frame(RooFit::Range(LeftLowLim[i-1], RightUpLim[i-1]));
        std::string histTitle = std::string(proj->GetTitle());
        frame->SetTitle("");
        frame->GetXaxis()->SetTitle("Invariant Mass [GeV]");
        data.plotOn(frame);
        model.plotOn(frame, RooFit::Components(background), RooFit::LineStyle(kDotted), RooFit::LineColor(gialloCMS), RooFit::LineWidth(5)); //fondo
        model.plotOn(frame, RooFit::Components(crystal), RooFit::LineStyle(kDashed), RooFit::LineColor(violaCMS), RooFit::LineWidth(5));
        model.plotOn(frame, RooFit::LineColor(rossoCMS), RooFit::LineWidth(5));
        frame->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");
        //box per chi^2, mu e sigma
        double chi2 = frame->chiSquare();
        // Stampare il chi-quadro sul plot
        TPaveText *paveText = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
        paveText->AddText(Form("#chi^{2} = %.2f", chi2));
        paveText->AddText(Form("#mu = %.4f +/- %.4f", mu_cb.getVal(), mu_cb.getError()));
        paveText->AddText(Form("#sigma = %.4f +/- %.4f", sigma_cb.getVal(), sigma_cb.getError()));
        //verifico che il fit converga
        if (fit_result->status() == 0) {
        paveText->AddText("Fit converged");
        } else {
        paveText->AddText("Fit did not converge");
        }

        frame->addObject(paveText);
        frame->Draw();
         // Aggiunta della legenda
        TLegend *legend = new TLegend(0.15, 0.72, 0.3, 0.9);
        legend->SetHeader(histTitle.c_str(), "C");  // Aggiunge l'intestazione con il titolo dell'istogramma
        legend->AddEntry(frame->getObject(0), "Data", "P");  // Dati
        legend->AddEntry(frame->getObject(1), "Background Fit", "L");  // Fit del fondo
        legend->AddEntry(frame->getObject(2), "Signal Fit", "L");  // Fit del segnale
        legend->AddEntry(frame->getObject(3), "Combined Fit", "L");  // Fit combinato
        legend->Draw();  // Disegna la legenda
        c->SaveAs(Form("FitVerifications/%s_fitbin_%d.png", invmass_vsrun->GetName(), i));
        delete c;

        c1->cd(i);
        proj->Draw();

        //Riempio il Tgrapherrors
        graph_mu_vsrun->SetPoint(i-1, bincenters[i-1], mu_cb.getVal());
        graph_mu_vsrun->SetPointError(i-1, binhalfwidths[i-1], mu_cb.getError());

    }
    c1->SaveAs(Form("FitVerifications/%s_ProjectionsInRunNumberBins.png", invmass_vsrun->GetName()));
    delete c1;
}