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
#define NbinsPt 6
#define Nsigma 1

///////////////Definisco i colori
    int bluCMS = TColor::GetColor("#5790fc");
    int rosaCMS = TColor::GetColor("#964a8b");
    int rossoCMS = TColor::GetColor("#e42536");
    int gialloCMS = TColor::GetColor("#f89c20");
    int violaCMS = TColor::GetColor("#7a21dd");
    int arancione = TColor::GetColor("#ca5200");
    ///////////////
void Fit_binPt(TH2D* invmass_vsPt, TGraphErrors *graph_mu_vsPt, TGraphErrors *graph_sigma_vsPt, double* bincenters, double* binhalfwidths, const std::string& filename); //restituisce l'andamento della media vs bin di Pt

//....oooOO0OOooo........oooOO0OOooo.... FUNZIONE PRINCIPALE ....oooOO0OOooo........oooOO0OOooo...
void ValidationVsPt(){
    TFile *filehisto = TFile::Open("outputHistograms_DATA_partF.root");
    TH2D *h2D_Pt_invM_check = (TH2D*)filehisto->Get("h_invMass_ECAL_check_Pt");
    TH2D *h2D_Pt_invM_check_offdiag = (TH2D*)filehisto->Get("h_invMass_ECAL_check_offdiag_Pt");
    TH2D *h2D_Pt_invM_corr = (TH2D*)filehisto->Get("h_invMass_ECAL_corrected_Pt");
    TH2D *h2D_Pt_invM_corr_diag = (TH2D*)filehisto->Get("h_invMass_ECAL_corr_diag_Pt");
    TH2D *h2D_Pt_invM_corr_offdiag = (TH2D*)filehisto->Get("h_invMass_ECAL_corr_offdiag_Pt");

    //file per salvare le projection
    TFile *fileprojections;
    fileprojections = new TFile("Validation_projections_histo.root", "RECREATE");
    std::cout << "File is writable: " << fileprojections->IsWritable() << std::endl;
    fileprojections->Close(); //lo riapro dopo nella funzione
    delete fileprojections;
    double bincenters[NbinsPt] = {5.5, 8, 10, 12.5, 17, 30};
    double binhalfwidths[NbinsPt] = {1.5, 1, 1, 1.5, 3, 10};

    //graph con la media dei dati diagonale
    TFile *filecorr = TFile::Open("scale_corrections.root");
    TGraphErrors *graph_mudata_diag = (TGraphErrors*)filecorr->Get("Graph_MeanData");
    
    // Extract h_smearing_inclusiveRun from the corrections file
    TH1D *h_smearing_inclusiveRun = (TH1D*)filecorr->Get("h_smearing_inclusiveRun");
    if (!h_smearing_inclusiveRun) {
        std::cerr << "Error: Could not find h_smearing_inclusiveRun in scale_corrections.root" << std::endl;
    }
    
    filecorr->Close();
    delete filecorr;

    //Dichiaro i TGraphErrors che conterranno gli andamenti vs Pt
    std::cout << "print 1" << std::endl;
    TGraphErrors *graph_mucheck = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_mucheck_offdiag = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_mucorr = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_mucorr_offdiag = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_mucorr_diag = new TGraphErrors(NbinsPt);

    // Create TGraphErrors for sigma values
    TGraphErrors *graph_sigmacheck = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_sigmacheck_offdiag = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_sigmacorr = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_sigmacorr_offdiag = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_sigmacorr_diag = new TGraphErrors(NbinsPt);

    //chiamo la funzione che fitta i profili su tutti e 4
    Fit_binPt(h2D_Pt_invM_check, graph_mucheck, graph_sigmacheck, bincenters, binhalfwidths, "InitialParamsFit/fit_parameters_check.txt");
    Fit_binPt(h2D_Pt_invM_corr_offdiag, graph_mucorr_offdiag, graph_sigmacorr_offdiag, bincenters, binhalfwidths, "InitialParamsFit/fit_parameters_corr_offdiag.txt");
    Fit_binPt(h2D_Pt_invM_corr, graph_mucorr, graph_sigmacorr, bincenters, binhalfwidths, "InitialParamsFit/fit_parameters_corr.txt");
    Fit_binPt(h2D_Pt_invM_check_offdiag, graph_mucheck_offdiag, graph_sigmacheck_offdiag, bincenters, binhalfwidths, "InitialParamsFit/fit_parameters_check_offdiag.txt");
    Fit_binPt(h2D_Pt_invM_corr_diag, graph_mucorr_diag, graph_sigmacorr_diag, bincenters, binhalfwidths, "InitialParamsFit/fit_parameters_corr_diag.txt");

    //memorizzo le y dei graph in degli array
    const double* mucorr = graph_mucorr->GetY();
    const double* mucheck = graph_mucheck->GetY();
    const double* mucheck_offdiag = graph_mucheck_offdiag->GetY();
    const double* mucheck_diag = graph_mudata_diag->GetY();
    const double* mucorr_offdiag = graph_mucorr_offdiag->GetY();
    const double* mucorr_diag = graph_mucorr_diag->GetY();

    //Riempio un TGraph anche per il Monte Carlo e i ratio plot
    TGraphErrors *graph_mu_mc = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_sigma_mc = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_mu_mc_offdiag = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_sigma_mc_offdiag = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_mu_mc_all = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_sigma_mc_all = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_ratio_check = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_ratio_corr = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_ratio_offdiag_beforecorr = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_ratio_diag_beforecorr = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_ratio_offdiag_aftercorr = new TGraphErrors(NbinsPt);
    TGraphErrors *graph_ratio_diag_aftercorr = new TGraphErrors(NbinsPt);
    TFile *file_param_mc = TFile::Open("fit_results.root", "READ"); //MC diagonal
    TFile *file_param_mc_offdiag = TFile::Open("fit_results_offdiag.root", "READ"); //MC offdiagonal
    TFile *file_param_mc_all = TFile::Open("fit_results_all.root", "READ"); //MC all events


    TTree *tree = (TTree*)file_param_mc->Get("fitResults");
    TTree *tree_offdiag = (TTree*)file_param_mc_offdiag->Get("fitResults");
    TTree *tree_all = (TTree*)file_param_mc_all->Get("fitResults");
    if (!tree || !tree_offdiag || !tree_all) {
        std::cerr << "Impossibile trovare il TTree fitResults." << std::endl;
        return;
    }
    double mu_mc, mu_mc_Err, mu_mc_offdiag, mu_mc_Err_offdiag, mu_mc_all, mu_mc_Err_all, sigma_mc, sigma_mc_Err, sigma_mc_offdiag, sigma_mc_Err_offdiag, sigma_mc_all, sigma_mc_Err_all;
    tree->SetBranchAddress("Mean", &mu_mc);
    tree->SetBranchAddress("MeanError", &mu_mc_Err);
    tree->SetBranchAddress("Sigma", &sigma_mc);
    tree->SetBranchAddress("SigmaError", &sigma_mc_Err);
    tree_offdiag->SetBranchAddress("Mean", &mu_mc_offdiag);
    tree_offdiag->SetBranchAddress("MeanError", &mu_mc_Err_offdiag);
    tree_offdiag->SetBranchAddress("Sigma", &sigma_mc_offdiag);
    tree_offdiag->SetBranchAddress("SigmaError", &sigma_mc_Err_offdiag);
    tree_all->SetBranchAddress("Mean", &mu_mc_all);
    tree_all->SetBranchAddress("MeanError", &mu_mc_Err_all);
    tree_all->SetBranchAddress("Sigma", &sigma_mc_all);
    tree_all->SetBranchAddress("SigmaError", &sigma_mc_Err_all);

    //Canvas per le distribuzioni sovrapposte
    TCanvas *c_proj = new TCanvas("c_proj", "distributions before and after corrections", 1800, 1200);
    c_proj->Divide(3, 2);
    TCanvas *c_proj_diag = new TCanvas("c_proj_diag", "distributions before and after corrections", 1800, 1200); //solo diagonali
    c_proj_diag->Divide(3,2);
    TCanvas *c_proj_offdiag = new TCanvas("c_proj_offdiag", "distributions before and after corrections", 1800, 1200); //solo off diagonal
    c_proj_offdiag->Divide(3, 2);
    fileprojections = new TFile("Validation_projections_histo.root", "READ");

    for(int i=0; i<NbinsPt; i++){
        tree->GetEntry(i);
        tree_offdiag->GetEntry(i);
        tree_all->GetEntry(i);

        graph_mu_mc->SetPoint(i, bincenters[i], mu_mc);
        graph_mu_mc->SetPointError(i, binhalfwidths[i], mu_mc_Err);
        graph_sigma_mc->SetPoint(i, bincenters[i], sigma_mc);
        graph_sigma_mc->SetPointError(i, binhalfwidths[i], sigma_mc_Err);
        graph_mu_mc_offdiag->SetPoint(i, bincenters[i], mu_mc_offdiag);
        graph_mu_mc_offdiag->SetPointError(i, binhalfwidths[i], mu_mc_Err_offdiag);
        graph_sigma_mc_offdiag->SetPoint(i, bincenters[i], sigma_mc_offdiag);
        graph_sigma_mc_offdiag->SetPointError(i, binhalfwidths[i], sigma_mc_Err_offdiag);
        graph_mu_mc_all->SetPoint(i, bincenters[i], mu_mc_all);
        graph_mu_mc_all->SetPointError(i, binhalfwidths[i], mu_mc_Err_all);
        graph_sigma_mc_all->SetPoint(i, bincenters[i], sigma_mc_all);
        graph_sigma_mc_all->SetPointError(i, binhalfwidths[i], sigma_mc_Err_all);

        //ratio
        graph_ratio_check->SetPoint(i, bincenters[i], mu_mc_all/mucheck[i]); // MC / dati non corretti
        graph_ratio_corr->SetPoint(i, bincenters[i], mu_mc_all/mucorr[i]); // MC / dati corretti

        // MC/dati post correzione - on and off diagonal
        graph_ratio_offdiag_aftercorr->SetPoint(i, bincenters[i], mu_mc_offdiag/mucorr_offdiag[i]);
        graph_ratio_diag_aftercorr->SetPoint(i, bincenters[i], mu_mc/mucorr_diag[i]);

        //  MC/dati pre correzione - on and off diagonal
        graph_ratio_offdiag_beforecorr->SetPoint(i, bincenters[i], mu_mc_offdiag/mucheck_offdiag[i]);
        graph_ratio_diag_beforecorr->SetPoint(i, bincenters[i], mu_mc/mucheck_diag[i]);
        //errori 
        graph_ratio_check->SetPointError(i, binhalfwidths[i], sqrt((mu_mc_Err_all/mu_mc_all)*(mu_mc_Err_all/mu_mc_all) + (graph_mucheck->GetErrorY(i)/mucheck[i])*(graph_mucheck->GetErrorY(i)/mucheck[i])) * mu_mc_all/mucheck[i]);
        graph_ratio_corr->SetPointError(i, binhalfwidths[i], sqrt((mu_mc_Err_all/mu_mc_all)*(mu_mc_Err_all/mu_mc_all) + (graph_mucorr->GetErrorY(i)/mucorr[i])*(graph_mucorr->GetErrorY(i)/mucheck[i])) * mu_mc_all/mucorr[i]);
        
        graph_ratio_offdiag_aftercorr->SetPointError(i, binhalfwidths[i], sqrt((graph_mucorr_offdiag->GetErrorY(i)/mucorr_offdiag[i])*(graph_mucorr_offdiag->GetErrorY(i)/mucorr_offdiag[i]) + (mu_mc_Err_offdiag/mu_mc_offdiag)*(mu_mc_Err_offdiag/mu_mc_offdiag)) * mu_mc_offdiag/mucorr_offdiag[i]);
        graph_ratio_diag_aftercorr->SetPointError(i, binhalfwidths[i], sqrt((graph_mucorr_diag->GetErrorY(i)/mucorr_diag[i])*(graph_mucorr_diag->GetErrorY(i)/mucorr_diag[i]) + (mu_mc_Err/mu_mc)*(mu_mc_Err/mu_mc)) * mu_mc/mucorr_diag[i]);
    
        graph_ratio_offdiag_beforecorr->SetPointError(i, binhalfwidths[i], sqrt((mu_mc_Err_offdiag/mu_mc_offdiag)*(mu_mc_Err_offdiag/mu_mc_offdiag) + (graph_mucheck_offdiag->GetErrorY(i)/mucheck_offdiag[i])*(graph_mucheck_offdiag->GetErrorY(i)/mucheck_offdiag[i])) * mu_mc_offdiag/mucheck_offdiag[i]);
        graph_ratio_diag_beforecorr->SetPointError(i, binhalfwidths[i], sqrt((mu_mc_Err/mu_mc)*(mu_mc_Err/mu_mc) + (graph_mudata_diag->GetErrorY(i)/mucheck_diag[i])*(graph_mudata_diag->GetErrorY(i)/mucheck_diag[i])) * mu_mc/mucheck_diag[i]);
        //Estraggo e disegno le distribuzioni da sovrapporre
        TH1D *proj_check = (TH1D*)fileprojections->Get(Form("h_invMass_ECAL_check_Pt_Ptbin_%d", i+1));
        TH1D *proj_corr = (TH1D*)fileprojections->Get(Form("h_invMass_ECAL_corrected_Pt_Ptbin_%d", i+1));
        TH1D *proj_corr_diag = (TH1D*)fileprojections->Get(Form("h_invMass_ECAL_corr_diag_Pt_Ptbin_%d", i+1));
        TH1D *proj_corr_offdiag = (TH1D*)fileprojections->Get(Form("h_invMass_ECAL_corr_offdiag_Pt_Ptbin_%d", i+1));

        c_proj->cd(i+1);
        proj_check->SetStats(kFALSE);
        proj_check->SetFillColorAlpha(rossoCMS, 0.3);
        proj_check->SetLineColor(rossoCMS);
        proj_check->Draw("HISTO");
        proj_corr->SetStats(kFALSE);
        proj_corr->SetFillColorAlpha(violaCMS, 0.3);
        proj_corr->SetLineColor(violaCMS);
        proj_corr->Draw("HISTO SAME");
        c_proj->Update();

        c_proj_diag->cd(i+1);
        proj_check->SetStats(kFALSE);
        proj_check->SetFillColorAlpha(rossoCMS, 0.3);
        proj_check->SetLineColor(rossoCMS);
        proj_check->Draw("HISTO");
        proj_corr_diag->SetStats(kFALSE);
        proj_corr_diag->SetFillColorAlpha(rosaCMS, 0.3);
        proj_corr_diag->SetLineColor(rosaCMS);
        proj_corr_diag->Draw("HISTO SAME");
        c_proj_diag->Update();

        c_proj_offdiag->cd(i+1);
        proj_check->SetStats(kFALSE);
        proj_check->SetFillColorAlpha(rossoCMS, 0.3);
        proj_check->SetLineColor(rossoCMS);
        proj_check->Draw("HISTO");
        proj_corr_offdiag->SetStats(kFALSE);
        proj_corr_offdiag->SetFillColorAlpha(violaCMS, 0.3);
        proj_corr_offdiag->SetLineColor(violaCMS);
        proj_corr_offdiag->Draw("HISTO SAME");
        c_proj_offdiag->Update();


    }
    //salvo distribuzioni sovrapposte
    c_proj->SaveAs("PlotConID2022/FitVerificationsPt/Overlapped_dist_corrected.png");
    c_proj_diag->SaveAs("PlotConID2022/FitVerificationsPt/Overlapped_dist_corrected_diag.png");
    c_proj_offdiag->SaveAs("PlotConID2022/FitVerificationsPt/Overlapped_dist_corrected_offdiag.png");
    delete c_proj;
    delete c_proj_diag;
    delete c_proj_offdiag;
    fileprojections->Close();

    //tolgo il primo punto dagli off-diagonal perché non è fisico
    graph_mucheck_offdiag->SetPoint(0, 0, 0);
    graph_mucorr_offdiag->SetPoint(0, 0, 0);

    // For the off-diagonal plots
    graph_mucheck_offdiag->GetXaxis()->SetLimits(0, 45);  // You already have this

    // Then match it in the ratio plot
    graph_ratio_offdiag_beforecorr->GetXaxis()->SetLimits(0, 45);  // Use same limits

    //plotto la media in funzione del bin per tutti e 3
    TCanvas *c = new TCanvas("c", "Mean vs p_{T}", 1600, 1200);
    //splitto il canvas in plot e ratio plot
    TPad *pad1 = new TPad("pad1", "Pad sinistra", 0.0, 0.35, 1.0, 1.0);
    pad1->SetBottomMargin(0.08); 
    pad1->SetLeftMargin(0.13);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2", "Pad inferiore (30%)", 0.0, 0.0, 1.0, 0.35);
    pad2->SetTopMargin(0.04);       
    pad2->Draw();

    pad1->cd();
    graph_mucheck->SetTitle("Mean invariant mass vs p_{T}");
    graph_mucheck->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_mucheck->GetYaxis()->SetTitle("m_{J/#psi} [GeV]");
    graph_mucheck->SetMinimum(3);
    graph_mucheck->SetMaximum(3.15);
    graph_mucheck->GetXaxis()->SetLimits(0, 45);  // Set explicit x-axis limits
    graph_mucheck->SetMarkerStyle(21);
    graph_mucheck->SetLineColor(rossoCMS);
    graph_mucheck->SetMarkerColor(rossoCMS);
    graph_mucheck->Draw("APE");
    graph_mucorr->SetMarkerStyle(21);
    graph_mucorr->SetLineColor(violaCMS);
    graph_mucorr->SetMarkerColor(violaCMS);
    graph_mucorr->Draw("PE SAME");
    graph_mucorr_offdiag->SetMarkerStyle(21);
    graph_mucorr_offdiag->SetLineColor(violaCMS);
    graph_mucorr_offdiag->SetMarkerColor(violaCMS);
    //graph_mucorr_offdiag->Draw("PE SAME");
    graph_mucorr_diag->SetMarkerStyle(21);
    graph_mucorr_diag->SetLineColor(violaCMS);
    graph_mucorr_diag->SetMarkerColor(violaCMS);
    //graph_mucorr_diag->Draw("PE SAME");
    //MC
    graph_mu_mc_all->SetMarkerStyle(21);
    graph_mu_mc_all->SetLineColor(bluCMS);
    graph_mu_mc_all->SetMarkerColor(bluCMS);
    graph_mu_mc_all->Draw("PE SAME");

    // Add a legend
    TLegend *legend = new TLegend(0.6, 0.8, 0.9, 0.9);  // Position the legend in the top-right corner
    legend->AddEntry(graph_mucheck, "m_{J/#psi} without corrections", "p");
    legend->AddEntry(graph_mucorr, "corrected m_{J/#psi}", "p");
    //legend->AddEntry(graph_mucorr_offdiag, "corrected m_{J/#psi} (only off-diag) ", "p");
    //legend->AddEntry(graph_mucorr_diag, "corrected m_{J/#psi} (only diagonal) ", "p");
    legend->AddEntry(graph_mu_mc_all, "m_{J/#psi} from MC", "p");
    legend->Draw();

    pad2->cd(); //ratio plot
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.13);
    graph_ratio_check->SetTitle("");
    graph_ratio_check->GetYaxis()->SetTitle("MC / data");
    graph_ratio_check->GetYaxis()->SetTitleSize(0.07); 
    graph_ratio_check->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_ratio_check->GetXaxis()->SetTitleSize(0.07);
    graph_ratio_check->GetXaxis()->SetLimits(0, 45);  // Use same limits as main plot
    graph_ratio_check->SetMarkerStyle(20);
    graph_ratio_check->SetLineColor(rossoCMS);
    graph_ratio_check->SetMarkerColor(rossoCMS);
    graph_ratio_check->SetMaximum(1.03);
    //graph_ratio_check->SetMinimum(0.996);
    graph_ratio_check->Draw("APE");

    graph_ratio_corr->SetTitle("");
    graph_ratio_corr->SetMarkerStyle(20);
    graph_ratio_corr->SetLineColor(violaCMS);
    graph_ratio_corr->SetMarkerColor(violaCMS);
    graph_ratio_corr->Draw("PE SAME");

    //riga per l'1
    TLine *Line = new TLine(graph_ratio_check->GetXaxis()->GetXmin(), 1,
                            graph_ratio_check->GetXaxis()->GetXmax(), 1);
    Line->SetLineColor(bluCMS);
    Line->SetLineWidth(2);
    Line->Draw();

    TLegend *legend_ratio = new TLegend(0.7, 0.75, 0.9, 0.92);  
    //legend_ratio->SetBorderSize(0);
    legend_ratio->AddEntry(graph_ratio_check, "before corrections", "p");
    legend_ratio->AddEntry(graph_ratio_corr, "after corrections", "p");
    legend_ratio->Draw();

    c->Update();
    c->SaveAs("PlotConID2022/CorrectionsValidation_vsPt.png");
    delete c;

    //Confronto diagonal e off diagonal
    TCanvas *c2 = new TCanvas("c2", "Mean vs p_{T}", 1600, 1200);
    TPad *padup = new TPad("padup", "Pad sinistra", 0.0, 0.35, 1.0, 1.0);
    padup->SetBottomMargin(0.08); 
    padup->Draw();
    TPad *paddown = new TPad("paddown", "Pad inferiore (30%)", 0.0, 0.0, 1.0, 0.35);
    paddown->SetTopMargin(0.04);       
    paddown->Draw();
    padup->cd();
    gPad->SetLeftMargin(0.1);
    graph_mucheck_offdiag->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_mucheck_offdiag->GetYaxis()->SetTitle("m_{J/#psi} [GeV]");
    graph_mucheck_offdiag->GetXaxis()->SetLimits(0, 45);
    graph_mucheck_offdiag->SetMinimum(3);
    graph_mucheck_offdiag->SetMaximum(3.15);
    graph_mucheck_offdiag->SetMarkerStyle(21);
    graph_mucheck_offdiag->SetLineColor(rossoCMS);
    graph_mucheck_offdiag->SetMarkerColor(rossoCMS);
    graph_mucheck_offdiag->Draw("APE");
    graph_mucorr_offdiag->Draw("PE SAME");
    graph_mu_mc_offdiag->SetMarkerStyle(21);
    graph_mu_mc_offdiag->SetLineColor(bluCMS);
    graph_mu_mc_offdiag->SetMarkerColor(bluCMS);
    graph_mu_mc_offdiag->Draw("PE SAME");

    TLegend *legend2 = new TLegend(0.6, 0.8, 0.9, 0.9);
    legend2->AddEntry(graph_mucheck_offdiag, "m_{J/#psi} without corrections", "p");
    legend2->AddEntry(graph_mucorr_offdiag, "corrected m_{J/#psi} (only off-diag) ", "p");
    legend2->AddEntry(graph_mu_mc_offdiag, "m_{J/#psi} from MC (only off-diag)", "p");
    legend2->Draw();

    paddown->cd(); //ratio plot
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.13);
    graph_ratio_offdiag_beforecorr->SetTitle("");
    graph_ratio_offdiag_beforecorr->GetYaxis()->SetTitle("MC / data");
    graph_ratio_offdiag_beforecorr->GetYaxis()->SetTitleSize(0.07);
    graph_ratio_offdiag_beforecorr->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_ratio_offdiag_beforecorr->GetXaxis()->SetTitleSize(0.07);
    graph_ratio_offdiag_beforecorr->GetXaxis()->SetLimits(0, 45);  // Use same limits
    graph_ratio_offdiag_beforecorr->SetMarkerStyle(20);
    graph_ratio_offdiag_beforecorr->SetLineColor(rossoCMS);
    graph_ratio_offdiag_beforecorr->SetMarkerColor(rossoCMS);
    graph_ratio_offdiag_beforecorr->SetMinimum(0.99);
    //graph_ratio_offdiag_beforecorr->SetMaximum(1.05);
    graph_ratio_offdiag_beforecorr->Draw("APE");

    graph_ratio_offdiag_aftercorr->SetTitle("");
    graph_ratio_offdiag_aftercorr->SetMarkerStyle(20);
    graph_ratio_offdiag_aftercorr->SetLineColor(violaCMS);
    graph_ratio_offdiag_aftercorr->SetMarkerColor(violaCMS);
    graph_ratio_offdiag_aftercorr->Draw("PE SAME");

    //riga per l'1
    TLine *Line2 = new TLine(graph_ratio_offdiag_beforecorr->GetXaxis()->GetXmin(), 1,
                            graph_ratio_offdiag_beforecorr->GetXaxis()->GetXmax(), 1);
    Line2->SetLineColor(bluCMS);
    Line2->SetLineWidth(2);
    Line2->Draw();

    TLegend *legend_ratio2 = new TLegend(0.7, 0.75, 0.9, 0.92);
    //legend_ratio2->SetBorderSize(0);
    legend_ratio2->AddEntry(graph_ratio_offdiag_beforecorr, "before corrections", "p");
    legend_ratio2->AddEntry(graph_ratio_offdiag_aftercorr, "after corrections", "p");
    legend_ratio2->Draw();

    c2->Update();
    c2->SaveAs("PlotConID2022/Offdiagonal_vsPt.png");
    delete c2;

    TCanvas *c3 = new TCanvas("c3", "Mean vs p_{T}", 1600, 1200);
    TPad *pad = new TPad("pad", "Pad sinistra", 0.0, 0.35, 1.0, 1.0);
    pad->SetBottomMargin(0.08); 
    pad->Draw();
    TPad *padratio = new TPad("padratio", "Pad inferiore (30%)", 0.0, 0.0, 1.0, 0.35);
    padratio->SetTopMargin(0.04);       
    padratio->Draw();
    pad->cd();
    gPad->SetLeftMargin(0.13);
    //graph_mucheck->Draw("APE");
    graph_mudata_diag->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_mudata_diag->GetYaxis()->SetTitle("m_{J/#psi} [GeV]");
    graph_mudata_diag->GetXaxis()->SetLimits(0, 45);  // Add this line
    graph_mudata_diag->SetLineColor(rossoCMS);
    graph_mudata_diag->SetMarkerStyle(21);
    graph_mudata_diag->SetMarkerColor(rossoCMS);
    graph_mudata_diag->SetMinimum(3);
    graph_mudata_diag->SetMaximum(3.15);
    graph_mudata_diag->Draw("APE");
    graph_mucorr_diag->Draw("PE SAME");
    graph_mu_mc->SetMarkerStyle(21);
    graph_mu_mc->SetLineColor(bluCMS);
    graph_mu_mc->SetMarkerColor(bluCMS);
    graph_mu_mc->Draw("PE SAME");

    TLegend *legend3 = new TLegend(0.6, 0.8, 0.9, 0.9);
    legend3->AddEntry(graph_mucheck, "m_{J/#psi} without corrections", "p");
    legend3->AddEntry(graph_mucorr_diag, "corrected m_{J/#psi} (only diagonal) ", "p");
    legend3->AddEntry(graph_mu_mc, "m_{J/#psi} from MC (only diagonal)", "p");
    legend3->Draw();
    padratio->cd();
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    graph_ratio_diag_beforecorr->SetTitle("");
    graph_ratio_diag_beforecorr->GetYaxis()->SetTitle("MC / data");
    graph_ratio_diag_beforecorr->GetYaxis()->SetTitleSize(0.07);
    graph_ratio_diag_beforecorr->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_ratio_diag_beforecorr->GetXaxis()->SetTitleSize(0.07);
    graph_ratio_diag_beforecorr->GetXaxis()->SetLimits(0, 45);  // Use same limits
    graph_ratio_diag_beforecorr->SetMarkerStyle(20);
    graph_ratio_diag_beforecorr->SetLineColor(rossoCMS);
    graph_ratio_diag_beforecorr->SetMarkerColor(rossoCMS);
    //graph_ratio_diag_beforecorr->SetMaximum(1.03);
    graph_ratio_diag_beforecorr->SetMinimum(0.996);
    graph_ratio_diag_beforecorr->Draw("APE");

    graph_ratio_diag_aftercorr->SetTitle("");
    graph_ratio_diag_aftercorr->SetMarkerStyle(20);
    graph_ratio_diag_aftercorr->SetLineColor(violaCMS);
    graph_ratio_diag_aftercorr->SetMarkerColor(violaCMS);
    graph_ratio_diag_aftercorr->Draw("PE SAME");

    //riga per l'1
    TLine *Line3 = new TLine(graph_ratio_diag_beforecorr->GetXaxis()->GetXmin(), 1,
                            graph_ratio_diag_beforecorr->GetXaxis()->GetXmax(), 1);
    Line3->SetLineColor(bluCMS);
    Line3->SetLineWidth(2);
    Line3->Draw();

    TLegend *legend_ratio3 = new TLegend(0.7, 0.75, 0.9, 0.92);
    //legend_ratio3->SetBorderSize(0);
    legend_ratio3->AddEntry(graph_ratio_diag_beforecorr, "before corrections", "p");
    legend_ratio3->AddEntry(graph_ratio_diag_aftercorr, "after corrections", "p");
    legend_ratio3->Draw();

    c3->Update();
    c3->SaveAs("PlotConID2022/Diagonal_vsPt.png");
    delete c3;

    //// Smearing
    TCanvas *c4 = new TCanvas("c4", "Smearing vs p_{T}", 1600, 900); 
    
    // Create a new TGraphErrors to hold the ratio of sigmas
    TGraphErrors *graph_ratio_sigma = new TGraphErrors(NbinsPt);
    
    // Fill the ratio graph with sigma_data/sigma_MC values
    for(int i=0; i<NbinsPt; i++) {
        double sigma_data = graph_sigmacorr->GetY()[i];
        double sigma_data_err = graph_sigmacorr->GetErrorY(i);
        double sigma_mc = graph_sigma_mc_all->GetY()[i];
        double sigma_mc_err = graph_sigma_mc_all->GetErrorY(i);
        
        // Calculate the ratio and its error
        double ratio = sigma_data/sigma_mc;
        double ratio_err = ratio * sqrt(
            (sigma_data_err/sigma_data)*(sigma_data_err/sigma_data) + 
            (sigma_mc_err/sigma_mc)*(sigma_mc_err/sigma_mc)
        );
        
        graph_ratio_sigma->SetPoint(i, bincenters[i], ratio);
        graph_ratio_sigma->SetPointError(i, binhalfwidths[i], ratio_err);
    }
    TGraphErrors *graph_ratio_sigma_before = new TGraphErrors(NbinsPt);
    
    // Fill the ratio graph with sigma_data/sigma_MC values before corrections
    for(int i=0; i<NbinsPt; i++) {
        // Instead of calculating the ratio from scratch, use values from h_smearing_inclusiveRun
        if (h_smearing_inclusiveRun) {
            // Find the bin in h_smearing_inclusiveRun that corresponds to bincenters[i]
            int smearing_bin = h_smearing_inclusiveRun->FindBin(bincenters[i]);
            double ratio = h_smearing_inclusiveRun->GetBinContent(smearing_bin);
            double ratio_err = h_smearing_inclusiveRun->GetBinError(smearing_bin);
            
            graph_ratio_sigma_before->SetPoint(i, bincenters[i], ratio);
            graph_ratio_sigma_before->SetPointError(i, binhalfwidths[i], ratio_err);
        } else {
            // Fallback to the original calculation if h_smearing_inclusiveRun is not found
            double sigma_data = graph_sigmacheck->GetY()[i];
            double sigma_data_err = graph_sigmacheck->GetErrorY(i);
            double sigma_mc = graph_sigma_mc_all->GetY()[i];
            double sigma_mc_err = graph_sigma_mc_all->GetErrorY(i);
            
            // Calculate the ratio and its error
            double ratio = sigma_data/sigma_mc;
            double ratio_err = ratio * sqrt(
                (sigma_data_err/sigma_data)*(sigma_data_err/sigma_data) + 
                (sigma_mc_err/sigma_mc)*(sigma_mc_err/sigma_mc)
            );
            
            graph_ratio_sigma_before->SetPoint(i, bincenters[i], ratio);
            graph_ratio_sigma_before->SetPointError(i, binhalfwidths[i], ratio_err);
        }
    }
    
    // Set up the graph appearance
    graph_ratio_sigma->SetTitle("Ratio of resolution in data to MC vs p_{T}");
    graph_ratio_sigma->GetXaxis()->SetTitle("p_{T} [GeV]");
    graph_ratio_sigma->GetYaxis()->SetTitle("#sigma_{data}/#sigma_{MC}");
    graph_ratio_sigma->GetXaxis()->SetLimits(0, 45);
    graph_ratio_sigma->SetMarkerStyle(21);
    graph_ratio_sigma->SetLineColor(violaCMS);
    graph_ratio_sigma->SetMarkerColor(violaCMS);
    graph_ratio_sigma_before->SetMarkerStyle(21);
    graph_ratio_sigma_before->SetLineColor(rossoCMS);
    graph_ratio_sigma_before->SetMarkerColor(rossoCMS);
    
    // Draw the graph
    graph_ratio_sigma->Draw("APE");
    graph_ratio_sigma_before->Draw("PE SAME");
    
    // Add a horizontal line at y=1 for reference
    TLine *smearingLine = new TLine(graph_ratio_sigma->GetXaxis()->GetXmin(), 1,
                           graph_ratio_sigma->GetXaxis()->GetXmax(), 1);
    smearingLine->SetLineColor(bluCMS);
    smearingLine->SetLineWidth(2);
    smearingLine->SetLineStyle(kDashed);
    smearingLine->Draw();
    
    // Add a legend - moved to lower right corner
    TLegend *legend_smearing = new TLegend(0.7, 0.8, 0.9, 0.95);
    legend_smearing->AddEntry(graph_ratio_sigma, "#sigma_{data}/#sigma_{MC} after corrections", "p");
    legend_smearing->AddEntry(graph_ratio_sigma_before, "#sigma_{data}/#sigma_{MC} before corrections", "p");
    legend_smearing->AddEntry(smearingLine, "Perfect agreement", "l");
    legend_smearing->Draw();
    
    c4->Update();
    c4->SaveAs("PlotConID2022/SmearingPostScalecorr_vsPt.png");
    delete c4;
    delete graph_ratio_sigma;

    double x_val;
    double ratiomax_before = (graph_ratio_check->GetPoint(1, x_val, ratiomax_before), ratiomax_before);
    double ratiomax_after = (graph_ratio_corr->GetPoint(1, x_val, ratiomax_after), ratiomax_after);
    std::cout << "Massima distanza da 1 pre-corr: " << ratiomax_before << std::endl;
    std::cout << "Massima distanza da 1 post-corr: " << ratiomax_after << std::endl;

    filehisto->Close();
    delete filehisto;
}

/////////////////////////////////////////////////////////////
void Fit_binPt(TH2D* invmass_vsPt, TGraphErrors *graph_mu_vsPt, TGraphErrors *graph_sigma_vsPt, double* bincenters, double* binhalfwidths, const std::string& filename){
    //leggo parametri iniziali dal tree
    TFile *file_param_mc = TFile::Open("fit_results.root", "READ");
    TTree *tree = (TTree*)file_param_mc->Get("fitResults");
    if (!tree) {
        std::cerr << "Impossibile trovare il TTree fitResults." << std::endl;
        return;
    }
    std::cout << "print 2" << std::endl;
    //apro il file su cui salvo le varie projections
    TFile *fileprojections = TFile::Open("Validation_projections_histo.root", "UPDATE");
    fileprojections->cd();
    // Variables to read parameters from file
    std::vector<double> LeftLowLim, LeftUpLim, RightLowLim, RightUpLim, BackgroundYlims;
    std::vector<double> Mucb_ini, Mucb_lowlim, Mucb_uplim, Sigmacb_ini, Sigmacb_lowlim, Sigmacb_uplim;
    std::vector<double> gauss_mu_init, gauss_mu_low, gauss_mu_up, gauss_sigma_init, gauss_sigma_low, gauss_sigma_up;

    std::ifstream paramFile(filename);
    if (!paramFile.is_open()) {
        std::cerr << "Unable to open parameter file!" << std::endl;
        return;
    }
    std::cout << "print 3" << std::endl;
    // Salta la prima riga (titoli delle colonne)
    std::string headerLine;
    std::getline(paramFile, headerLine);
    //Read parameters from file
    double temp;
    while (paramFile >> temp) {
        LeftLowLim.push_back(temp);
        paramFile >> temp; LeftUpLim.push_back(temp);
        paramFile >> temp; RightLowLim.push_back(temp);
        paramFile >> temp; RightUpLim.push_back(temp);
        paramFile >> temp; BackgroundYlims.push_back(temp);
        paramFile >> temp; Mucb_ini.push_back(temp);
        paramFile >> temp; Mucb_lowlim.push_back(temp);
        paramFile >> temp; Mucb_uplim.push_back(temp);
        paramFile >> temp; Sigmacb_ini.push_back(temp);
        paramFile >> temp; Sigmacb_lowlim.push_back(temp);
        paramFile >> temp; Sigmacb_uplim.push_back(temp);
        paramFile >> temp; gauss_mu_init.push_back(temp);
        paramFile >> temp; gauss_mu_low.push_back(temp);
        paramFile >> temp; gauss_mu_up.push_back(temp);
        paramFile >> temp; gauss_sigma_init.push_back(temp);
        paramFile >> temp; gauss_sigma_low.push_back(temp);
        paramFile >> temp; gauss_sigma_up.push_back(temp);
    }
    paramFile.close();
    
    std::cout << "print 4" << std::endl;
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

    double compatibility_psi2s[NbinsPt];

    for(int i=1; i<=NbinsPt; i++){
        //Estraggo le projection binwise e faccio il fit
        std::cout << "print 5" << std::endl;
        std::cout << "print 6" << std::endl;
        std::cout << "invmass_vsPt pointer: " << invmass_vsPt << std::endl;
        TString projName = Form("%s_Ptbin_%d", invmass_vsPt->GetName(), i);
        TH1D* proj = invmass_vsPt->ProjectionY(projName, i, i);
        proj->Write();
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
        RooFitResult *fit_result = background.fitTo(data, RooFit::Range("range1,range2"), RooFit::Save(), RooFit::SumW2Error(true), RooFit::PrintLevel(-1));

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
        TCanvas *canvas_bg = new TCanvas(Form("canvas_bg_%d", i), Form("Fit for proj_bin_%d - Background Only", i), 1200, 1200);
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
        TPaveText *paveText_bg = new TPaveText(0.8, 0.75, 0.9, 0.9, "NDC");
        paveText_bg->AddText(Form("#chi^{2} = %.2f", chi2_bg));
        paveText_bg->AddText(Form("gauss #mu = %.4f +/- %.4f", gauss_mu_val, gauss_mu_err));
        paveText_bg->AddText(Form("gauss #sigma = %.4f +/- %.5f", gauss_sigma_val, gauss_sigma_err));
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
        canvas_bg->SaveAs(Form("PlotConID2022/FitVerificationsPt/Background/FitBackgroundOnly_proj_bin_%d.png", i));
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
        RooRealVar frac("frac", "fraction of background", 0.5, 0.0, 1.0);
        RooAddPdf model("model", "signal + background", RooArgList(crystal, background), RooArgList(frac));

        // Esegui il fit segnale + fondo
        fit_result = model.fitTo(data, RooFit::Range("range_full"), RooFit::Save(), RooFit::SumW2Error(true), RooFit::PrintLevel(-1));

        // Store the mean and sigma values in the TGraphErrors objects
        graph_mu_vsPt->SetPoint(i-1, bincenters[i-1], mu_cb.getVal());
        graph_mu_vsPt->SetPointError(i-1, binhalfwidths[i-1], mu_cb.getError());
        
        graph_sigma_vsPt->SetPoint(i-1, bincenters[i-1], sigma_cb.getVal());
        graph_sigma_vsPt->SetPointError(i-1, binhalfwidths[i-1], sigma_cb.getError());
        
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
        paveText->AddText(Form("#sigma = %.4f +/- %.5f", sigma_cb.getVal(), sigma_cb.getError()));
        //verifico che il fit converga
        /*if (fit_result->status() == 0) {
        paveText->AddText("Fit converged");
        } else {
        paveText->AddText("Fit did not converge");
        }*/

        frame->addObject(paveText);
        frame->Draw();
         // Aggiunta della legenda
        TLegend *legend = new TLegend(0.15, 0.72, 0.35, 0.9);
        //legend->SetHeader(histTitle.c_str(), "C");  // Aggiunge l'intestazione con il titolo dell'istogramma
        legend->AddEntry(frame->getObject(0), "Data", "P");  // Dati
        legend->AddEntry(frame->getObject(1), "Background Fit", "L");  // Fit del fondo
        legend->AddEntry(frame->getObject(2), "Signal Fit", "L");  // Fit del segnale
        legend->AddEntry(frame->getObject(3), "Combined Fit", "L");  // Fit combinato
        legend->Draw();  // Disegna la legenda
        c->SaveAs(Form("PlotConID2022/FitVerificationsPt/%s_fitbin_%d.png", invmass_vsPt->GetName(), i));
        delete c;


        //verifico se la differenza tra il centro della gaussiana e il centro della crystal ball è compatibile con 0.5
        compatibility_psi2s[i-1] = (gauss_mu_val - mu_cb.getVal() - 0.5) / sqrt(gauss_mu_err*gauss_mu_err + mu_cb.getError()*mu_cb.getError());

    }

    //stampo compatibilità con psi2s
    for(int i=0; i< NbinsPt; i++){
    std::cout << "compatibilità tra Deltamu e 0.5: Bin " << i+1 << "  valore:" << compatibility_psi2s[i] << std::endl;
    }
    file_param_mc->Close();
    fileprojections->Close();
    }











