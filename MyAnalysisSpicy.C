#include "MyAnalysisSpicy.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <TH2D.h>
#include <TProfile.h>
#include <TColor.h>

///////////////////////////////////// Funzioni globali /////////////////////////////////
void SetPadMargins(TCanvas *canvas, int padIndex) {
    if (!canvas) {
        std::cerr << "Error: Canvas pointer is null!" << std::endl;
        return;
    }

    TPad *pad = dynamic_cast<TPad*>(canvas->GetPad(padIndex));
    if (!pad) {
        std::cerr << "Error: TPad at index " << padIndex << " not found!" << std::endl;
        return;
    }

    // Imposta i margini per il TPad specificato
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0.05);
    pad->SetTopMargin(0.1);
    pad->SetBottomMargin(0.15);
}


TH1F* Weights(TH1F *histo1, TH1F *histo2) {  // ripesa l'istogramma 1 (MC) sull'istogramma 2 (dati)

    /////////////////////Definisco i colori
    int bluCMS = TColor::GetColor("#5790fc");
    int rosaCMS = TColor::GetColor("#964a8b");
    int rossoCMS = TColor::GetColor("#e42536");
    int gialloCMS = TColor::GetColor("#f89c20");
    //Plot pre-ripesamento
    TCanvas *c_rew = new TCanvas("c_rew", "Pre and post reweighting", 1200, 600);
    c_rew->Divide(2,1);
    c_rew->cd(1);
    SetPadMargins(c_rew, 1);

    histo1->SetLineColor(kBlack);
    histo1->SetFillColorAlpha(bluCMS, 0.7); // Colore di riempimento
    histo1->SetFillStyle(1001); 
    histo1->SetStats(kFALSE);
    histo1->SetLineWidth(2);
    histo1->SetMaximum(0.6);

    histo2->SetLineColor(kBlack);
    histo2->SetFillColorAlpha(rossoCMS, 0.7); // Colore di riempimento
    histo2->SetFillStyle(1001); 
    histo2->SetStats(kFALSE);
    histo2->SetLineWidth(2);
    histo1->Draw("HISTO F");
    histo2->Draw("HISTO F SAME");
    TLegend *leg1 = new TLegend(0.65, 0.75, 0.95, 0.9);
    //leg1->SetTextSize(0.05);
    leg1->AddEntry(histo1, "MC (pre-reweighting)", "f");
    leg1->AddEntry(histo2, "Data (pre-reweighting)", "f");
    leg1->Draw();

    ////////////////////////////////////////////////////////////
    TH1F *histoRatio = (TH1F*)histo2->Clone("histoRatio");
    histoRatio->SetTitle("Ratio of histo2 / histo1");
    histoRatio->Divide(histo1);

    TH1F *histo1_bis = (TH1F*)histo1->Clone("Copy of MC");
    histo1_bis->Multiply(histoRatio);
    ////////////////////////////////////////////////////////////
    //Plot post ripesamento
    c_rew->cd(2);
    SetPadMargins(c_rew, 2);
    histo1_bis->SetLineColor(kBlack);
    histo2->SetMaximum(0.6);
    histo2->Draw("HISTO F");
    histo1_bis->Draw("HISTO F SAME");

    TLegend *leg2 = new TLegend(0.65, 0.75, 0.95, 0.9);
    //leg2->SetTextSize(0.05);
    leg2->AddEntry(histo1_bis, "MC (post-reweighting)", "f");
    leg2->AddEntry(histo2, "Data (post-reweighting)", "f");
    leg2->Draw();
    ////////////////////////////////////////////////////////////
    c_rew->SaveAs("Reweighting_closure.png");
    delete c_rew;
    return histoRatio;
}

//////////////////////////////////////////////////////////////////////////////////

// Costruttore
MyAnalysisSpicy::MyAnalysisSpicy(TTree *tree, Int_t n) : fChain(0), ntupla(n) {
    if (tree == 0) {
        return;
    }
    Init(tree);
    
    h_invMass_ECAL_ele = new TH1F("h_invMass_ECAL_ele", "Invariant Mass (ECAL ele);Mass [GeV];A.U.", 120, 0, 6);
    h_dxyEle = new TH1F("h_dxyEle", "Displacement on transverse plane [mm];dxy [mm];A.U.", 30, 0, 1.5);
    h_Pt_JPsi = new TH1F("h_Pt_JPsi", "Transverse momentum of J/#Psi; P_{T} [GeV]; A.U.", 30, 0, 60); //verrà riempito solo nella regione di segnale
}

  // n=0 per i dati, n=1 per il MC


// Distruttore
MyAnalysisSpicy::~MyAnalysisSpicy() {
    if (!fChain) return;
    //delete fChain->GetCurrentFile();  // Libera il file corrente della TChain
    delete h_invMass_ECAL_ele;
    delete h_dxyEle;
    delete h_Pt_JPsi;
}

// Carica una specifica entry dal TTree
Int_t MyAnalysisSpicy::GetEntry(Long64_t entry) {
    // Se fChain è nullo, ritorna 0
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

// Carica un albero specifico dalla TChain
Long64_t MyAnalysisSpicy::LoadTree(Long64_t entry) {
    // Se fChain è nullo, ritorna -5
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

// Inizializza la classe associando i branch alle variabili membro
void MyAnalysisSpicy::Init(TTree *tree) {
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

// Collega le variabili ai branch

    if(ntupla == 0){ //solo per i dati
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
    //fChain->SetBranchAddress("lumiBlock", lumiBlock, &b_lumiBlock);
    }

    fChain->SetBranchAddress("eleID", eleID, &b_eleID);
    fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
    fChain->SetBranchAddress("ptEle", ptEle, &b_ptEle);
    fChain->SetBranchAddress("etaEle", etaEle, &b_etaEle);
    fChain->SetBranchAddress("phiEle", phiEle, &b_phiEle);
    fChain->SetBranchAddress("dxyEle", dxyEle, &b_dxyEle);
    fChain->SetBranchAddress("dzEle", dzEle, &b_dzEle);
    fChain->SetBranchAddress("ipEle", ipEle, &b_ipEle);  //controllare con stampe /////////
    fChain->SetBranchAddress("dxyErrEle", dxyErrEle, &b_dxyErrEle);
    fChain->SetBranchAddress("dzErrEle", dzErrEle, &b_dzErrEle);
    fChain->SetBranchAddress("ipErrEle", ipErrEle, &b_ipErrEle);
    fChain->SetBranchAddress("etaSCEle", etaSCEle, &b_etaSCEle);
    fChain->SetBranchAddress("phiSCEle", phiSCEle, &b_phiSCEle);
    fChain->SetBranchAddress("R9Ele", R9Ele, &b_R9Ele);
    fChain->SetBranchAddress("triggeringEle", triggeringEle, &b_triggeringEle);
    fChain->SetBranchAddress("LowE_pfmvaIdEle", LowE_pfmvaIdEle, &b_LowE_pfmvaIdEle);
    fChain->SetBranchAddress("pfmvaIdEle", pfmvaIdEle, &b_pfmvaIdEle);
    fChain->SetBranchAddress("pfRelIsoEle", pfRelIsoEle, &b_pfRelIsoEle);
    fChain->SetBranchAddress("rawEnergySCEle", rawEnergySCEle, &b_rawEnergySCEle);
    fChain->SetBranchAddress("esEnergySCEle", esEnergySCEle, &b_esEnergySCEle);
    fChain->SetBranchAddress("energyEle", energyEle, &b_energyEle);
    fChain->SetBranchAddress("energy_5x5SC", energy_5x5SC, &b_energy_5x5SC);
    fChain->SetBranchAddress("energy_ECAL_ele", energy_ECAL_ele, &b_energy_ECAL_ele);
    fChain->SetBranchAddress("energy_ECAL_pho", energy_ECAL_pho, &b_energy_ECAL_pho);
    fChain->SetBranchAddress("BMass", &BMass, &b_BMass);
    fChain->SetBranchAddress("invMass", &invMass, &b_invMass);
    fChain->SetBranchAddress("invMass_5x5SC", &invMass_5x5SC, &b_invMass_5x5SC);
    fChain->SetBranchAddress("invMass_ECAL_ele", &invMass_ECAL_ele, &b_invMass_ECAL_ele);
    fChain->SetBranchAddress("invMass_ECAL_pho", &invMass_ECAL_pho, &b_invMass_ECAL_pho);
    fChain->SetBranchAddress("invMass_rawSC", &invMass_rawSC, &b_invMass_rawSC);
    fChain->SetBranchAddress("invMass_rawSC_esSC", &invMass_rawSC_esSC, &b_invMass_rawSC_esSC);
    
    /*fChain->SetBranchAddress("Gen_Pt", &Gen_Pt, &b_Gen_Pt);
    fChain->SetBranchAddress("Gen_Eta", &Gen_Eta, &b_Gen_Eta);
    fChain->SetBranchAddress("Gen_Phi", &Gen_Phi, &b_Gen_Phi);
    fChain->SetBranchAddress("Gen_E", &Gen_E, &b_Gen_E); */
    
    Notify();
}

// Notifica quando cambia il file nella TChain
Bool_t MyAnalysisSpicy::Notify() {
    return kTRUE;
}

// Mostra i dati di una entry specifica
void MyAnalysisSpicy::Show(Long64_t entry) {
    if (!fChain) return;
    // Se l'entry è -1, si usa l'entry corrente
    if (entry == -1) entry = fChain->GetReadEntry();
    // Accedi all'entry
    GetEntry(entry);
    // Stampa i valori delle variabili per quella entry
    printf("Entry %lld:\n", entry);
    printf("chargeEle[0] = %d, ptEle[0] = %f, etaEle[0] = %f\n", chargeEle[0], ptEle[0], etaEle[0]);
    // questo è un esempio, puoi stampare qualunque cosa

}

// Filtra gli eventi in base a qualche criterio
Int_t MyAnalysisSpicy::Cut(Long64_t entry) {
    //selection options: 0 = no cuts; 1 = L0 cut; 2 = LL cut; 3 = LM cut; 4 = MM cut; 5 = lowE_ID cut
    int selection = 5; 
    bool iscut = false;

    if(selection == 0){
        iscut = ptEle[0] > 2 && ptEle[1] > 2 && triggeringEle[0] ==1 &&  chargeEle[0]* chargeEle[1] < 0 && triggeringEle[1] ==1;
    }else if(selection ==1){
        iscut = ptEle[0] > 2 && ptEle[1] > 2 && triggeringEle[0] ==1.0 && triggeringEle[1] ==1.0 &&  chargeEle[0]* chargeEle[1] < 0 && (eleID[0] + eleID[1] >= 4) ; //L0
    }else if(selection ==2){
        iscut = ptEle[0] > 2 && ptEle[1] > 2 && triggeringEle[0] ==1.0 && triggeringEle[1] ==1.0 &&  chargeEle[0]* chargeEle[1] < 0 && (eleID[0]*eleID[1]>= 16) ; //LL
    }else if(selection ==3){
        iscut = ptEle[0] > 2 && ptEle[1] > 2 && triggeringEle[0] ==1.0 && triggeringEle[1] ==1.0 && chargeEle[0]* chargeEle[1] < 0 && (eleID[0]*eleID[1]>= 48) ; //LM
    }else if(selection == 4){
        iscut = ptEle[0] > 2 && ptEle[1] > 2 && triggeringEle[0] ==1.0 && triggeringEle[1] ==1.0 && eleID[0]>=12 && eleID[1] >=12; //MM
    }else if(selection == 5){
        iscut = ptEle[0] > 2 && ptEle[1] > 2 && triggeringEle[0] ==1.0 && triggeringEle[1] ==1.0 && LowE_pfmvaIdEle[0] > -1.0  && LowE_pfmvaIdEle[1] > -1.0 ;

    }

    if(iscut){
        return 1;
    }else{
        return 0;
    }

}


// Metodo principale dell'analisi, iterando su tutte le entry
void MyAnalysisSpicy::Loop() {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast(); //numero tot di eventi

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    //                                                          DICHIARAZIONE ISTOGRAMMI
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


    //ora faccio gli istogrammi
    TH1F *h_ptEle = new TH1F("h_ptEle", "P_t of Electron;p_{T} (GeV);A.U.", 60, 0, 40);
    TH1F *h_ptEle_2 = new TH1F("h_ptEle_2", "P_t of second Electron;p_{T} (GeV);A.U.", 60, 0, 40); //Pt del secondo elettrone
    TH1F *h_etaEle = new TH1F("h_etaEle", "Eta of Electron;#eta;A.U.", 30, -1.5, 1.5);
    TH1F *h_phiEle = new TH1F("h_phiEle", "Phi of Electron;#phi;A.U.", 30, -3.14, 3.14);

    TH1F *h_rawEnergySCEle = new TH1F("h_rawEnergySCEle", "Raw energy of Electron SuperCluster", 100, 0, 100);
    TH1F *h_energyEle = new TH1F("h_EnergyEle", " energyEle", 100, 0, 100);
    TH1F *h_energy_5x5SC = new TH1F("h_energy_5x5SC", " energy_5x5SC", 100, 0, 100);
    TH1F *h_energy_ECAL_ele = new TH1F("h_energy_ECAL_ele", "energy_ECAL_ele", 100, 0, 100);
    TH1F *h_energy_ECAL_pho = new TH1F("h_energy_ECAL_pho", "energy_ECAL_pho", 100, 0, 100);

    TH1F *h_BMass = new TH1F("h_BMass", "Invariant Mass of B;Mass [GeV];A.U.", 100, 0, 10);
    TH1F *h_invMass = new TH1F("h_invMass", "Invariant Mass with LM eleID cuts;Mass [GeV];A.U.", 100, 0, 6);
    TH1F *h_invMass_5x5SC = new TH1F("h_invMass_5x5SC", "Invariant Mass (5x5 Supercluster);Mass [GeV];A.U.", 100, 0, 6);
    //TH1F *h_invMass_ECAL_ele = new TH1F("h_invMass_ECAL_ele", "Invariant Mass (ECAL ele);Mass [GeV];A.U.", 60, 0, 6); //bin da 0.1
    TH1F *h_invMass_ECAL_pho = new TH1F("h_invMass_ECAL_pho", "Invariant Mass (ECAL pho);Mass [GeV];A.U.", 100, 0, 6);
    TH1F *h_invMass_rawSC = new TH1F("h_invMass_rawSC", "Invariant Mass (raw Supercluster);Mass [GeV];A.U.", 100, 0, 6);
    TH1F *h_invMass_rawSC_esSC = new TH1F("h_invMass_rawSC_esSC", "Invariant Mass (raw SC + ES SC);Mass [GeV];A.U.", 100, 0, 6);
    TH1F *h_my_invMass = new TH1F("h_my_invMass", "Invariant Mass reconstructed from kinematics; Mass [GeV];A.U.", 100, 0, 6);

    //TH1F *h_dxyEle = new TH1F("h_dxyEle", "Displacement on transverse plane [mm];dxy [mm];A.U.", 30, 0, 1.5); 
    TH1F *h_dzEle = new TH1F("h_dzEle", "Displacement along z [mm];dz [mm];A.U.", 30, -1.2, 1.2);  
    TH1F *h_ipEle = new TH1F("h_ipEle", "Impact parameter [mm];ip [mm];A.U.", 30, -1.2, 1.2); 

    TH1F *h_dxySignificance = new TH1F("h_dxySignificance", "Significance of dxy; ;A.U.", 40, -10, 10);   
    TH1F *h_dzSignificance = new TH1F("h_dzSignificance", "Significance of dz; ;A.U.", 40, -10, 10);   
    TH1F *h_ipSignificance = new TH1F("h_ipSignificance", "Significance of IP; ;A.U.", 40, -10, 10);

    //creo degli istogrammi con tagli diversi
    TH1F *h_invMass_nocuts = new TH1F("h_invMass_nocuts", "Invariant Mass without cuts;Mass [GeV];A.U.", 100, 0, 6);
    TH1F *h_invMass_LLcut = new TH1F("h_invMass_LLcut", "Invariant Mass with LL eleID cuts;Mass [GeV];A.U.", 100, 0, 6); // cut with only loose electrons (at least)
    TH1F *h_invMass_ECAL_ele_nocuts = new TH1F("h_invMass_ECAL_ele_nocuts", "ECAL invariant Mass without eleID cuts;Mass [GeV];A.U.", 100, 0, 6);
    TH1F *h_invMass_ECAL_ele_LLcut = new TH1F("h_invMass_ECAL_ele_LLcut", "ECAL invariant Mass with LL eleID cuts;Mass [GeV];A.U.", 50, 0, 6); // cut with only loose electrons (at least)
    TH1F *h_invMass_rawSC_nocuts = new TH1F("h_invMass_rawSC_nocuts", "Raw invariant Mass without cuts;Mass [GeV];A.U.", 100, 0, 6);

    //Istogrammi per le variabili della J_Psi con selezione
    TH1F *h_Beta_JPsi = new TH1F("h_Beta_JPsi", "#beta value of J/#psi; #beta; A.U.", 30, 0, 1);
    TH1F *h_BetaGamma_JPsi = new TH1F("h_BetaGamma_JPsi", "#beta#gamma value of J/#psi; #beta#gamma; A.U.", 30, 0, 20);
    TH1F *h_Rapidity_JPsi = new TH1F("h_Rapidity_JPsi", "Rapidity value of J/#psi; Rapidity; A.U.", 30, -5, 5);
    TH1F *h_phi_JPsi = new TH1F("h_phi_JPsi", "Phi of J/#psi; #phi; A.U.", 30, -3.14, 3.14);
    TH1F *h_eta_JPsi = new TH1F("h_eta_JPsi", "Pseudorapidity of J/#psi; Pseudorapidity; A.U.", 30, -5, 5);
    TH1F *h_Pt_JPsi_all = new TH1F("h_Pt_JPsi_all", "Pt of J/#psi ; p_{T}; A.U.", 30, 0 , 60); //quello dichiarato nel costruttore contiene solo regione di segnale

    //Istogrammi per le variabili della J_Psi senza selezione
    TH1F *h_Beta_JPsi_nocuts = new TH1F("h_Beta_JPsi_nocuts", "#beta value of J/#psi; #beta; A.U.", 30, 0, 1);
    TH1F *h_BetaGamma_JPsi_nocuts = new TH1F("h_BetaGamma_JPsi_nocuts", "#beta#gamma value of J/#psi; #beta#gamma; A.U.", 30, 0, 20);
    TH1F *h_Rapidity_JPsi_nocuts = new TH1F("h_Rapidity_JPsi_nocuts", "Rapidity value of J/#psi; Rapidity; A.U.", 30, -5, 5);
    TH1F *h_phi_JPsi_nocuts = new TH1F("h_phi_JPsi_nocuts", "Phi of J/#psi; #phi; A.U.", 30, -3.14, 3.14);
    TH1F *h_eta_JPsi_nocuts = new TH1F("h_eta_JPsi_nocuts", "Pseudorapidity of J/#psi; Pseudorapidity; A.U.", 30, -5, 5);
    TH1F *h_Pt_JPsi_all_nocuts = new TH1F("h_Pt_JPsi_all_nocuts", "p_{T} of J/#psi; p_{T}; A.U.", 30, 0 , 60);

    //Istogrammi differenza displacement
    TH1F *h_delta_dxy = new TH1F("h_delta_dxy", "|dxy[0] - dxy[1]|;|dxy[0] - dxy[1]|;A.U.", 40, 0, 0.2);
    TH1F *h_delta_dz = new TH1F("h_delta_dz", "|dz[0] - dz[1]|;|dz[0] - dz[1]|;A.U.", 40, 0, 0.3);
    TH1F *h_delta_ip = new TH1F("h_delta_ip", "|ip[0] - ip[1]|;|ip[0] - ip[1]|;A.U.", 40, 0, 0.4);

    //LowE ID nella regione di fondo e di segnale
    TH1F *h_LowE_pfmvaIdEle = new TH1F("h_LowE_pfmvaIdEle", "Multivariate ID; LowE_pfmvaIdEle; A.U.", 50, -12, 7); // fondo + segnale
    TH1F *h_LowE_pfmvaIdEle_sign = new TH1F("h_LowE_pfmvaIdEle_sign", "Multivariate ID in the Signal region; LowE_pfmvaIdEle; A.U.", 50, -12, 7);
    TH1F *h_LowE_pfmvaIdEle_back = new TH1F("h_LowE_pfmvaIdEle_back", "Multivariate ID in the Background region; LowE_pfmvaIdEle; A.U.", 50, -12, 7);

    //Fondo stimato dai dati
    TH1F *h_invmass_background = new TH1F("h_invmass_background", "data driven background;Invariant mass [Gev];A.U.", 100, 0, 6);

    //Bin personalizzati per l'asse di Pt
    double Ptbins[] = {4, 7, 8, 10, 14, 40};  
    int nbinsPt = sizeof(Ptbins)/sizeof(double) - 1; 
    TH2D *h_scale_vs_pt = new TH2D("h_scale_vs_pt", "TH2D of invariant mass and p_{T}[0];p_{T}[0];Invariant Mass (ECAL) [GeV]", nbinsPt, Ptbins, 150, 0, 6);
    TH2D *h_rawSC_vs_pt = new TH2D("h_rawSC_vs_pt", "TH2D of raw invMass and p_{T}[0]; p_{T}[0]; Raw SC invMass [GeV]", nbinsPt, Ptbins, 150, 0, 6);

    //Istogramma 2D di massa invariante vs run number
    TH2D *h_invmass_runNumber_eta0_06 = new TH2D("h_invmass_runNumber_eta0_06", "TH2D of invariant mass and run Number, 0 < eta < 0.6; runNumber ;Invariant Mass (ECAL) [GeV]", 15, 360000, 362500, 150, 2.6, 3.1);
    TH2D *h_invmass_runNumber_eta06_12 = new TH2D("h_invmass_runNumber_eta06_12", "TH2D of invariant mass and run Number, 0.6 < eta < 1.22; runNumber ;Invariant Mass (ECAL) [GeV]", 15, 360000, 362500, 150, 2.6, 3.1);
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    //                                                          LOOP SULLE ENTRIES (con riempimento istogrammi)
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break; //se è negativo c'è stato un problema con il caricamento del file
        fChain->GetEntry(jentry);

        // Filtro gli eventi se necessario, con cut

        //calcolo la massa invariante
        double my_invMass;
        my_invMass = sqrt(2*ptEle[0]*ptEle[1]*(cosh(etaEle[0]-etaEle[1]) - cos(phiEle[0] - phiEle[1])));

        //ricostruisco il quadrimpulso dei due elettroni e della JPsi
        TLorentzVector P_ele1, P_ele2, P_JPsi;
        P_ele1.SetPtEtaPhiE(ptEle[0],etaEle[0],phiEle[0],energy_ECAL_ele[0]);
        P_ele2.SetPtEtaPhiE(ptEle[1],etaEle[1],phiEle[1],energy_ECAL_ele[1]);
        P_JPsi = P_ele1 + P_ele2;
        double Beta_jpsi, BetaGamma_jpsi;



        // Riempio gli istogrammi con i valori delle variabili e normalizzo

        if(Cut(jentry)){
        h_ptEle->Fill(ptEle[0]);
        h_ptEle_2->Fill(ptEle[1]);
        h_etaEle->Fill(etaEle[0]);
        h_phiEle->Fill(phiEle[0]);
        
        h_rawEnergySCEle->Fill(rawEnergySCEle[0]);
        h_energyEle->Fill(energyEle[0]);
        h_energy_5x5SC->Fill(energy_5x5SC[0]);
        h_energy_ECAL_ele->Fill(energy_ECAL_ele[0]);
        h_energy_ECAL_pho->Fill(energy_ECAL_pho[0]);


        h_BMass->Fill(BMass);
        h_invMass->Fill(invMass);
        h_invMass_5x5SC->Fill(invMass_5x5SC);
        h_invMass_ECAL_ele->Fill(invMass_ECAL_ele);
        h_invMass_ECAL_pho->Fill(invMass_ECAL_pho);
        h_invMass_rawSC->Fill(invMass_rawSC);
        h_invMass_rawSC_esSC->Fill(invMass_rawSC_esSC);
        h_my_invMass->Fill(my_invMass);

        
        //riempio l'istogramma 2d solo se i due elettroni sono nello stesso bin di pt
        int bin_index_pt0, bin_index_pt1;
        bin_index_pt0 = h_scale_vs_pt->GetXaxis()->FindBin(ptEle[0]);
        bin_index_pt1 = h_scale_vs_pt->GetXaxis()->FindBin(ptEle[1]);
        if(bin_index_pt0 == bin_index_pt1){
        h_scale_vs_pt->Fill(ptEle[0],invMass_ECAL_ele);
        h_rawSC_vs_pt->Fill(ptEle[0], invMass_rawSC);

        //riempio l'istogramma massa invariante vs runNumber
        if(ntupla == 0 && invMass_ECAL_ele > 2.7 && invMass_ECAL_ele < 3.3){
            if(fabs(etaEle[0]) < 0.6 && fabs(etaEle[1]) < 0.6){
            h_invmass_runNumber_eta0_06->Fill(runNumber, invMass_ECAL_ele);  //Selezione spacca e pesa della regione di segnale
            }else if(fabs(etaEle[0]) > 0.6 && fabs(etaEle[1]) > 0.6){
            h_invmass_runNumber_eta06_12->Fill(runNumber, invMass_ECAL_ele);   
            }
        }
        }

        //variabili di displacement
        h_dxyEle->Fill(fabs(dxyEle[0]*10)); //passo in mm
        h_dzEle->Fill(dzEle[0]*10);
        h_ipEle->Fill(ipEle[0]*10);

        h_dxySignificance->Fill(dxyEle[0]/dxyErrEle[0]);
        h_dzSignificance->Fill(dzEle[0]/dxyEle[0]);
        h_ipSignificance->Fill(ipEle[0]/ipErrEle[0]);

        //spread nel displacement dei due elettroni
        h_delta_dxy->Fill(fabs(dxyEle[0] -  dxyEle[1]));
        h_delta_dz->Fill(fabs(dzEle[0] -  dzEle[1]));
        h_delta_ip->Fill(fabs(ipEle[0] -  ipEle[1]));


        //Riempio istogrammi nella regione di segnale

        if(invMass_ECAL_ele > 2.7 && invMass_ECAL_ele < 3.3){
             h_Pt_JPsi->Fill(P_JPsi.Pt());
        }

        //riempio gli istogrammi delle variabili della JPsi 

        h_Beta_JPsi->Fill(P_JPsi.P() / P_JPsi.E());
        h_BetaGamma_JPsi->Fill(P_JPsi.P() / P_JPsi.M());
        h_Rapidity_JPsi->Fill(P_JPsi.Rapidity());
        h_phi_JPsi->Fill(P_JPsi.Phi());
        h_eta_JPsi->Fill(P_JPsi.Eta());
        h_Pt_JPsi_all->Fill(P_JPsi.Pt());
       
        }

        //istogramma del fondo con !Cut

        if(chargeEle[0]*chargeEle[1] > 0 && ptEle[0] > 2 && ptEle[1] > 2 && triggeringEle[0] ==1.0 && triggeringEle[1] ==1.0 && (eleID[0]*eleID[1]>= 48) ){
            h_invmass_background->Fill(invMass_ECAL_ele);
        }

        h_invMass_nocuts->Fill(invMass);
        h_invMass_ECAL_ele_nocuts->Fill(invMass_ECAL_ele);
        h_invMass_rawSC_nocuts->Fill(invMass_rawSC);

        //Istogrammi delle variabili della JPsi senza nessun taglio
        h_Beta_JPsi_nocuts->Fill(P_JPsi.P() / P_JPsi.E());
        h_BetaGamma_JPsi_nocuts->Fill(P_JPsi.P() / P_JPsi.M());
        h_Rapidity_JPsi_nocuts->Fill(P_JPsi.Rapidity());
        h_phi_JPsi_nocuts->Fill(P_JPsi.Phi());
        h_eta_JPsi_nocuts->Fill(P_JPsi.Eta());
        h_Pt_JPsi_all_nocuts->Fill(P_JPsi.Pt());

        //LowE ID senza tagli
        h_LowE_pfmvaIdEle->Fill(LowE_pfmvaIdEle[0]);
        //LowE ID nella regione di segnale e nelle sidebands
        if(invMass_ECAL_ele > 2.6 && invMass_ECAL_ele < 3.1){
            h_LowE_pfmvaIdEle_sign->Fill(LowE_pfmvaIdEle[0]);
        }else{
            h_LowE_pfmvaIdEle_back->Fill(LowE_pfmvaIdEle[0]);
        }

       
    } //fine loop sulle entries


    h_invMass_ECAL_ele->Scale(1.0/ h_invMass_ECAL_ele->Integral());
    TProfile *prof = h_scale_vs_pt->ProfileX("mean inv Mass vs p_{T} without reweight");
    TH1D* proj_full = h_scale_vs_pt->ProjectionY("invmass_fullrange");

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    //                                                                 NORMALIZZAZIONE ISTOGRAMMI
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



   //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    //                                                                PLOT
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

    //istogramma 2d massa invariante vs runNumber
    if(ntupla == 0){
    /*TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    h_invmass_runNumber->SetStats(kFALSE);
    h_invmass_runNumber->Draw("COLZ");
    c1->SaveAs("invMass_runNumber.png");*/
    //Profilo della massa invariante vs il run number
    //primo bin su eta
    TProfile *prof_rNumber_1 = h_invmass_runNumber_eta0_06->ProfileX("prof_rNumber 0 < #eta < 0.6");  // Nome dell'oggetto
    prof_rNumber_1->SetTitle("Mean invariant mass vs run Number");  // Imposta il titolo
    prof_rNumber_1->GetXaxis()->SetTitle("Run Number");  // Titolo asse X per il profilo
    prof_rNumber_1->GetYaxis()->SetTitle("Mean Invariant Mass [GeV]");  // Titolo asse Y per il profilo
    prof_rNumber_1->GetYaxis()->SetRangeUser(2.85, 3);
    prof_rNumber_1->SetLineColor(kMagenta);        // Imposta il colore della linea (kBlue è un esempio)
    prof_rNumber_1->SetLineWidth(2);            // Imposta lo spessore della linea
    prof_rNumber_1->SetMarkerColor(kMagenta);       // Imposta il colore dei marker (kRed è un esempio)
    prof_rNumber_1->SetMarkerStyle(20);

    //secondo bin su eta
    TProfile *prof_rNumber_2 = h_invmass_runNumber_eta06_12->ProfileX("prof_rNumber 0.6 < #eta < 1.2");  // Nome dell'oggetto
    prof_rNumber_2->SetTitle("Mean invariant mass vs run Number");  // Imposta il titolo
    prof_rNumber_2->GetXaxis()->SetTitle("Run Number");  // Titolo asse X per il profilo
    prof_rNumber_2->GetYaxis()->SetTitle("Mean Invariant Mass [GeV]");  // Titolo asse Y per il profilo
    prof_rNumber_2->GetYaxis()->SetRangeUser(2.85, 3);
    prof_rNumber_2->SetLineColor(kGreen);        // Imposta il colore della linea (kBlue è un esempio)
    prof_rNumber_2->SetLineWidth(2);            // Imposta lo spessore della linea
    prof_rNumber_2->SetMarkerColor(kGreen);       // Imposta il colore dei marker (kRed è un esempio)
    prof_rNumber_2->SetMarkerStyle(20);

    TCanvas *c2 = new TCanvas("c2", "Canvas 2", 800, 600);
    prof_rNumber_1->SetStats(kFALSE);
    prof_rNumber_1->Draw("EP");  // Disegna il profilo con punti ed errori
    prof_rNumber_2->SetStats(kFALSE);
    prof_rNumber_2->Draw("EP SAME");
    // Aggiungi la legenda
    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);  // Definisce la posizione della legenda
    legend->AddEntry(prof_rNumber_1, "0 < #eta < 0.6", "l");  // Aggiunge la prima curva alla legenda
    legend->AddEntry(prof_rNumber_2, "0.6 < #eta < 1.2", "l");  // Aggiunge la seconda curva
    legend->Draw();
    c2->SaveAs("MeanInvmass_vs_runNumber.png");
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    //                                                          SCRITTURA SUI FILE
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

    //salvo gli istogrammi per poter plottare insieme dati e montecarlo
    TFile *outputFile;
    TFile *cutFile;
    if(ntupla == 0){
         outputFile = new TFile("outputHistograms_DATA_partF.root", "RECREATE");
         cutFile = new TFile("cutcomparison_DATA_partF.root",  "RECREATE");
    }else{
         outputFile = new TFile("outputHistograms_MC_LLcut.root", "RECREATE");
         cutFile = new TFile("cutcomparison_MC_LLcut.root",  "RECREATE");
    }
    outputFile->cd();
    //massa invariante (tracker, ECAL, Raw)
    h_invMass->Write();
    h_invMass_ECAL_ele->Write();
    h_invMass_rawSC->Write();
    //P_t dei due elettroni
    h_ptEle->Write();
    h_ptEle_2->Write();
    //variabili cinematiche non mass-related
    h_phiEle->Write();
    h_etaEle->Write();
    h_dxyEle->Write();
    h_dzEle->Write();
    h_ipEle->Write();
    //significanza del displacement
    h_dxySignificance->Write();
    h_dzSignificance->Write();
    h_ipSignificance->Write();

    //spread nel displacement dei due elettroni
    h_delta_dxy->Write();
    h_delta_dz->Write();
    h_delta_ip->Write();

    //Variabili Jpsi
    h_Beta_JPsi->Write();
    h_BetaGamma_JPsi->Write();
    h_Rapidity_JPsi->Write();
    h_phi_JPsi->Write();
    h_eta_JPsi->Write();
    h_Pt_JPsi_all->Write();

    //Transverse momentum della JPsi
    h_Pt_JPsi->Write();
    prof->Write();

    //Istogrammi per l'ID lowE
    h_LowE_pfmvaIdEle->Write();
    h_LowE_pfmvaIdEle_sign->Write();
    h_LowE_pfmvaIdEle_back->Write();

    //Projections della massa invariante (raw e regressed) nei vari bin di Pt
    for(int i=0; i < nbinsPt; i++){
    TH1D* proj = h_scale_vs_pt->ProjectionY(Form("proj_bin_%d", i+1), i+1, i+1);
    TH1D* proj_rawSC = h_rawSC_vs_pt->ProjectionY(Form("RawSC_bin_%d", i+1), i+1, i+1);
    proj->SetTitle(Form("%2.1lf GeV < p_{T} < %2.1lf GeV", h_scale_vs_pt->GetXaxis()->GetBinLowEdge(i+1), (h_scale_vs_pt->GetXaxis()->GetBinLowEdge(i+1) + h_scale_vs_pt->GetXaxis()->GetBinWidth(i+1))));
    //faccio la projection anche della massa raw SC
    proj_rawSC->SetTitle(Form("Raw Mass - %2.1lf GeV < p_{T} < %2.1lf GeV", h_rawSC_vs_pt->GetXaxis()->GetBinLowEdge(i+1), (h_rawSC_vs_pt->GetXaxis()->GetBinLowEdge(i+1) + h_rawSC_vs_pt->GetXaxis()->GetBinWidth(i+1))));
    proj->Write();
    proj_rawSC->Write();
    }
    proj_full->Write();
    h_invmass_background->Write();

    outputFile->Close();
    delete outputFile;
    //////////////////////////////////////////////////////
    cutFile->cd();  //salvo gli istogrammi con e senza i tagli per fare un confronto

    h_invMass->Write();
    h_invMass_ECAL_ele->Write();
    h_invMass_nocuts->Write();
    h_invMass_ECAL_ele_nocuts->Write();
    h_invMass_LLcut->Write();
    h_invMass_ECAL_ele_LLcut->Write();
    h_invMass_rawSC_nocuts->Write();
    h_invMass_rawSC->Write();

    //isto JPsi con tagli
    h_Beta_JPsi->Write();
    h_BetaGamma_JPsi->Write();
    h_Rapidity_JPsi->Write();
    h_phi_JPsi->Write();
    h_eta_JPsi->Write();
    h_Pt_JPsi_all->Write();

    //isto JPsi senza tagli
    h_Beta_JPsi_nocuts->Write();
    h_BetaGamma_JPsi_nocuts->Write();
    h_Rapidity_JPsi_nocuts->Write();
    h_phi_JPsi_nocuts->Write();
    h_eta_JPsi_nocuts->Write();
    h_Pt_JPsi_all_nocuts->Write();

    cutFile->Close();
    delete cutFile;

    //deallocare tutti gli istogrammi

    Reweighted = 0; //ecal invmass è stato riempito senza pesi

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void MyAnalysisSpicy::ReweightOnPt(){ 
    TFile *dataHistFile = TFile::Open("outputHistograms_DATA_partF.root", "UPDATE"); //outputHistograms_DATA_LLcut.root
    TH1F *h_ptjpsi_data = (TH1F*)dataHistFile->Get("h_Pt_JPsi");

    TFile *mcHistFile = TFile::Open("outputHistograms_MC_LLcut.root","UPDATE");
    TH1F *h_ptjpsi_mc = (TH1F*)mcHistFile->Get("h_Pt_JPsi");

    TH2D *h_scale_vs_pt = new TH2D("h_scale_vs_pt", "TH2D of invariant mass and P_{T}[0];P_{T}[0];Invariant Mass (ECAL)", 20, 0, 40, 150, 0, 6);
    //// Istogrammi massa invariante ecal vs track

    // Copio l'istogramma del MC che voglio ripesare
    TH1F* h_pt_mc_copy = (TH1F*)h_ptjpsi_mc->Clone("h_pt_mc_copy");

    h_pt_mc_copy->Scale(1.0/ h_pt_mc_copy->Integral());
    h_ptjpsi_data->Scale(1.0/h_ptjpsi_data->Integral());

    TH1F *histoWeights;
    histoWeights = Weights(h_pt_mc_copy, h_ptjpsi_data);
    h_ptjpsi_mc->Multiply(histoWeights); //posso usarlo per controllare che il Reweighting si chiuda

    //ora devo loopare sulle entries per applicare i pesi sia alla massa invariante che a dxy

    Long64_t nentries = fChain->GetEntriesFast(); 

    //svuoto gli istogrammi che dovrò ripesare
    h_invMass_ECAL_ele->Reset();
    h_dxyEle->Reset();

     for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break; //se è negativo c'è stato un problema con il caricamento del file
        fChain->GetEntry(jentry);
        if(Cut(jentry)){
        TLorentzVector P_ele1, P_ele2, P_JPsi;
        P_ele1.SetPtEtaPhiE(ptEle[0],etaEle[0],phiEle[0],energy_ECAL_ele[0]);
        P_ele2.SetPtEtaPhiE(ptEle[1],etaEle[1],phiEle[1],energy_ECAL_ele[1]);
        P_JPsi = P_ele1 + P_ele2;
        //Riempio l'istogramma con la massa ripesata e anche dxy post reweighting
    int binidx = h_Pt_JPsi->FindBin(P_JPsi.Pt());
    int nBins_w = histoWeights->GetNbinsX(); 
    if(binidx != 0 && binidx != nBins_w + 1) { // no underflow o overflow
    double weight = histoWeights->GetBinContent(binidx);

    int bin_index_pt0, bin_index_pt1;
    // Riempimento degli istogrammi con il valore pesato
    h_invMass_ECAL_ele->Fill(invMass_ECAL_ele, weight);
    h_dxyEle->Fill(dxyEle[0] * 10, weight);
    //riempio l'istogramma 2d solo se i due elettroni sono nello stesso bin di pt
    bin_index_pt0 = h_scale_vs_pt->GetXaxis()->FindBin(ptEle[0]);
    bin_index_pt1 = h_scale_vs_pt->GetXaxis()->FindBin(ptEle[1]);
    if(bin_index_pt0 == bin_index_pt1)h_scale_vs_pt->Fill(ptEle[0],invMass_ECAL_ele, weight);


            } 
        }  
    }
    
    h_dxyEle->Scale(1.0/ h_dxyEle->Integral());
    h_invMass_ECAL_ele->Scale(1.0/h_invMass_ECAL_ele->Integral());

    TProfile *prof = h_scale_vs_pt->ProfileX("mean inv Mass vs P_t after 1 reweight"); //voglio usarlo per vedere il profilo della massa invariante
    //aggiorno il contenuto degli istogrammi nel file .root e notifico che il ripesamento è stato effettuato
    mcHistFile->cd();
    TH1F *h_invMass_rew1 = (TH1F*)h_invMass_ECAL_ele->Clone("ECAL invMass reweighted on Pt_JPsi");
    h_invMass_rew1->Write();
    prof->Write("", TObject::kOverwrite);
    h_dxyEle->Write("", TObject::kOverwrite);
    Reweighted = 1;
    //chiudo i file che ho aperto
    mcHistFile->Close();
    dataHistFile->Close();
    delete mcHistFile;
    delete dataHistFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MyAnalysisSpicy::ReweightOnDxy(){
    
    TFile *dataHistFile = TFile::Open("outputHistograms_DATA_partF.root", "UPDATE");
    TH1F *h_dxy_data = (TH1F*)dataHistFile->Get("h_dxyEle");

    TFile *mcHistFile = TFile::Open("outputHistograms_MC_LLcut.root", "UPDATE");
    TH1F *h_dxy_mc = (TH1F*)mcHistFile->Get("h_dxyEle");

    //Bin personalizzati per l'asse di Pt
    double Ptbins[] = {4, 7, 8, 10, 14, 40};  
    int nbinsPt = sizeof(Ptbins)/sizeof(double) - 1; 
    TH2D *h_scale_vs_pt = new TH2D("h_scale_vs_pt", "TH2D of invariant mass and p_{T}[0];p_{T}[0];Invariant Mass (ECAL) [GeV]", nbinsPt, Ptbins, 150, 0, 6);
    TH2D *h_rawSC_vs_pt = new TH2D("h_rawSC_vs_pt", "TH2D of raw invMass and p_{T}[0]; p_{T}[0]; Raw SC invMass [GeV]", nbinsPt, Ptbins, 150, 0, 6);
    vector<TH1F*> histo_track_intervals; //massa invariante ricostruita con track info binnata in pt
    for(int i=0; i<nbinsPt ; i++){
        //std::string histName_ecal = "invmass_ecal_" + std::to_string(i);
        std::string histName_track = "invmass_track_" + std::to_string(i);
        //TH1F *hist_ecal = new TH1F(histName_ecal.c_str(), histName_ecal.c_str(), 100, 0, 6);
        TH1F *hist_track = new TH1F(histName_track.c_str(), histName_track.c_str(), 100, 0, 6);
        //histo_ecal_intervals.push_back(hist_ecal);
        histo_track_intervals.push_back(hist_track);
    }

    // Copio l'istogramma del MC che voglio ripesare
    TH1F* h_dxy_mc_copy = (TH1F*)h_dxy_mc->Clone("h_dxy_mc_copy");

    h_dxy_mc_copy->Scale(1.0/ h_dxy_mc_copy->Integral());
    h_dxy_data->Scale(1.0/h_dxy_data->Integral());

    TH1F *histoWeights;
    histoWeights = Weights(h_dxy_mc_copy, h_dxy_data);
    h_dxy_mc->Multiply(histoWeights); //posso usarlo per controllare che il Reweighting si chiuda

    //ora devo loopare sulle entries per applicare i pesi sia alla massa invariante che a dxy

    Long64_t nentries = fChain->GetEntriesFast(); 

    //svuoto gli istogrammi che dovrò ripesare
    h_invMass_ECAL_ele->Reset();

     for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break; //se è negativo c'è stato un problema con il caricamento del file
        fChain->GetEntry(jentry);
        if(Cut(jentry)){
        //Riempio l'istogramma con la massa ripesata e anche dxy post reweighting
    int binidx = h_dxyEle->FindBin(dxyEle[0]);
    int nBins_w = histoWeights->GetNbinsX(); 
    if(binidx != 0 && binidx != nBins_w + 1) { // no underflow o overflow
    double weight = histoWeights->GetBinContent(binidx);

    // Riempimento degli istogrammi con il valore pesato
    h_invMass_ECAL_ele->Fill(invMass_ECAL_ele, weight);


    //riempio l'istogramma 2d solo se i due elettroni sono nello stesso bin di pt
    int bin_index_pt0, bin_index_pt1;
    bin_index_pt0 = h_scale_vs_pt->GetXaxis()->FindBin(ptEle[0]);
    bin_index_pt1 = h_scale_vs_pt->GetXaxis()->FindBin(ptEle[1]);
    if(bin_index_pt0 == bin_index_pt1 && bin_index_pt0 != 0 && bin_index_pt0 != (nbinsPt+1)){
        //std::cout << ptEle[0] << " " << ptEle[1] << "bin ele0: "<< h_scale_vs_pt->GetXaxis()->FindBin(ptEle[0]) << "bin ele1: " << h_scale_vs_pt->GetXaxis()->FindBin(ptEle[0]) << std::endl;
        h_scale_vs_pt->Fill(ptEle[0],invMass_ECAL_ele, weight);
        //riempio gli istogrammi della massa invariante con track info
        histo_track_intervals[bin_index_pt0 -1]->Fill(invMass);
        //riempio anche l'istogramma con la massa raw SC
        h_rawSC_vs_pt->Fill(ptEle[0], invMass_rawSC, weight);
    }


            } 
        }  
    }

    h_invMass_ECAL_ele->Scale(1.0/h_invMass_ECAL_ele->Integral());
    TProfile *prof = h_scale_vs_pt->ProfileX("mean inv Mass vs p_{T} after 2 reweights"); //voglio usarlo per vedere il profilo della massa invariante

    //Intervalli in Pt su cui binnare la projection
    ////////////////////////////////////
    Int_t nBins = 5; 
    mcHistFile->cd();

   for(int i=0; i < nbinsPt; i++){
    TH1D* proj = h_scale_vs_pt->ProjectionY(Form("proj_bin_%d", i+1), i+1, i+1);
    TH1D* proj_rawSC = h_rawSC_vs_pt->ProjectionY(Form("rawSC_bin_%d", i+1), i+1, i+1);
    proj->SetTitle(Form("%2.1lf GeV < p_{T} < %2.1lf GeV", h_scale_vs_pt->GetXaxis()->GetBinLowEdge(i+1), (h_scale_vs_pt->GetXaxis()->GetBinLowEdge(i+1) + h_scale_vs_pt->GetXaxis()->GetBinWidth(i+1))));
    //faccio la projection anche della massa raw SC
    proj_rawSC->SetTitle(Form("Raw Mass - %2.1lf GeV < p_{T} < %2.1lf GeV", h_rawSC_vs_pt->GetXaxis()->GetBinLowEdge(i+1), (h_rawSC_vs_pt->GetXaxis()->GetBinLowEdge(i+1) + h_rawSC_vs_pt->GetXaxis()->GetBinWidth(i+1))));
    mcHistFile->Write("", TObject::kOverwrite);
    histo_track_intervals[i]->Write("", TObject::kOverwrite);
    
   }
   //histo_track_intervals[nbinsPt-1]->Write("", TObject::kOverwrite); //mancava l'ultimo
   TH1D* proj_full = h_scale_vs_pt->ProjectionY("invmass_fullrange");

    
    ///////////////////////////////////

    //aggiorno il contenuto degli istogrammi nel file .root e notifico che il ripesamento è stato effettuato
    TH1F *h_invMass_rew2 = (TH1F*)h_invMass_ECAL_ele->Clone("ECAL invMass reweighted on Pt_JPsi and dxy");
    h_invMass_rew2->Write("", TObject::kOverwrite);
    prof->Write("", TObject::kOverwrite);
    h_scale_vs_pt->Write("", TObject::kOverwrite);
    Reweighted = 2;
    //chiudo i file che ho aperto
    mcHistFile->Close();
    dataHistFile->Close();
    delete mcHistFile;
    delete dataHistFile;

}