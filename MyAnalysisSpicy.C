#include "MyAnalysisSpicy.h"
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

    
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
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
    TH1F *h_etaEle = new TH1F("h_etaEle", "Eta of Electron;#eta;A.U.", 30, -1.5, 1.5);
    TH1F *h_phiEle = new TH1F("h_phiEle", "Phi of Electron;#phi;A.U.", 30, -3.14, 3.14);

    TH1F *h_invMass = new TH1F("h_invMass", "Invariant Mass with LM eleID cuts;Mass [GeV];A.U.", 100, 0, 6);
    TH1F *h_invMass_rawSC = new TH1F("h_invMass_rawSC", "Invariant Mass (raw Supercluster);Mass [GeV];A.U.", 100, 0, 6);
    TH1F *h_invMass_nocuts = new TH1F("h_invMass_nocuts", "Invariant Mass without cuts;Mass [GeV];A.U.", 100, 0, 6);

    TH1F *h_Pt_JPsi_all = new TH1F("h_Pt_JPsi_all", "Pt of J/#psi ; p_{T}; A.U.", 30, 0 , 60); //quello dichiarato nel costruttore contiene solo regione di segnale

    TH1F *h_Pt_JPsi_all_nocuts = new TH1F("h_Pt_JPsi_all_nocuts", "p_{T} of J/#psi; p_{T}; A.U.", 30, 0 , 60);

    //Istogrammi differenza displacement
    TH1F *h_delta_dxy = new TH1F("h_delta_dxy", "|dxy[0] - dxy[1]|;|dxy[0] - dxy[1]|;A.U.", 40, 0, 0.2);

    //LowE ID nella regione di fondo e di segnale
    TH1F *h_LowE_pfmvaIdEle = new TH1F("h_LowE_pfmvaIdEle", "Multivariate ID; LowE_pfmvaIdEle; A.U.", 50, -12, 7); // fondo + segnale
    TH1F *h_LowE_pfmvaIdEle_sign = new TH1F("h_LowE_pfmvaIdEle_sign", "Multivariate ID in the Signal region; LowE_pfmvaIdEle; A.U.", 50, -12, 7);
    TH1F *h_LowE_pfmvaIdEle_back = new TH1F("h_LowE_pfmvaIdEle_back", "Multivariate ID in the Background region; LowE_pfmvaIdEle; A.U.", 50, -12, 7);


    //Bin personalizzati per l'asse di Pt
    double Ptbins[] = {4, 7, 9, 11, 14, 20, 40}; 
    int nbinsPt = sizeof(Ptbins)/sizeof(double) - 1; 
    //Bin personalizzati per l'asse run Number
    int nbinsRunN = 5; //numero di bin (simmetrici) sul run Number
    double yMin = 360000;
    double yMax = 362500;
    double yBins[] = {360000, 360800, 361200, 361600, 362050, 362500};

    //Bin per la massa invariante
    int nbinsZ = 150;
    double zMin = 0;
    double zMax = 6;
    double zBins[nbinsZ + 1];
    double zStep = (zMax - zMin) / nbinsZ;
    for (int i = 0; i <= nbinsZ; ++i) {
        zBins[i] = zMin + i * zStep;
    }

    TH2D *h_scale_vs_pt = new TH2D("h_scale_vs_pt", "TH2D of invariant mass and p_{T}[0];p_{T}[0];Invariant Mass (ECAL) [GeV]", nbinsPt, Ptbins, 150, 0, 6);
    TH2D *h_rawSC_vs_pt = new TH2D("h_rawSC_vs_pt", "TH2D of raw invMass and p_{T}[0]; p_{T}[0]; Raw SC invMass [GeV]", nbinsPt, Ptbins, 150, 0, 6);

    TH3D *h_scalevs_pt_runN = new TH3D("h_scalevs_pt_runN", "Istogramma 3D invMass vs pT and run Number; p_{T}[0]; Run Number; m(e^{+}e^{-}) [GeV]",nbinsPt, Ptbins, nbinsRunN, yBins, nbinsZ, zBins);


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    //                                                          LOOP SULLE ENTRIES (con riempimento istogrammi)
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break; //se è negativo c'è stato un problema con il caricamento del file
        fChain->GetEntry(jentry);


        //ricostruisco il quadrimpulso dei due elettroni e della JPsi
        TLorentzVector P_ele1, P_ele2, P_JPsi;
        P_ele1.SetPtEtaPhiE(ptEle[0],etaEle[0],phiEle[0],energy_ECAL_ele[0]);
        P_ele2.SetPtEtaPhiE(ptEle[1],etaEle[1],phiEle[1],energy_ECAL_ele[1]);
        P_JPsi = P_ele1 + P_ele2;
        double Beta_jpsi, BetaGamma_jpsi;

        // Riempio gli istogrammi con i valori delle variabili e normalizzo

        if(Cut(jentry)){  //Applica il filtro sugli eventi
        h_ptEle->Fill(ptEle[0]);
        h_etaEle->Fill(etaEle[0]);
        h_phiEle->Fill(phiEle[0]);
        
        h_invMass->Fill(invMass);
        h_invMass_ECAL_ele->Fill(invMass_ECAL_ele);
        h_invMass_rawSC->Fill(invMass_rawSC);

        
        //riempio l'istogramma 2d solo se i due elettroni sono nello stesso bin di pt
        int bin_index_pt0, bin_index_pt1;
        bin_index_pt0 = h_scale_vs_pt->GetXaxis()->FindBin(ptEle[0]);
        bin_index_pt1 = h_scale_vs_pt->GetXaxis()->FindBin(ptEle[1]);
        if(bin_index_pt0 == bin_index_pt1){
        h_scale_vs_pt->Fill(ptEle[0],invMass_ECAL_ele); //qui siamo inclusivi nel run number
        h_rawSC_vs_pt->Fill(ptEle[0], invMass_rawSC);

        //riempio l'histo 3D, tanto i 2 elettroni dello stesso evento sono sempre nello stesso bin di run number
        h_scalevs_pt_runN->Fill(ptEle[0], runNumber, invMass_ECAL_ele);
        }

        //variabili di displacement
        h_dxyEle->Fill(fabs(dxyEle[0]*10)); //passo in mm

        //spread nel displacement dei due elettroni
        h_delta_dxy->Fill(fabs(dxyEle[0] -  dxyEle[1]));

        if(invMass_ECAL_ele > 2.7 && invMass_ECAL_ele < 3.3){
             h_Pt_JPsi->Fill(P_JPsi.Pt());
        }

        //riempio gli istogrammi delle variabili della JPsi 
        h_Pt_JPsi_all->Fill(P_JPsi.Pt());
       
        }

        h_invMass_nocuts->Fill(invMass);

        //Istogrammi delle variabili della JPsi senza nessun taglio
        h_Pt_JPsi_all_nocuts->Fill(P_JPsi.Pt());

        //LowE ID senza tagli
        h_LowE_pfmvaIdEle->Fill(LowE_pfmvaIdEle[0]);
        //LowE ID nella regione di segnale e nelle sidebands
        if(invMass_ECAL_ele > 2.7 && invMass_ECAL_ele < 3.3){
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
         outputFile = new TFile("outputHistograms_MC.root", "RECREATE");
         cutFile = new TFile("cutcomparison_MC.root",  "RECREATE");
    }
    outputFile->cd();
    //massa invariante (tracker, ECAL, Raw)
    h_invMass->Write();
    h_invMass_ECAL_ele->Write();
    h_invMass_rawSC->Write();
    //P_t dei due elettroni
    h_ptEle->Write();
    //variabili cinematiche non mass-related
    h_phiEle->Write();
    h_etaEle->Write();
    h_dxyEle->Write();

    //spread nel displacement dei due elettroni
    h_delta_dxy->Write();

    //Variabili Jpsi
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
    ////////////////////////////////////////////////////// ora binno sia in pt che sul Run number (Questa cosa solo nei dati!)
    if(ntupla == 0){
    for(int i=0; i < nbinsPt; i++){ //scorre sui bin di pt
        for(int j=0; j < nbinsRunN; j++){  //scorre sui bin della run number
          TH1D* proj_invmass = h_scalevs_pt_runN->ProjectionZ(Form("proj_bins_%d_%d", i+1, j+1), i+1, i+1, j+1, j+1);  //proietto sulla massa invariante binwise
          proj_invmass->SetTitle(Form("%2.1lf GeV < p_{T} < %2.1lf GeV - Run Number in [%f, %f]", h_scalevs_pt_runN->GetXaxis()->GetBinLowEdge(i+1), (h_scalevs_pt_runN->GetXaxis()->GetBinLowEdge(i+1) + h_scalevs_pt_runN->GetXaxis()->GetBinWidth(i+1)), h_scalevs_pt_runN->GetYaxis()->GetBinLowEdge(j+1), (h_scalevs_pt_runN->GetYaxis()->GetBinLowEdge(j+1) + h_scalevs_pt_runN->GetYaxis()->GetBinWidth(j+1))));
          proj_invmass->Write();
        }
    }
    }
    outputFile->Close();
    delete outputFile;
    //////////////////////////////////////////////////////
    cutFile->cd();  //salvo gli istogrammi con e senza i tagli per fare un confronto

    h_invMass->Write();
    h_invMass_ECAL_ele->Write();
    h_invMass_nocuts->Write();
    h_invMass_rawSC->Write();

    //isto JPsi con tagli
    h_Pt_JPsi_all->Write();

    //isto JPsi senza tagli
    h_Pt_JPsi_all_nocuts->Write();

    cutFile->Close();
    delete cutFile;

    //deallocare tutti gli istogrammi

    Reweighted = 0; //ecal invmass è stato riempito senza pesi

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void MyAnalysisSpicy::ReweightOnPt(){ 
    TFile *dataHistFile = TFile::Open("outputHistograms_DATA_partF.root", "UPDATE"); 
    TH1F *h_ptjpsi_data = (TH1F*)dataHistFile->Get("h_Pt_JPsi");

    TFile *mcHistFile = TFile::Open("outputHistograms_MC.root","UPDATE");
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

    Reweighted = 1;  //Notifica ripesamento
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

    TFile *mcHistFile = TFile::Open("outputHistograms_MC.root", "UPDATE");
    TH1F *h_dxy_mc = (TH1F*)mcHistFile->Get("h_dxyEle");

    //Bin personalizzati per l'asse di Pt
    double Ptbins[] = {4, 7, 9, 11, 14, 20, 40};  
    int nbinsPt = sizeof(Ptbins)/sizeof(double) - 1; 
    TH2D *h_scale_vs_pt = new TH2D("h_scale_vs_pt", "TH2D of invariant mass and p_{T}[0];p_{T}[0];Invariant Mass (ECAL) [GeV]", nbinsPt, Ptbins, 150, 0, 6);
    TH2D *h_rawSC_vs_pt = new TH2D("h_rawSC_vs_pt", "TH2D of raw invMass and p_{T}[0]; p_{T}[0]; Raw SC invMass [GeV]", nbinsPt, Ptbins, 150, 0, 6);

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
        h_scale_vs_pt->Fill(ptEle[0],invMass_ECAL_ele, weight);
        //riempio anche l'istogramma con la massa raw SC
        h_rawSC_vs_pt->Fill(ptEle[0], invMass_rawSC, weight);
    }


            } 
        }  
    }

    h_invMass_ECAL_ele->Scale(1.0/h_invMass_ECAL_ele->Integral());
    TProfile *prof = h_scale_vs_pt->ProfileX("mean inv Mass vs p_{T} after 2 reweights"); //voglio usarlo per vedere il profilo della massa invariante

    ////////////////////////////////////
    mcHistFile->cd();

   for(int i=0; i < nbinsPt; i++){
    TH1D* proj = h_scale_vs_pt->ProjectionY(Form("proj_bin_%d", i+1), i+1, i+1);
    TH1D* proj_rawSC = h_rawSC_vs_pt->ProjectionY(Form("rawSC_bin_%d", i+1), i+1, i+1);
    proj->SetTitle(Form("%2.1lf GeV < p_{T} < %2.1lf GeV", h_scale_vs_pt->GetXaxis()->GetBinLowEdge(i+1), (h_scale_vs_pt->GetXaxis()->GetBinLowEdge(i+1) + h_scale_vs_pt->GetXaxis()->GetBinWidth(i+1))));
    //faccio la projection anche della massa raw SC
    proj_rawSC->SetTitle(Form("Raw Mass - %2.1lf GeV < p_{T} < %2.1lf GeV", h_rawSC_vs_pt->GetXaxis()->GetBinLowEdge(i+1), (h_rawSC_vs_pt->GetXaxis()->GetBinLowEdge(i+1) + h_rawSC_vs_pt->GetXaxis()->GetBinWidth(i+1))));
    mcHistFile->Write("", TObject::kOverwrite);
    
   }
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