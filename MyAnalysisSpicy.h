#ifndef MyAnalysisSpicy_h
#define MyAnalysisSpicy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <TH1F.h>
#include <TCanvas.h>

using namespace std;

const Int_t kMaxTracks = 100;  //dimensione massima degli array

class MyAnalysisSpicy {
public:
    TTree *fChain; //puntatore al tree o alla chain che sto analizzando  (dato)

    Int_t fCurrent; //numero del tree corrente in una TChain

    Int_t Reweighted = 0; //Verrà messa a 1 se è stato fatto il ripesamento su Pt, a 2 se è stato fatto anche su dxy

    //istogrammi da ripesare
    TH1F *h_invMass_ECAL_ele;
    TH1F *h_dxyEle;
    TH1F *h_Pt_JPsi;
    TH1F *h_Pt_ele1;
    TH1F *h_Pt_ele2;
    TH1F *h_nPV;
    TH1F *h_nPU;



    const Int_t ntupla; // =0 per i dati e =1 per il MC, così il codice sa su cosa sta loopando

                                       //leaf types
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    Char_t          nPV;
    Char_t           nPU;
    Int_t           runNumber;
    Int_t           eleID[3];
    Short_t         chargeEle[3];
    Float_t         ptEle[3];
    Float_t         etaEle[3];
    Float_t         phiEle[3];
    Float_t         dxyEle[3];
    Float_t         dzEle[3];
    Float_t         ipEle[3];
    Float_t         dxyErrEle[3];
    Float_t         dzErrEle[3];
    Float_t         ipErrEle[3];
    Float_t         etaSCEle[3];
    Float_t         phiSCEle[3];
    Float_t         R9Ele[3];
    Float_t         triggeringEle[3];
    Float_t         LowE_pfmvaIdEle[3]; //nuova ID
    Float_t         pfmvaIdEle[3]; //nuova ID n.2
    Float_t         pfmvaIdEle22[3]; //nuova ID retraining 2022
    Float_t         pfRelIsoEle[3]; //variabile di isolamento


    //energie
    Float_t         rawEnergySCEle[3];
    Float_t         esEnergySCEle[3];
    Float_t         energyEle[3];
    Float_t         energy_5x5SC[3];
    Float_t         energy_ECAL_ele[3];
    Float_t         energy_ECAL_pho[3];

    //masse invarianti varie 
    Float_t         BMass;
    Float_t         invMass;
    Float_t         invMass_5x5SC;
    Float_t         invMass_ECAL_ele;
    Float_t         invMass_ECAL_pho;
    Float_t         invMass_rawSC;
    Float_t         invMass_rawSC_esSC;

    //variabili gen matching
    vector<float>   *Gen_Pt;
    vector<float>   *Gen_Eta;
    vector<float>   *Gen_Phi;
    vector<float>   *Gen_motherPdgId;         //scommentare per mc prompt
    vector<float>   *Gen_E; 


    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
     
                            //lista dei branches

    //branch solo dati
    TBranch         *b_runNumber;
    //TBranch         *b_lumiBlock;

    //
    TBranch        *b_nPV;
    TBranch       *b_nPU;
    TBranch        *b_eleID;
    TBranch        *b_chargeEle;
    TBranch        *b_ptEle;      
    TBranch        *b_etaEle;      
    TBranch        *b_phiEle;      
    TBranch        *b_dxyEle;      
    TBranch        *b_dzEle;       
    TBranch        *b_ipEle;       
    TBranch        *b_dxyErrEle;   
    TBranch        *b_dzErrEle;    
    TBranch        *b_ipErrEle; 
    TBranch        *b_etaSCEle;    //!
    TBranch        *b_phiSCEle;
    TBranch        *b_R9Ele; 
    TBranch        *b_triggeringEle;
    TBranch        *b_LowE_pfmvaIdEle;
    TBranch        *b_pfmvaIdEle;
    TBranch        *b_pfmvaIdEle22;
    TBranch        *b_pfRelIsoEle;

    TBranch        *b_rawEnergySCEle; //!
   TBranch        *b_esEnergySCEle; //!
   TBranch        *b_energyEle;   //!
   TBranch        *b_energy_5x5SC; //!
   TBranch        *b_energy_ECAL_ele; //!
   TBranch        *b_energy_ECAL_pho; 

    TBranch        *b_BMass;       //!
   TBranch        *b_invMass;     //!
   TBranch        *b_invMass_5x5SC; //!
   TBranch        *b_invMass_ECAL_ele; //!
   TBranch        *b_invMass_ECAL_pho; //!
   TBranch        *b_invMass_rawSC; //!
   TBranch        *b_invMass_rawSC_esSC;

   TBranch        *b_Gen_Pt;      //!
   TBranch        *b_Gen_Eta;     //!
   TBranch        *b_Gen_Phi;     //!
   TBranch        *b_Gen_E; 
   TBranch        *b_Gen_motherPdgId; //scommentare per mc prompt

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    MyAnalysisSpicy(TTree *tree=0, Int_t n =0);
   virtual ~MyAnalysisSpicy();
   virtual Int_t    Cut(Long64_t entry);  //per filtrare
   virtual Int_t    GenMatching(Long64_t entry); // per fare il gen matching
   virtual Int_t GetElectronGenMatch(Long64_t entry, int electronIndex);
   virtual Int_t    GetEntry(Long64_t entry);  //per accedere ai dati di una specifica entry del tree
   virtual Long64_t LoadTree(Long64_t entry);  //carica una nuova entry dal TTree e ritorna l'indice di file locale in cui l'entry è memorizzata
   virtual void     Init(TTree *tree);  //associa i branch del tree alle variabili membro della classe
   virtual void     Loop(); //codice di analisi, itera su tutte le entry
   //virtual void     Reweight(); 
   virtual Bool_t   Notify(); //per gestire cambiamenti quando si passa da un file all'altro della chain
   virtual void     Show(Long64_t entry = -1); //visualizza i dati di una entry, quella corrente se non specificato

   virtual void ApplyCorrectionsVsPtandRun();
   virtual void MonteCarloReweighting();
   virtual void ApplySmearingCorrectionsFullEle();
   virtual void ApplyCorrectionsFullEle();

    //N.B entry = evento


};

#endif