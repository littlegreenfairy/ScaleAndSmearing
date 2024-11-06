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

void VerifyCorrections() {
    TFile *file = TFile::Open("outputHistograms_DATA_partF.root");

    // Crea un canvas
    TCanvas *c1 = new TCanvas("c1", "Profile Plots", 800, 600);

    // Recupera i TProfile dal file
    TProfile *invMass_vsrun_check = (TProfile*)file->Get("invMass_vsrun_check");
    TProfile *invMass_vsrun_corr = (TProfile*)file->Get("invMass_vsrun_corr");
    TProfile *invMass_vsrun_corr_offdiag = (TProfile*)file->Get("invMass_vsrun_corr_offdiag");

    // Disegna il primo profilo
    /*invMass_vsrun_check->SetLineColor(kRed);
    invMass_vsrun_check->SetMinimum(2.5);
    invMass_vsrun_check->SetMaximum(3.3);
    invMass_vsrun_check->SetTitle("Invariant Mass Profiles; Run Number; Invariant Mass");
    invMass_vsrun_check->Draw();*/

    // Disegna il secondo profilo in sovrapposizione
    invMass_vsrun_corr->SetLineColor(kBlue);
    invMass_vsrun_corr->SetStats(kFALSE);
    invMass_vsrun_corr->Draw();

    // Disegna il terzo profilo in sovrapposizione
    invMass_vsrun_corr_offdiag->SetLineColor(kGreen);
    invMass_vsrun_corr_offdiag->SetStats(kFALSE);
    invMass_vsrun_corr_offdiag->Draw("SAME");

    // Aggiungi la leggenda
    TLegend *legend = new TLegend(0.7, 0.3, 0.9, 0.5);
    legend->AddEntry(invMass_vsrun_check, "Invariant Mass Check (no correction)", "l");
    legend->AddEntry(invMass_vsrun_corr, "Invariant Mass Corrected", "l");
    legend->AddEntry(invMass_vsrun_corr_offdiag, "Invariant Mass Corrected Off-Diagonal", "l");
    legend->Draw();

    // Salva il canvas come immagine (opzionale)
     c1->SaveAs("correction_verification_vsrunN.png");

    // Chiudi il file
    file->Close();
}
