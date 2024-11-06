#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>

// Function to load initial parameters for the Crystal Ball function from a file
void LoadCrystalBallInitialParams(const char* filePath, double& mean, double& sigma, double& alpha, double& n) {
    TFile* file = TFile::Open(filePath);
    if (file && !file->IsZombie()) {
        // Example assumes parameters are stored as TParameter<double> or similar
        // Adapt the variable names as needed based on your fit_results.root structure
        mean = ((TParameter<double>*)file->Get("mean"))->GetVal();
        sigma = ((TParameter<double>*)file->Get("sigma"))->GetVal();
        alpha = ((TParameter<double>*)file->Get("alpha"))->GetVal();
        n = ((TParameter<double>*)file->Get("n"))->GetVal();
        file->Close();
    } else {
        std::cerr << "Error loading parameters from file!" << std::endl;
    }
}

// Function to fit each histogram
TF1* FitHistogram(TH1D* hist, double mean, double sigma, double alpha, double n, double sidebandLow, double sidebandHigh) {
    // Fit the background in the sideband region with a linear function
    TF1* background = new TF1("background", "pol1", sidebandLow, sidebandHigh);
    hist->Fit(background, "R");  // Fit only in the sideband region

    // Combined fit function: Crystal Ball + background
    TF1* combinedFit = new TF1("combinedFit", "[0]*TMath::CrystalBall(x, [1], [2], [3], [4]) + pol1(5)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Set initial parameters for Crystal Ball
    combinedFit->SetParameter(0, hist->GetMaximum());  // Scale
    combinedFit->SetParameter(1, mean);  // Mean
    combinedFit->SetParameter(2, sigma);  // Sigma
    combinedFit->SetParameter(3, alpha);  // Alpha
    combinedFit->SetParameter(4, n);  // n

    // Set initial parameters for the background from the sideband fit
    combinedFit->SetParameter(5, background->GetParameter(0));
    combinedFit->SetParameter(6, background->GetParameter(1));

    // Perform the combined fit
    hist->Fit(combinedFit, "R");

    delete background;  // Clean up background fit
    return combinedFit;
}

// Main macro function to process the file and perform fits on all histograms
void FitData(const char* inputFilePath, const char* paramsFilePath, double sidebandLow, double sidebandHigh) {
    // Load initial Crystal Ball parameters
    double mean, sigma, alpha, n;
    LoadCrystalBallInitialParams(paramsFilePath, mean, sigma, alpha, n);

    // Open input ROOT file
    TFile* inputFile = TFile::Open(inputFilePath);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Failed to open input file: " << inputFilePath << std::endl;
        return;
    }

    // Create an output file to store fitted histograms
    TFile* outputFile = new TFile("fit_results_output.root", "RECREATE");

    // Loop over histograms in the input file
    TIter next(inputFile->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TH1D* hist = dynamic_cast<TH1D*>(key->ReadObj());
        if (!hist) continue;  // Skip non-histogram objects

        std::cout << "Processing histogram: " << hist->GetName() << std::endl;

        // Perform the fits
        TF1* fitResult = FitHistogram(hist, mean, sigma, alpha, n, sidebandLow, sidebandHigh);

        // Write the histogram and fit result to the output file
        outputFile->cd();
        hist->Write();
    }

    // Clean up
    outputFile->Close();
    inputFile->Close();
}

