#include <iostream>
#include <random>
#include <cmath>
#include "TCanvas.h"
#include "TH2D.h"

// Define the size of the lattice and the number of iterations
const int N = 3000;
const int iterations = 100000;

// Define the range of kT and H values to test
const int n_kT = 40;
const double kT_min = 1;
const double kT_max = 3;
const int n_H = 20;
const double H_min = -0.2;
const double H_max = 0.2;

// Define the lattice as a 2D array of random spins (+1 or -1)
std::vector<std::vector<int>> lattice(N, std::vector<int>(N));

// Define a random number generator
std::mt19937 rng;

// Define a function to calculate the energy change of flipping a spin at a given site
std::pair<double, double> energy(int i, int j, double H)
{
    double sum_neighbors = (double)(lattice[(i-1+N)%N][j]
                                + lattice[(i+1)%N][j]
                                + lattice[i][(j-1+N)%N]
                                + lattice[i][(j+1)%N])
                            + 0.5 * (double)(lattice[(i-1+N)%N][(j-1+N)%N]
                                            + lattice[(i+1)%N][(j-1+N)%N]
                                            + lattice[(i-1+N)%N][(j+1)%N]
                                            + lattice[(i+1)%N][(j+1)%N]);
    double E_plus = -sum_neighbors - H;
    double E_minus = sum_neighbors + H;
    return std::make_pair(E_plus, E_minus);
}

// Define a function to run the simulation for a single value of kT and H
double simulate(double kT, double H)
{
    // Initialize the lattice
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            lattice[i][j] = rng() % 2 == 0 ? 1 : -1;
        }
    }
    
    // Run the simulation until it reaches a steady state
    for (int n = 0; n < iterations; n++) {
        int i = rng() % N;
        int j = rng() % N;
        auto [E_plus, E_minus] = energy(i, j, H);
        double x = (E_plus - E_minus) / kT;
        double prob = 1.0 / (1.0 + std::exp(-x));
        if (rng() / double(rng.max()) < prob) {
            lattice[i][j] = 1;
        } else {
            lattice[i][j] = -1;
        }
    }
    
    // Calculate the ratio of black to white spins
    int count_white = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (lattice[i][j] == 1) {
                count_white++;
            }
        }
    }
    double ratio = double(count_white) / (N*N);
    return ratio;
}

void percolation_map()
{
    // Seed the random number generator
    std::random_device rd;
    rng.seed(rd());
    
    // Create a 2D histogram to store the results
    TH2D* h_ratios = new TH2D("ratios", "", n_kT, kT_min, kT_max, n_H, H_min, H_max);
    
    // Loop over each value of kT and H
    for (int i = 1; i <= n_kT; i++)
    {
        double kT = h_ratios->GetXaxis()->GetBinCenter(i);
        for(int j = 1; j <= n_H; j++)
        {
            double H = h_ratios->GetYaxis()->GetBinCenter(j);
            double ratio = simulate(kT, H);
            h_ratios->SetBinContent(i, j, ratio);
        }
    }

    // Create a canvas and draw the histogram as a 2D map
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    h_ratios->Draw("colz");
    h_ratios->GetXaxis()->SetTitle("kT");
    h_ratios->GetYaxis()->SetTitle("H");
    h_ratios->GetZaxis()->SetTitle("Ratio of White Spins");
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.12);
    canvas->SetBottomMargin(0.12);
    canvas->SetTopMargin(0.1);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->Update();

    // // Run the ROOT event loop
    // gSystem->Run();
}