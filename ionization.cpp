#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <bits/stdc++.h>
#include "physicalconst.h"
#include "mathstuff.h"
#include "functions.h"
#include <string>
using namespace std;

int main()
{

    vector<double> frequency;
    generate_n(back_inserter(frequency), num, Logspace<>(log10(v_ion), log10(v_max), num));
    vector<double> frequency_sig;
    generate_n(back_inserter(frequency_sig), num, Logspace<>(log10(v_ion), log10(v_max), num));
    vector<double> radius; //= linspace(R_in, 10e18, 100);


    double PHI_initial, PHI_initial_s;
   
    std::cout << PHI_initial << "\n";
    double n_e = 100;
    double n__H = 1.001e2;
    double n__S = n__H * 1.e-4;
    double eps = 1, dens, alph;
    vector<double> n_h_grid, n_e_grid;
    vector<double> n_h_i_grid, xi_h_grid, n_s_grid, n_s_i_grid, xi_s_grid;
    double ni_to_nh, n_h, n_i, n_s, n_s_i, ni_to_ns;
    double x = n_i / n_h;
    double xi_new, xi_old;
    vector<vector<double>> absorbtion_grid;
    vector<vector<double>> tau_grid;
    ofstream couth("ionizationh_s_2.txt");
    ofstream couts("ionizations_2.txt");
    ofstream coutsig("sigmas.txt");
    ofstream coutsigh("sigmah.txt");
    vector<double> n_H[numr];
    double sig[num];
    double sigS[num];
    for (int i = 0; i < num; i++)
    {
        sig[i] = sigma(frequency[i]);
       // coutsigh<<frequency[i]<<" "<<sig[i]<<"\n";
    }
    for (int i = 0; i < num; i++)
    {
        sigS[i] = sigmaS(frequency_sig[i]);
        //coutsig<<frequency_sig[i]<<" "<<sigS[i]<<"\n";
    }
//return 0;
    absorbtion_grid.resize(numr);
    tau_grid.resize(numr);
    for (size_t k = 0; k < numr; k++)
    {
        absorbtion_grid[k].resize(num);
        tau_grid[k].resize(num);
    }

    for (int i = 0; i < numr; i++)
    {

        if (i == 0)
        {
            radius.push_back(R_in);
            n_e_grid.push_back(n_e);
            for (int j = 0; j < num; j++)
            {

                absorbtion_grid[i][j] = 0.0;
                tau_grid[i][j] = 0.0;
            }
            n_h_grid.push_back(0);
            n_h_i_grid.push_back(0);
            xi_h_grid.push_back(0);
            PHI_initial = simpson_rule(v_ion, v_max, num, radius.back(), tau_grid, i, PHI);

            ni_to_nh = PHI_initial / (n_e_grid.back() * alpha_B());
            n_i = n__H * ni_to_nh / ((1 + ni_to_nh));
            n_h = n__H / (1 + ni_to_nh);

            PHI_initial_s = simpson_rule(v_ion, v_max, num, radius.back(), tau_grid, i, PHIS);
            ni_to_ns = PHI_initial_s / (n_e_grid.back() * alpha_B_Sulphur());
            n_s_i = n__S * ni_to_ns / ((1 + ni_to_ns));
            n_s = n__S / (1 + ni_to_ns);
            dens = sqrt((n_i + n_s_i) * n_e_grid.back());
            for (int j = 0; j < num; j++)
            {
                absorbtion_grid[i][j] = n_h * sig[j] + n_s * sigS[j];
                tau_grid[i][j] = 0.0;
                absorbtion_grid[i + 1][j] = absorbtion_grid[i][j];
                tau_grid[i + 1][j] = tau_grid[i][j];
            }
            continue;
        }
        else if (i == 1)
        {
            radius.push_back(radius.back() + dr);
            // for (int j = 0; j < num; j++)
            // {
            //       absorbtion_grid[i + 1][j] = absorbtion_grid[i][j];
            //     tau_grid[i][j] = absorbtion_grid[i][j] * (radius.back() - radius[radius.size() - 2]);
            //}
            n_e_grid.push_back(n_e_grid.back());
            PHI_initial = simpson_rule(v_ion, v_max, num, radius.back(), tau_grid, i, PHI);
            ni_to_nh = PHI_initial / (n_e_grid.back() * alpha_B());
            n_i = n__H * ni_to_nh / (1 + ni_to_nh);
            n_h = n__H / (1 + ni_to_nh);

            PHI_initial_s = simpson_rule(v_ion, v_max, num, radius.back(), tau_grid, i, PHIS);
            ni_to_ns = PHI_initial_s / (n_e_grid.back() * alpha_B_Sulphur());
            n_s_i = n__S * ni_to_ns / ((1 + ni_to_ns));
            n_s = n__S / (1 + ni_to_ns);
            dens = sqrt((n_i + n_s_i) * n_e_grid.back());
            // dens = sqrt(n_i * n_e_grid.back());
            for (int j = 0; j < num; j++)
            {
                absorbtion_grid[i][j] = n_h * sig[j] + n_s * sigS[j];
                tau_grid[i][j] = absorbtion_grid[i][j] * (radius.back() - radius[radius.size() - 2]);
            }

            continue;
        }

        radius.push_back(radius.back() + dr); 
        if (i < numr - 1)
        {
            for (int j = 0; j < num; j++)
            {
                absorbtion_grid[i][j] = absorbtion_grid[i - 1][j] + (absorbtion_grid[i - 1][j] - absorbtion_grid[i - 2][j]) * (radius.back() - radius[radius.size() - 2]) / (radius[radius.size() - 2] - radius[radius.size() - 3]);

                tau_grid[i][j] = tau_grid[i - 1][j] + 0.5 * (absorbtion_grid[i][j] + absorbtion_grid[i - 1][j]) * (radius.back() - radius[radius.size() - 2]);
            }
        }
        n_e_grid.push_back(n_e_grid.back());

       

        while (true)
        {
            PHI_initial = simpson_rule(v_ion, v_max, num, radius.back(), tau_grid, i, PHI);
            ni_to_nh = PHI_initial / (n_e_grid.back() * alpha_B());
            n_i = n__H * ni_to_nh / (1 + ni_to_nh);
            n_h = n__H / (1 + ni_to_nh);

            PHI_initial_s = simpson_rule(v_ion, v_max, num, radius.back(), tau_grid, i, PHIS);
            ni_to_ns = PHI_initial_s / (n_e_grid.back() * alpha_B_Sulphur());
            n_s_i = n__S * ni_to_ns / ((1 + ni_to_ns));
            n_s = n__S / (1 + ni_to_ns);
            dens = sqrt((n_i + n_s_i) * n_e_grid.back());

            
            for (int j = 0; j < num; j++)
            {
                absorbtion_grid[i][j] = n_h * sig[j] + n_s * sigS[j];
                tau_grid[i][j] = tau_grid[i - 1][j] + 0.5 * (absorbtion_grid[i][j] + absorbtion_grid[i - 1][j]) * (radius.back() - radius[radius.size() - 2]);
            }

            if (abs(dens / n_e_grid.back() - 1) > 1e-4)
            {
                n_e_grid.back() = dens;
                
                continue;
            }
            n_e_grid.push_back(dens);
            n_h_grid.push_back(n_h);
            n_h_i_grid.push_back(n_i);
            xi_h_grid.push_back(n_i / (n_i + n_h));

            n_s_grid.push_back(n_s);
            n_s_i_grid.push_back(n_s_i);
            xi_s_grid.push_back(n_s_i / (n_s_i + n_s));

            couth << radius.back() / 3.086e+18 << " " << xi_h_grid.back() << " " << xi_s_grid.back() << "\n";
            std::cout << absorbtion_grid[i][2] << " " << i << " " << radius.back() / 3.086e+18 << " " << tau_grid[i][2] << " " << xi_h_grid.back() << " " << n_h << "\n";
            break;
        }
        if (xi_h_grid.back() < 1e-10 || tau_grid[i][num - 1] > 100)
        {
            break;
        }
    }

    return 0;
}
