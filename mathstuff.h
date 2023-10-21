#pragma once
#include"physicalconst.h"
#include<vector>
#include <algorithm>
#include <functional>
#include <iterator>
#include <math.h>
#include <cmath>
#include <numeric>

using namespace std;

template<typename T = double>
class Logspace {
private:
    T curValue, base, step;

public:
    Logspace(T first, T last, int num, T base = 10.0) : curValue(first), base(base){
       step = (last - first)/(num-1);
    }

    T operator()() {
        T retval = pow(base, curValue);
        curValue += step;
        return retval;
    }
};

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); 
                            
  return linspaced;
}
template <typename func_type>
double simpson_rule(double a, double b,
                    int n, double R_in, vector <vector <double>> & tau, int ir, // Number of intervals
                    func_type f)
{
    double h = (b - a) / n;

    // Internal sample points, there should be n - 1 of them
    double sum_odds = 0.0;
    for (int i = 1; i < n; i += 2)
    {
        sum_odds += f(R_in, tau[ir][i], a + i * h);
    }
    double sum_evens = 0.0;
    for (int i = 2; i < n; i += 2)
    {
        sum_evens += f(R_in,tau[ir][i],a + i * h);
    }

    return (f(R_in,tau[ir][0],a) + f(R_in,tau[ir][n-1],b) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}
template <typename func_type>
double simps(double a, double b,
                    int n,double n_h,double sig, // Number of intervals
                    func_type f)
{
    double h = (b - a) / n;

    // Internal sample points, there should be n - 1 of them
    double sum_odds = 0.0;
    for (int i = 1; i < n; i += 2)
    {
        sum_odds += f(n_h,sig,a + i * h);
    }
    double sum_evens = 0.0;
    for (int i = 2; i < n; i += 2)
    {
        sum_evens += f(n_h,sig,a + i * h);
    }

    return (f(n_h,sig,a) + f(n_h,sig,b) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}

