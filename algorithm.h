#include <vector>
#include <functional>
#include <random>
#include <ctime>
using namespace std;
// 傳入目標函數
vector<double> differential_evolution(
    int D, int NP, int G, double p, double c, double minVal, double maxVal,
    function<double(const vector<double>&)> func
);
void select_pbest_and_parents(
    int i, int NP, int D, double p_i,
    const vector<double>& fitness,
    const vector<vector<double>>& P,
    const vector<vector<double>>& Archive,
    mt19937& gen,
    uniform_int_distribution<>& randNP,
    uniform_real_distribution<>& rand01,
    int& pBestIdx, int& r1, vector<double>& xr2
);