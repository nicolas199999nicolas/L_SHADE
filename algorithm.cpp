    return ui;
}
// mutation/crossover operator forward declarations
std::vector<double> mutation_op(
    const std::vector<double>& xi,
    const std::vector<double>& xpbest,
    const std::vector<double>& xr1,
    const std::vector<double>& xr2,
    double Fi,
    double minVal,
    double maxVal
);
std::vector<double> crossover_op(
    const std::vector<double>& xi,
    const std::vector<double>& vi,
    double CRi,
    std::uniform_real_distribution<>& rand01,
    std::mt19937& gen
);
#include "algorithm.h"
#include "functions.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <cstdlib>

using namespace std;
vector<double> differential_evolution(
    int D, int NP, int Gmax, double pb, double c, double minVal, double maxVal,
    function<double(const vector<double>&)> f
) {
    // Initialization phase
    int H = NP; // 固定 H 為 NP
    vector<vector<double>> P(NP); // Population
    vector<double> MCR(H, 0.5), MF(H, 0.5); // 記錄交叉率與縮放因子
    vector<vector<double>> Archive; // 歷史存檔
    int k = 1; // Index counter 從1開始
    std::mt19937 gen(static_cast<unsigned>(time(nullptr)));
    std::uniform_real_distribution<> rand01(0.0, 1.0);
    std::uniform_int_distribution<> randH(1, H); // 1 to H
    std::uniform_int_distribution<> randNP(0, NP-1);
    
    // 初始化族群
    vector<double> fitness(NP);
    for (int i = 0; i < NP; ++i) {
        P[i] = generateRandomIndividual(D, minVal, maxVal, gen);
        fitness[i] = f(P[i]);
    }
    
    // Main loop
    int MAX_NFE = D*10000;
    int NFE = 0 ;
    int Ng = NP; // initial Ng = Ninit = NP
    int G = 1; // G=1
    while (G <= Gmax && NFE < MAX_NFE) {
        vector<double> CR(NP), F(NP);
        vector<double> S_CR, S_F, S_delta;
        vector<vector<double>> newP(NP);
        vector<double> newFitness(NP);
        cout << "Generation " << G << " (NFE: " << NFE << "): " << endl;
        // mutate & crossover
        for (int i = 0; i < NP; ++i) {
            int r = randH(gen); // 隨機選取歷史索引 1 to H
            
            //CRi
            normal_distribution<> randCR(MCR[r-1], 0.1); 
            if (MCR[r-1] == -1.0) {
                CR[i] = 0.0;
            } else {
                CR[i] = min(1.0, max(0.0, randCR(gen)));
            }

            // Fi
            cauchy_distribution<> randF(MF[r-1], 0.1);
            double Fi;
            do {
                Fi = randF(gen);
            } while (Fi <= 0.0);
            if (Fi > 1.0) Fi = 1.0;
            F[i] = Fi;

            // p_i uniform [pb, 0.2]
            double pmin = pb;
            double pmax = 0.2;
            double p_i = pmin + (pmax - pmin) * rand01(gen);

            int pBestIdx, r1;
            vector<double> xr2;
            choose_pbest_and_r1r2(i, NP, D, p_i, fitness, P, Archive, gen, randNP, rand01, pBestIdx, r1, xr2);

            // mutation
            vector<double> vi = mutation_op(P[i], P[pBestIdx], P[r1], xr2, F[i], minVal, maxVal);

            // crossover
            vector<double> ui = crossover_op(P[i], vi, CR[i], rand01, gen);
            newP[i] = ui;
            newFitness[i] = f(ui);
            NFE++;
        }
    }
    // 回傳最佳解
    int bestIdx = min_element(fitness.begin(), fitness.end()) - fitness.begin();
    return P[bestIdx];
}

// mutation operator
vector<double> mutation_op(
    const vector<double>& xi,
    const vector<double>& xpbest,
    const vector<double>& xr1,
    const vector<double>& xr2,
    double Fi,
    double minVal,
    double maxVal
) {
    int D = xi.size();
    vector<double> vi(D);
    for (int j = 0; j < D; ++j) {
        vi[j] = xi[j] + Fi * (xpbest[j] - xi[j]) + Fi * (xr1[j] - xr2[j]);
        if (vi[j] < minVal) vi[j] = (minVal + xi[j]) / 2;
        if (vi[j] > maxVal) vi[j] = (maxVal + xi[j]) / 2;
    }
    return vi;
}

// crossover operator
vector<double> crossover_op(
    const vector<double>& xi,
    const vector<double>& vi,
    double CRi,
    std::uniform_real_distribution<>& rand01,
    std::mt19937& gen
) {
    int D = xi.size();
    vector<double> ui(D);
    int jrand = static_cast<int>(rand01(gen) * D);
    for (int j = 0; j < D; ++j) {
        if (rand01(gen) < CRi || j == jrand) ui[j] = vi[j];
        else ui[j] = xi[j];
    }
    return ui;
}
        }
        // selection
        for (int i = 0; i < NP; ++i) {
            if (newFitness[i] <= fitness[i]) {
                if (newFitness[i] < fitness[i]) {
                    Archive.push_back(P[i]); // 加入 archive
                    S_delta.push_back(abs(fitness[i] - newFitness[i])); // 記錄 delta f
                    S_CR.push_back(CR[i]);
                    S_F.push_back(F[i]);
                }
                P[i] = newP[i];
                fitness[i] = newFitness[i];
            } else {
                // xi,G+1 = xi,G (do nothing)
            }
        }
        // 控制 archive 大小
        uniform_int_distribution<> randArchive(0, Archive.size() - 1);
        while (Archive.size() > NP) {       
            int eraseIdx = randArchive(gen);
            Archive.erase(Archive.begin() + eraseIdx);
        }
        // 更新 MCR, MF
        if (!S_CR.empty() && !S_F.empty()) {
            double max_SCR = *max_element(S_CR.begin(), S_CR.end());
            if (max_SCR == 0.0) {
                MCR[k-1] = -1.0; // ⊥, C++ 索引從0
            } else {
                MCR[k-1] = meanWL(S_CR, S_delta); // arithmetic
            }
            MF[k-1] = meanWL(S_F, S_delta); // lehmer
            k++;
            if (k > H) k = 1;
        }
        
        // Optional LPSR strategy
        int Ng_next = N_G(NP, 4, MAX_NFE, NFE); // Eq. (10)
        int delta_NG = Ng_next - Ng;
        //cout << "delta_NG: " << delta_NG << endl;
        //if (NG+1<NG)  此處疑似論文有誤 
        if (delta_NG < 0) {
            // 依 fitness 排序 population，並刪除最差的個體
            vector<size_t> idx(P.size());
            iota(idx.begin(), idx.end(), 0);
            sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return fitness[a] < fitness[b]; });

            int remove_count = abs(delta_NG);
            if (remove_count > Ng - 4) remove_count = Ng - 4; // Nmin=4
            vector<size_t> remove_idx(idx.end() - remove_count, idx.end());
            sort(remove_idx.rbegin(), remove_idx.rend());
            for (size_t r_idx : remove_idx) {
                P.erase(P.begin() + r_idx);
                fitness.erase(fitness.begin() + r_idx);
            }
            // 更新 Ng, NP
            Ng = P.size();
            NP = Ng;
            // 更新 randNP
            uniform_int_distribution<> new_randNP(0, NP-1);
            randNP = new_randNP;
            // Resize Archive to new |P|
            uniform_int_distribution<> randArchive(0, Archive.size() - 1);
            while (Archive.size() > NP) {
                
                int eraseIdx = randArchive(gen);
                Archive.erase(Archive.begin() + eraseIdx);
            }
        }
        G++; 
    }
    // 回傳最佳解
    int bestIdx = min_element(fitness.begin(), fitness.end()) - fitness.begin();
    return P[bestIdx];
}
// mutation
void choose_pbest_and_r1r2(
    int i, int NP, int D, double p_i,
    const vector<double>& fitness,
    const vector<vector<double>>& P,
    const vector<vector<double>>& Archive,
    mt19937& gen,
    uniform_int_distribution<>& randNP,
    uniform_real_distribution<>& rand01,
    int& pBestIdx, int& r1, vector<double>& xr2
) {
    vector<int> sortedIdx(NP);
    iota(sortedIdx.begin(), sortedIdx.end(), 0);
    sort(sortedIdx.begin(), sortedIdx.end(), [&](int a, int b){ return fitness[a] < fitness[b]; });
    int num_p = max(2, static_cast<int>(NP * p_i));
    pBestIdx = sortedIdx[static_cast<int>(rand01(gen) * num_p)];
    do { r1 = randNP(gen); } while (r1 == i);
    vector<int> candidates;
    for (int c = 0; c < NP; ++c) {
        if (c != i && c != r1) candidates.push_back(c);
    }
    for (int a = 0; a < static_cast<int>(Archive.size()); ++a) {
        candidates.push_back(NP + a);
    }
    int xr2_idx = candidates[static_cast<int>(rand01(gen) * candidates.size())];
    xr2.resize(D);
    if (xr2_idx < NP) {
        xr2 = P[xr2_idx];
    } else {
        xr2 = Archive[xr2_idx - NP];
    }
}