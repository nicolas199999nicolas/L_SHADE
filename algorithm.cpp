#include "algorithm.h"
#include "functions.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <cstdlib>

using namespace std;

vector<double> differential_evolution(
    int D, int NP, int G, double pb, double c, double minVal, double maxVal,
    function<double(const vector<double>&)> f
) {
    // Initialization phase
    int H = NP; // 固定 H 為 NP, or set to a fixed value like 100
    vector<vector<double>> P(NP); // Population
    vector<double> MCR(H, 0.5), MF(H, 0.5); // 記錄交叉率與縮放因子
    vector<vector<double>> Archive; // 歷史存檔
    int k = 0; // Index counter (虛擬碼中 k=1, C++ 索引從0)
    std::mt19937 gen(static_cast<unsigned>(time(nullptr)));
    std::uniform_real_distribution<> rand01(0.0, 1.0);
    std::uniform_int_distribution<> randH(0, H-1);
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
    while ( NFE < MAX_NFE) {
        vector<double> CR(NP), F(NP);
        vector<double> S_CR, S_F, S_delta; // S_delta for weighted mean
        vector<vector<double>> newP(NP);
        vector<double> newFitness(NP);
        cout << "NFE " << NFE << ": "<< endl;
        // mutate & crossover
        for (int i = 0; i < NP; ++i) {
            int r = randH(gen); // 隨機選取歷史索引

            normal_distribution<> randCR(MCR[r], 0.1);
            cauchy_distribution<> randF(MF[r], 0.1);

            CR[i] = min(1.0, max(0.0, randCR(gen)));
            //當Fi<=0時，重新生成Fi直到它大於0
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

            // p-best selection
            vector<int> sortedIdx(NP);
            iota(sortedIdx.begin(), sortedIdx.end(), 0);
            sort(sortedIdx.begin(), sortedIdx.end(), [&](int a, int b){ return fitness[a] < fitness[b]; });
            int num_p = max(2, static_cast<int>(NP * p_i)); // 至少 2 個
            int pBestIdx = sortedIdx[static_cast<int>(rand01(gen) * num_p)];

            // 隨機選取 r1 from P excluding i
            int r1;
            do { r1 = randNP(gen); } while (r1 == i);

            // 選取 xr2 from P U A excluding i and r1
            vector<int> candidates;
            for (int c = 0; c < NP; ++c) {
                if (c != i && c != r1) candidates.push_back(c);
            }
            for (int a = 0; a < Archive.size(); ++a) {
                candidates.push_back(NP + a); // 偏移索引
            }
            int xr2_idx = candidates[static_cast<int>(rand01(gen) * candidates.size())];
            vector<double> xr2(D);
            if (xr2_idx < NP) {
                xr2 = P[xr2_idx];
            } else {
                xr2 = Archive[xr2_idx - NP];
            }

            // 變異
            vector<double> vi(D);
            for (int j = 0; j < D; ++j) {
                vi[j] = P[i][j] + F[i] * (P[pBestIdx][j] - P[i][j]) + F[i] * (P[r1][j] - xr2[j]);
                // 邊界處理
                if (vi[j] < minVal) vi[j] = (minVal + P[i][j]) / 2;
                if (vi[j] > maxVal) vi[j] = (maxVal + P[i][j]) / 2;
            }
                                
            // 交叉
            vector<double> ui(D);
            int jrand = static_cast<int>(rand01(gen) * D);
            for (int j = 0; j < D; ++j) {
                if (rand01(gen) < CR[i] || j == jrand) ui[j] = vi[j];
                else ui[j] = P[i][j];
            }
            newP[i] = ui;
            newFitness[i] = f(ui);
            NFE++;
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
            }
            
            /*
            if (newFitness[i] <= fitness[i]) {
                P[i] = newP[i];
                fitness[i] = newFitness[i];
            }else{
                //P[i+1] = P[i];
                //fitness[i+1] = fitness[i];
            }
            if (newFitness[i] < fitness[i]) {
                Archive.push_back(P[i]); // 加入 archive
                S_delta.push_back(abs(fitness[i] - newFitness[i])); // 記錄 delta f
                S_CR.push_back(CR[i]);
                S_F.push_back(F[i]);
            }
            */
        }
        // 控制 archive 大小
        while (Archive.size() > NP) {
            uniform_int_distribution<> randArchive(0, Archive.size() - 1);
            int eraseIdx = randArchive(gen);
            Archive.erase(Archive.begin() + eraseIdx);
        }
        // 更新 MCR, MF
        if (!S_CR.empty() && !S_F.empty()) {
            double max_SCR = *max_element(S_CR.begin(), S_CR.end());
            if (max_SCR == 0.0) {
                // 代表 ⊥, 虛擬碼中不更新, 這裡 skip or set to previous
                // 例如 MCR[k] = MCR[(k - 1 + H) % H];
            } else {
                MCR[k] = meanWL(S_CR, S_delta); // 平均
            }
            MF[k] = meanWL(S_F, S_delta); // 平均
            k = (k + 1) % H;
            if (k >= H) k = 1; // 虛擬碼中 If k > H, k=1
        }
        
        //Calculate NG+1
        int delta_NG = N_G(NP, 4, MAX_NFE, NFE) - NP; // 假設 N_G 返回 NG+1, delta = NG+1 - NP
        cout << "delta_NG: " << delta_NG << endl;//會跑出0 待解決
        if(delta_NG < 0){
            // 1. 依 fitness 排序 population，並刪除最差的個體
            vector<size_t> idx(P.size());
            iota(idx.begin(), idx.end(), 0);
            sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return fitness[a] < fitness[b]; });

            int remove_count = abs(delta_NG);
            if (remove_count > NP) remove_count = NP - 1; // 避免刪除過多
            // 取出最差的 remove_count 個 index，並從大到小排序
            vector<size_t> remove_idx(idx.end() - remove_count, idx.end());
            sort(remove_idx.rbegin(), remove_idx.rend());
            // 從大到小刪除，避免索引改變
            for (size_t r_idx : remove_idx) {
                P.erase(P.begin() + r_idx);
                fitness.erase(fitness.begin() + r_idx);
            }
            // 更新 NP
            NP = P.size();
            // 更新 randNP
            uniform_int_distribution<> new_randNP(0, NP-1);
            randNP = new_randNP;
        }
    }
    // 回傳最佳解
    int bestIdx = min_element(fitness.begin(), fitness.end()) - fitness.begin();
    return P[bestIdx];
}