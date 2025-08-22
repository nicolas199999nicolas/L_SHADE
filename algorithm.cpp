#include "algorithm.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <functional>
#include <cmath>
#include "functions.h"

using namespace std;

vector<double> differential_evolution(
    int D, int NP, int G, double pb, double c, double minVal, double maxVal,
    function<double(const vector<double>&)> f
) {
    // Initialization phase
    int H = NP; // 固定 H 為 NP
    vector<vector<double>> P(NP); // Population
    vector<double> MCR(H, 0.5), MF(H, 0.5); // 記錄交叉率與縮放因子
    vector<vector<double>> Archive; // 歷史存檔
    int k = 0; // Index counter
    mt19937 gen(static_cast<unsigned>(time(nullptr)));
    uniform_real_distribution<> rand01(0.0, 1.0);
    uniform_int_distribution<> randH(0, H-1);
    uniform_int_distribution<> randNP(0, NP-1);
    
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
        //成功的CR和F
        vector<double> S_CR, S_F, F_delta; // 新增 F_delta 用於 weighted mean
        vector<vector<double>> newP(NP);
        vector<double> newFitness(NP);
        cout << "Generation " << NFE << ": "<< endl;
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
            double pmin = 2.0 / NP;
            double pmax = 0.2;
            double p_i = pmin + (pmax - pmin) * rand01(gen);

            // p-best selection
            vector<int> sortedIdx(NP);
            iota(sortedIdx.begin(), sortedIdx.end(), 0);
            sort(sortedIdx.begin(), sortedIdx.end(), [&](int a, int b){ return fitness[a] < fitness[b]; });//排序
            //NFE+=2;
            int num_p = max(2, static_cast<int>(NP * p_i)); // 至少 2 個
            int pBestIdx = sortedIdx[static_cast<int>(rand01(gen) * num_p)];//隨機選一個

            // 隨機選取 r1 from P excluding i
            int r1;
            do { r1 = randNP(gen); } while (r1 == i);

            // 選取 xr2 from P U A excluding i and r1
            vector<int> candidates;
            //把族群 P 中除了自己（i）和 r1 以外的所有個體索引都加入 candidates
            for (int c = 0; c < NP; ++c) {
                if (c != i && c != r1) candidates.push_back(c);
            }
            //把 Archive（歷史存檔）裡的所有個體也加入 candidates
            for (int a = 0; a < Archive.size(); ++a) {
                candidates.push_back(NP + a); // 偏移索引
            }
            // 隨機選取 xr2 from candidates
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
                // 邊界處理，依照論文做法
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
                    F_delta.push_back(abs(fitness[i] - newFitness[i])); // 記錄 delta f
                    S_CR.push_back(CR[i]); // 記錄 CR
                    S_F.push_back(F[i]);   // 記錄 F
                }
                P[i] = newP[i];
                fitness[i] = newFitness[i];
            }
        }
        /*
        // 控制 archive 大小
        while (Archive.size() > NP+1) {
            //cout<<Archive.size() << " > " << NP << ", removing oldest individual from archive." << endl;
            int eraseIdx = static_cast<int>(rand01(gen) * Archive.size());
            Archive.erase(Archive.begin() + eraseIdx);
        }
            */
        // 更新 MCR, MF (使用 weighted mean)
        if (!S_CR.empty() && !S_F.empty()) {
            // 判斷 max(S_CR) == 0
            double max_SCR = *max_element(S_CR.begin(), S_CR.end());
            if (max_SCR == 0.0) {
                MCR[k] = 0.0; // 這裡用 0.0 代表 ⊥
            } else {
                MCR[k] = meanWL(S_CR, F_delta);     
            }
            MF[k]  = meanWL(S_F, F_delta);
            k++;
            if (k > H) k = 1;
        } else {
            // S_CR 或 S_F 為空，沿用上一代
            int prev = (k - 1 + H) % H;
            MCR[k] = MCR[prev];
            MF[k]  = MF[prev];
        }
        
        //Calculate NG+1
        int delta_NG = N_G(NP, 4, MAX_NFE, NFE)- N_G(NP, 4, MAX_NFE, NFE);
        cout << "delta_NG: " << delta_NG << endl;
        if(delta_NG<0){
            // 1. 依 fitness 排序 population，並刪除最差的個體
            vector<size_t> idx(P.size());
            iota(idx.begin(), idx.end(), 0);
            sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return fitness[a] < fitness[b]; });

            int remove_count = abs(delta_NG);
            // 取出最差的 remove_count 個 index，並從大到小排序
            vector<size_t> remove_idx(idx.end() - remove_count, idx.end());
            sort(remove_idx.rbegin(), remove_idx.rend());
            //for (size_t j : remove_idx) {
                //P.erase(P.begin());
                //fitness.erase(fitness.begin());
            //}
        }

    }
    // 回傳最佳解
    int bestIdx = min_element(fitness.begin(), fitness.end()) - fitness.begin();
    return P[bestIdx];
}