#include <iostream>
#include <vector>
#include "algorithm.h"
#include "functions.h"
#include "cec14_test_func.h"
using namespace std;

int main(int argc, char *argv[]) {
    // 選擇目標函數
    //auto targetFunction = F2; // 只需更改這一行即可切換目標函數
    
    // 1. 預計的測試函數列表
    // 包裝 cec14_wrapper 以符合 std::function<double(const vector<double>&)>
    auto cec14_func = [](const vector<double>& x) {
        // 這裡以 func_num=1 為例，可根據需求調整
        return cec14_wrapper(x, 5);
    };
    vector<function<double(const vector<double>&)>> funcs = {
        F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, cec14_func
    };
    
    // 2. 對應的 G (gen) 數值
    vector<int> gens = {
        1500, 2000, 5000, 5000, 3000, 100, 3000, 1000, 1000, 500, 500, 500, 500 , 1000
    };
    
    
    // 參數設置
    const int D = 50; // 維度
    const int NP = 100; // 種群大小
    const int G = 500; // 迭代次數
    const double pb = 0.05; //取前幾%的個體
    const double c = 0.1; //自適應參數
    

    cout << "Initializing parameters:\n";
    cout << "D: " << D << ", NP: " << NP << ", G: " << G << ", p: " << pb << ", c: " << c << "\n";

    
    //指定開頭和結尾函式
    //14=cec14_wrapper
    int st = 14  ,ed = 14;
    min(int(funcs.size()),ed);
    // 3. 每個函數的邊界最大值 F1 ~ F13
    vector<double> bounds = {
        100, 10, 100, 100, 30, 100, 1.28, 500, 5.12, 32, 600, 50, 50, 32
    };

    for (int i = st; i <= ed; ++i) {
        double bound = bounds[i-1];
        double minVal = -bound;
        double maxVal = bound;
        //cout << "==== F" << (i) << " (Gen=" << gens[i-1] << endl;
        vector<double> best = differential_evolution(D, NP, gens[i-1], pb, c, minVal, maxVal, funcs[i-1]);
        cout << "Best solution: ";
        for (double val : best) cout << val << " ";
        cout << "\nFitness: " << funcs[i-1](best) << "\n";
        cout<<"F"<<i<<" result:"<<"\n";
        cout << "========================\n";
    }
    
    system("pause");
    return 0;
}