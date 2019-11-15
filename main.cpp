#include <iostream>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <cmath>
#include <random>
static std::mt19937 generator{std::random_device{}()};

std::vector<double> params_x;
std::vector<double> params_y;
std::vector<double> params_y_noise;
std::vector<double> params_y_filtered;
//params
double x_min = 0;
double x_max = M_PI;
int K = 100;
double a = 0.25;
//k = [0,K]
double P = 0.95;
int L = 10;
//l = [0,L]
int r1 = 3, r2 = 5;
double eps = 0.01;
double bestX, bestY;

double lyamda(int l){
    return l/L;
}

int M(int r){
    return (r-1)/2;
}

double x_k(int k){
    return x_min + k*(x_max - x_min) / K;
}

double signal(double xk){
    return sin(xk) + 0.5;
}

double badsignal(double xk){
    std::uniform_real_distribution<> num{-a, a};
    double bad = num(generator);
    return signal(xk) + bad;
}

double chebishev(){
    double max_distance = 0;
    double cur_distance;
    for(int k = 0; k <= K; k++){
        cur_distance = fabs(params_y_filtered[k] - params_y_noise[k]);
        if(cur_distance > max_distance){
            max_distance = cur_distance;
        }
    }
    return max_distance;
}

void harmonic_mean(std::vector<double> weights) {
    if(weights.size() == 3){
        params_y_filtered[0] = params_y_noise[0];
        params_y_filtered[K] = params_y_noise[K];
        for (int k = 1; k <= 99; k++){
            double local_sum = 0;
            for(int j = k - 1, key = 0; j <= k +1, key<3; j++, key++){
                local_sum += weights[key] / params_y_noise[j];
            }
            local_sum = 1/ local_sum;
            params_y_filtered[k] = local_sum;
        }
    }
    else {
        params_y_filtered[0] = params_y_noise[0];
        params_y_filtered[1] = params_y_noise[1];
        params_y_filtered[K] = params_y_noise[K];
        params_y_filtered[K-1] = params_y_noise[K-1];
        for (int k = 2; k<= 98; k++){
            double local_sum = 0;
            for(int j = k-2, key = 0; j <= k + 2, key < 5; j++, key++){
                local_sum += weights[key] / params_y_noise[j];
            }
            local_sum = 1/ local_sum;
            params_y_filtered[k] = local_sum;
        }
    }
}

std::vector<double> RandomWeights3(){
    std::uniform_real_distribution<> number{0, 1};

    std::vector<double> weights; weights.resize(3);
    weights[1] = number(generator);

    double a = 0.5 * (1- weights[1]);
    weights[0] = a;
    weights[2] = a;

    return weights;
}

std::vector<double> RandomWeights5(){
    std::uniform_real_distribution<> number{0, 1};

    int middle_index = 2;
    std::vector<double> weights; weights.resize(5);
    weights[2] = number(generator);
    std::uniform_real_distribution<> middle{0, 1 - weights[2]};
    double mid = 0.5* middle(generator);
    weights[1] = mid;
    weights[3] = mid;

    double a = 0.5 * (1- weights[1]);
    weights[0] = a;
    weights[4] = a;

    return weights;
}


void RandomSearch() {
    double N = (log(1-P))/log(1-(eps/(x_max - x_min)));
    double best_y = 0, best_x = 0;
    std::uniform_real_distribution<> dist{x_min, x_max};
    std::vector<double> current_x_storage;
    for (int i = 0; i < N; i++){
        current_x_storage.push_back(dist(generator));
    }
    for (auto index : current_x_storage) {
        if (signal(index) > best_y) {
            best_y = signal(index);
            best_x = index;
        }
    }

    bestX = best_x; bestY = best_y;
}

int main() {

    double xk;
    for(int k = 0; k <= K; k++){
        xk = x_k(k);
        params_x.push_back(xk);
        params_y.push_back(signal(xk));
        params_y_noise.push_back(badsignal(xk));
        params_y_filtered.push_back(0);
    }
    std::vector<double> weights = RandomWeights3();
    std::cout << weights[0] << " " << weights[1] << " " << weights[2] << "\n\n";
    harmonic_mean(weights);
    for (int k = 0; k <= K; k++){
        std::cout << "(" << params_x[k] << ";" << params_y_noise[k] << ")";
    }
    std::cout << "\n\n";
    for (int k = 0; k <= K; k++){
        std::cout << "(" << params_x[k] << ";" << params_y_filtered[k] << ")";
    }
    return 0;
}