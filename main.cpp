#include <iostream>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <cmath>
#include <random>
#include <iomanip>
static std::mt19937 generator{std::random_device{}()};

std::vector<double> params_x;
std::vector<double> params_y;
std::vector<double> params_y_noise;
std::vector<double> params_y_filtered;

std::vector<double> best_params_y_filtered;
std::vector<double> best_weights;
std::vector<double> best_weights_5;
double min_J = 100;
double best_w = 100;
double best_d = 100;
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


double h_values[11];
double dis_values[11];
double weights_3_values[11][3];
double weights_5_values[11][5];
double w_values[11];
double d_values[11];
double J_values[11];

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
double chebishev_noise(){
    double max_noise = 0;
    double cur_noise;
    for(int k = 1; k <= K; k++){
        cur_noise = fabs(params_y_filtered[k] - params_y_filtered[k-1]);
        if(cur_noise > max_noise){
            max_noise = cur_noise;
        }
    }
    return max_noise;
}

double chebishev_difference(){
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

double max_distance(double w, double d){
    if(w > d){
        return w;
    }
    else if (d > w){
        return d;
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
    std::uniform_real_distribution<> middleweight{0, 1};

    int middle_index = 2;
    std::vector<double> weights; weights.resize(5);
    weights[2] = middleweight(generator);
    std::uniform_real_distribution<> middle{0, 1 - weights[2]};
    double mid = 0.5* middle(generator);
    weights[1] = mid;
    weights[3] = mid;

    double a = 0.5 * (1- weights[2] - weights[1]);
    weights[0] = a;
    weights[4] = a;

    return weights;
}
void harmonic_mean(std::vector<double> weights, int M_) {
    if(M_ == 1){
        params_y_filtered[0] = params_y_noise[0];
        params_y_filtered[K] = params_y_noise[K];
    } else{
        params_y_filtered[0] = params_y_noise[0];
        params_y_filtered[1] = params_y_noise[1];
        params_y_filtered[K] = params_y_noise[K];
        params_y_filtered[K-1] = params_y_noise[K-1];
    }
    double local_sum = 0;
    for (int k = M_; k <= K - M_; k++){
        local_sum = 0;
        for(int j = k - M_; j <= k + M_; j++){
            local_sum += weights[j-k+M_] / params_y_noise[j];
        }
        params_y_filtered[k] = 1 / local_sum;
    }

}



void PrintResult3() {
    std::cout << "+-----+--------+-------------------------------+--------+--------+\n"
              << "|  h  |   dis  |            weights            |    w   |    d   |\n"
              << "+-----+--------+-------------------------------+--------+--------+\n";
    for(int i = 0; i < 11; i++){
        std::cout << "|" << std::fixed <<std::setprecision(1) << std::setw(5)   << h_values[i]
                  << "|" << std::fixed <<std::setprecision(4) << std::setw(8)   << dis_values[i]
                  << "|" << std::fixed <<std::setprecision(4) << std::setw(31)  << "[" + std::to_string(weights_3_values[i][0]) +
                  ", " + std::to_string(weights_3_values[i][1])
                  + ", " + std::to_string(weights_3_values[i][2]) + "]"
                  << "|" << std::fixed <<std::setprecision(4) << std::setw(8)   << w_values[i]
                  << "|" << std::fixed <<std::setprecision(4) << std::setw(8)   << d_values[i]
                  << "|\n";
    }
    std::cout << "+-----+--------+-------------------------------+--------+--------+\n";

    size_t best_index= 0;
    double distance = 1;
    for (size_t it = 0; it < 11; it++){
        if(dis_values[it] < distance){
            distance = dis_values[it];
            best_index = it;
        }
    }

    std::cout << "+-----+---------+--------+--------+\n"
              << "|  h  |    J    |    w   |    d   |\n"
              << "+-----+---------+--------+--------+\n";

    std::cout << "|" << std::fixed <<std::setprecision(1) << std::setw(5)   << h_values[best_index]
    << "|" << std::fixed <<std::setprecision(4) << std::setw(9)   << J_values[best_index]
    << "|" << std::fixed <<std::setprecision(4) << std::setw(8)   << w_values[best_index]
    << "|" << std::fixed <<std::setprecision(4) << std::setw(8)   << d_values[best_index]
    << "|\n";

    std::cout << "+-----+---------+--------+--------+\n";

    std::vector<double> weig;
    weig.push_back(weights_3_values[best_index][0]);
    weig.push_back(weights_3_values[best_index][1]);
    weig.push_back(weights_3_values[best_index][2]);

    harmonic_mean(weig, 1);
    for (int k = 0; k <= K; k++){
        std::cout << "(" << params_x[k] << ";" << params_y_noise[k] << ")";
    }
    std::cout << "\n\n";
    for (int k = 0; k <= K; k++){
        std::cout << "(" << params_x[k] << ";" << params_y_filtered[k] << ")";
    }
    std::cout << "\n\n";
    for (int k = 0; k < 11; k++){
        std::cout << "(" << w_values[k] << "," << d_values[k] << ")";
    }
    std::cout << "\n\n";

}
void PrintResult5() {
    std::cout << "+-----+--------+--------------------------------------------------+--------+--------+\n"
              << "|  h  |   dis  |                  weights                         |    w   |    d   |\n"
              << "+-----+--------+--------------------------------------------------+--------+--------+\n";
    for (int i = 0; i < 11; i++) {
        std::cout << "|" << std::fixed << std::setprecision(1) << std::setw(5) << h_values[i]
                  << "|" << std::fixed << std::setprecision(4) << std::setw(8) << dis_values[i]
                  << "|" << std::fixed << std::setprecision(4) << std::setw(31)
                  << "[" + std::to_string(weights_5_values[i][0]) +
                     ", " + std::to_string(weights_5_values[i][1]) +
                     ", " + std::to_string(weights_5_values[i][2]) +
                     ", " + std::to_string(weights_5_values[i][3]) +
                     ", " + std::to_string(weights_5_values[i][4]) + "]"
                  << "|" << std::fixed << std::setprecision(4) << std::setw(8) << w_values[i]
                  << "|" << std::fixed << std::setprecision(4) << std::setw(8) << d_values[i]
                  << "|\n";
    }
    std::cout << "+-----+--------+--------------------------------------------------+--------+--------+\n";

    size_t best_index = 0;
    double distance = 1;
    for (size_t it = 0; it < 11; it++) {
        if (dis_values[it] < distance) {
            distance = dis_values[it];
            best_index = it;
        }
    }

    std::cout << "+-----+---------+--------+--------+\n"
              << "|  h  |    J    |    w   |    d   |\n"
              << "+-----+---------+--------+--------+\n";

    std::cout << "|" << std::fixed << std::setprecision(1) << std::setw(5) << h_values[best_index]
              << "|" << std::fixed << std::setprecision(4) << std::setw(9) << J_values[best_index]
              << "|" << std::fixed << std::setprecision(4) << std::setw(8) << w_values[best_index]
              << "|" << std::fixed << std::setprecision(4) << std::setw(8) << d_values[best_index]
              << "|\n";

    std::cout << "+-----+---------+--------+--------+\n";

    std::vector<double> weig;
    weig.push_back(weights_5_values[best_index][0]);
    weig.push_back(weights_5_values[best_index][1]);
    weig.push_back(weights_5_values[best_index][2]);
    weig.push_back(weights_5_values[best_index][3]);
    weig.push_back(weights_5_values[best_index][4]);

    harmonic_mean(weig, 2);
    for (int k = 0; k <= K; k++) {
        std::cout << "(" << params_x[k] << ";" << params_y_noise[k] << ")";
    }
    std::cout << "\n\n";
    for (int k = 0; k <= K; k++) {
        std::cout << "(" << params_x[k] << ";" << params_y_filtered[k] << ")";
    }
    std::cout << "\n\n";
    for (int k = 0; k < 11; k++) {
        std::cout << "(" << w_values[k] << "," << d_values[k] << ")";
    }
}
void RandomSearch(double l, int iterator, int R) {
    double N = (log(1-P))/log(1-(eps/(x_max - x_min)));
    double J = 0;
    for(int i = 0; i < N; i++) {
        std::vector<double> weights;
        if(R ==3) {
            weights = RandomWeights3();
            harmonic_mean(weights, 1);
        } else {
            weights = RandomWeights5();
            harmonic_mean(weights, 2);
        }


        double w = chebishev_noise();
        double d = chebishev_difference();
        J = l*w + (1 - l)*d;
        if(J < min_J){
            min_J=J;
            best_w = w;
            best_d = d;
            if(weights.size() == 3) {
                best_weights = weights;
            } else {
                best_weights_5 = weights;
            }
        }
    }
    h_values[iterator] = l;
    dis_values[iterator] = max_distance(best_w, best_d);

    if(R == 3) {
        weights_3_values[iterator][0] = best_weights[0];
        weights_3_values[iterator][1] = best_weights[1];
        weights_3_values[iterator][2] = best_weights[2];
    }
    else {
        weights_5_values[iterator][0] = best_weights_5[0];
        weights_5_values[iterator][1] = best_weights_5[1];
        weights_5_values[iterator][2] = best_weights_5[2];
        weights_5_values[iterator][3] = best_weights_5[3];
        weights_5_values[iterator][4] = best_weights_5[4];
    }
    w_values[iterator] = best_w;
    d_values[iterator] = best_d;
    J_values[iterator] = min_J;
    
    min_J = 100;
    best_w = 100;
    best_d = 100;

}

int main() {

    srand(time(nullptr));
    double xk;
    for(int j = 0; j <= K; j++){
        best_params_y_filtered.push_back(0);
    }
    for(int k = 0; k <= K; k++){
        xk = x_k(k);
        params_x.push_back(xk);
        params_y.push_back(signal(xk));
        params_y_noise.push_back(badsignal(xk));
        params_y_filtered.push_back(0);
    }
    int i = 0;
    for (double l = 0.0; l <=1; l+= 0.1) {

        RandomSearch(l, i, 3);
        i++;
    }
    i = 0;
    PrintResult3();

    for(int jter = 0; jter <= K; jter++){
        params_y_filtered[jter] = 0;
    }
    int j = 0;
    for (double l = 0.0; l <=1; l+= 0.1) {

        RandomSearch(l, j, 5);
        j++;
    }
    j = 0;
    PrintResult5();

    return 0;
}
