#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <fstream>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <vector>


struct Vec3 { 
    double x, y, z; 
    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
};

class NBodySystem {
    std::vector<double> masses;
    std::vector<Vec3> positions, velocities, accelerations;
    const double G = 6.67e-11;
    const double epsilon = 1e-3; // шаг 
    const int num_threads;

public:
    NBodySystem(const std::string& filename, int threads = 1) : num_threads(threads) {
        std::ifstream file(filename);
        int N;
        file >> N;
        
        masses.resize(N);
        positions.resize(N);
        velocities.resize(N);
        accelerations.resize(N);
        
        for (int i = 0; i < N; ++i) {
            file >> masses[i] 
                 >> positions[i].x >> positions[i].y >> positions[i].z
                 >> velocities[i].x >> velocities[i].y >> velocities[i].z;
        }
    }

    void computeStep(double tau) {
        const int N = positions.size();
        omp_set_num_threads(num_threads);
        
        // обнуление ускорений
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            accelerations[i] = Vec3(0, 0, 0);
        }

        // расчет гравитационных ускорений
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                
                double dx = positions[j].x - positions[i].x;
                double dy = positions[j].y - positions[i].y;
                double dz = positions[j].z - positions[i].z;
                
                double dist_sq = dx*dx + dy*dy + dz*dz;
                double dist = sqrt(dist_sq + epsilon*epsilon);
                double dist_cubed = dist * dist * dist;
                
                double force_mag = G * masses[j] / dist_cubed;
                
                accelerations[i].x += force_mag * dx;
                accelerations[i].y += force_mag * dy;
                accelerations[i].z += force_mag * dz;
            }
        }

        // интегрирование методом Эйлера
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            velocities[i].x += accelerations[i].x * tau;
            velocities[i].y += accelerations[i].y * tau;
            velocities[i].z += accelerations[i].z * tau;
            
            positions[i].x += velocities[i].x * tau;
            positions[i].y += velocities[i].y * tau;
            positions[i].z += velocities[i].z * tau;
        }
    }
    
    void clean_up()
    {
    for (int i = 0; i < positions.size(); ++i) 
       std::ofstream file("traj" + std::to_string(i+1) + ".txt", std::ios::trunc);
    }

    void saveTrajectories(double t) {
        for (int i = 0; i < positions.size(); ++i) {
            std::ofstream file("traj" + std::to_string(i+1) + ".txt", std::ios::app);
            file << t << "\t" << positions[i].x << "\t" << positions[i].y << "\t" << positions[i].z << std::endl;
        }
    }
};


char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

int main(int argc, char* argv[]) {
    if (!cmdOptionExists(argv, argv+argc, "--th") || !cmdOptionExists(argv, argv+argc, "--file"))
    {
        std::cout<<"require --file and --th args\n--dump for dumping res in file (optional)"<<std::endl;
        return 1;
    }
    bool dump = cmdOptionExists(argv, argv+argc, "--dump");
    int NUM_THREADS = std::stoi(getCmdOption(argv, argv+argc, "--th"));
    char* filename = getCmdOption(argv, argv+argc, "--file");
    const double tau = 0.1;
    const double total_time = 10.0;


    NBodySystem system(filename, NUM_THREADS);
    if(dump)
        system.clean_up();

    int steps = total_time / tau;
    std::vector<double> time_diffs;

    for (int step = 0; step <= steps; ++step) {
        double t = step * tau;
        if ((step % 10 == 0) && dump) {
            system.saveTrajectories(t);
        }
        double start = omp_get_wtime();
        system.computeStep(tau);
        double end = omp_get_wtime();
        time_diffs.push_back(end-start);
    }

    double sum = std::accumulate(time_diffs.begin(), time_diffs.end(), 0.0);

    std::cout <<sum / time_diffs.size()<<std::endl;
    return 0;
    // g++ -fopenmp -o program example.cpp
}