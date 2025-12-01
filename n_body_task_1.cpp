#include <iostream>
#include <vector>
#include <fstream>
#include <omp.h>
#include <cmath>

struct Vec3 { 
    double x, y, z; 
    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
};

class NBodySystem {
    std::vector<double> masses;
    std::vector<Vec3> positions, velocities, accelerations;
    const double G = 6.67e-11;
    const double epsilon = 1e-3;
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
        clean_up();
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

int main() {
    const int NUM_THREADS = 4;  // константное количество потоков
    const double tau = 0.01;
    const double total_time = 10.0;
    
    NBodySystem system("input.txt", NUM_THREADS);
    int steps = total_time / tau;
    
    for (int step = 0; step <= steps; ++step) {
        double t = step * tau;
        if (step % 10 == 0) {
            system.saveTrajectories(t);
        }
        system.computeStep(tau);
    }
    
    return 0;
}