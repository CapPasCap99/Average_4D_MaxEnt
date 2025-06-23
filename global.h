#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

template <typename T> void hash_combine(std::size_t& seed, const T& v) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& pair) const noexcept {
        std::size_t seed = 0;
        hash_combine(seed, pair.first);
        hash_combine(seed, pair.second);
        return seed;
    }
};
struct vec_hash {
    std::size_t operator()(const std::vector<int> vec) const noexcept {
        std::size_t seed = 0;
        for (auto el : vec) {
            hash_combine(seed,el);
        }
        return seed;
    }
};

extern std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts;
extern std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations;

extern std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts_inter;  //order: <polymer, lin_polymer>
extern std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts_lin;
extern std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations_lin;


extern const int total_simulations;
extern const int number_of_stages;
extern const int stage_count;
extern const int replicates_per_stage;
extern const int number_of_threads;
extern std::vector<RandomGenerator> generators;
extern const std::vector<int> stages;
extern const std::vector<int> reps;
extern const std::vector<double> cell_len;
extern const int mc_moves;
extern double learning_rate;
extern double learning_rate_far;
extern double learning_rate_close;
extern double learning_rate_far_var;
extern double learning_rate_close_var;

extern double learning_rate_means;
extern double learning_rate_separations;

extern std::vector<int> sites_constrained_mean;
extern std::vector<int> sites_constrained_separation;
extern std::unordered_map<int, int> sites_constrained_mean_map; // GG: keeps index of a site in the "sites_constrained_mean" list
extern std::unordered_map<int, int> sites_constrained_separation_map;
extern bool include_replicated_in_mean;

extern std::vector<std::vector<bool>> is_constrained_ori;
extern std::vector<std::vector<bool>> is_constrained_mean;
extern std::vector<std::vector<bool>> is_constrained_separation;

extern std::vector<std::vector<double>> target_means;
extern std::vector<std::vector<double>> target_separations;

extern std::vector<std::vector<double>> z_mean_data;
extern std::vector<std::vector<double>> z_separation_data;
extern std::vector<std::vector<double>> z_mean_data_tot;
extern std::vector<std::vector<double>> z_separation_data_tot;

extern std::string dir;
extern std::string output_folder;
extern std::string HiC_file;
extern bool initConfig;
extern std::string configuration_data_folder;

extern const std::string bacteria_name;
extern const int bin_num;
extern const int pol_length;
extern const int reduction_factor;
extern const int fine;
extern const std::pair<int, int> ends;
extern std::vector<std::vector<Eigen::Vector3i>> polymer;
extern std::vector<std::vector<Eigen::Vector3i>> lin_polymer;
extern std::vector<std::vector<bool>> is_replicated;
extern std::vector< std::vector<double>> Interaction_E;
extern std::vector<double> alpha;
extern std::vector<double> beta;
extern std::vector<double> alpha2;
extern std::vector<double> beta2;
extern std::vector<int> pole;

extern std::vector<std::vector<double>> energ_coeff_mean;
extern std::vector<std::vector<double>> energ_coeff_separation;
extern std::vector<std::unordered_map<int, double>> energy_mean_map;
extern std::vector<std::unordered_map<int, double>> energy_separation_map;

extern const int res;
extern std::vector<double> z_close;
extern std::vector<double> z_far;
extern std::vector<double> z_close_tot;
extern std::vector<double> z_far_tot;
extern std::vector<double> z_close_squared;
extern std::vector<double> z_far_squared;
extern std::vector<double> z_close_squared_tot;
extern std::vector<double> z_far_squared_tot;
extern std::vector<double> z;
extern std::vector<double> z_lin;
extern std::vector<double> z_tot;
extern std::vector<double> z_lin_tot;
extern std::vector<double> xp_z_close;
extern std::vector<double> xp_z_far;
extern std::vector<double> xp_z_close_var;
extern std::vector<double> xp_z_far_var;

extern std::vector<double> xp_z_close_simunits;
extern std::vector<double> xp_z_far_simunits;

extern std::vector<double> xp_close_max;
extern std::vector<double> xp_far_min;
extern std::vector<double> z_far_var;
extern std::vector<double> z_close_var;

extern bool boundary_cond;
extern double radius;
extern std::vector<double> length;
extern const int x_axis;
extern const int y_axis;
extern const int z_axis;
extern std::vector<double> offset;
extern std::vector<double> offset_z;
extern bool constrain_pol;
extern bool use_fork_distribution;

extern std::vector<int> lin_length;
extern const int oriC;

extern std::vector<std::vector< std::vector<double>>> total_contacts;
extern std::vector< std::vector<double>> final_contacts;
extern std::vector<std::vector<double>> xp_contacts;

#endif
