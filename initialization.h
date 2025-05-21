#ifndef INIT_H
#define INIT_H

#include "global.h"
#include "moves.h"
using namespace Eigen;

/*
All of these maps store only physical monomers, ie polymer index/reduction_factor
Makes maps 4*smaller, meaning searches are a lot faster.
Note that polymer and lin_polymer arrays use full indeces; these keep track of unphysical monomers
*/
std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts(number_of_threads);
std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations(number_of_threads);

std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts_inter(number_of_threads);
std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts_lin(number_of_threads);
std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations_lin(number_of_threads);


void initialize(int thread_num, int step) {
    std::cout << "Initializing thread " << thread_num << " with config " << step << std::endl;


//    int stage = thread_num % number_of_stages;
    int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)

    // Set the values of is_replicated//
    for (int i = 0; i < pol_length; i++) {
        if(lin_length[thread_num]==0){
            is_replicated[thread_num][i] = false;
            continue;
        }
        if (oriC + lin_length[stage] >= pol_length) {
            if (i < oriC - lin_length[stage] && i > (oriC + lin_length[stage]) % pol_length) {
                is_replicated[thread_num][i] = false;
            }
        }
        else {
            if (i > oriC + lin_length[stage]) {
                is_replicated[thread_num][i] = false;
            }
            else if (i < oriC - lin_length[stage]) {
                is_replicated[thread_num][i] = false;
            }
        }
    }

    // Construct initial configurations if these are not taken from saved configurations//
    if(!initConfig) {
        int xi{0}, yi{0}, zi{0};
        Vector3i monomer;
        for (int i = 0; i < pol_length; i++) {
            if (i % 2 == 0) {
                zi -= 1;
            } else {
                zi += 1;
            }
            monomer = {xi, yi, zi};
            polymer[thread_num].push_back(monomer);
            lin_polymer[thread_num].push_back(monomer);

            //JM: The values of the linear polymer are set to zero if they are past the replication fork
            if (oriC + lin_length[stage] >= pol_length) {
                if (i < oriC - lin_length[stage] && i > (oriC + lin_length[stage]) % pol_length) {
                    lin_polymer[thread_num][i] = {0, 0, 0};
                }
            } else {
                if (i > oriC + lin_length[stage]) {
                    lin_polymer[thread_num][i] = {0, 0, 0};
                } else if (i < oriC - lin_length[stage]) {
                    lin_polymer[thread_num][i] = {0, 0, 0};
                }
            }
        }
    }

    // Construct initial configurations from saved configurations//
    if (initConfig) {
        check_input_configuration_compatibility(configuration_data_folder);
        read_configuration("ring",  thread_num, step);
        read_configuration("strand", thread_num, step);
    }


    // set contacts and locations. ONLY TRACK PHYSICAL CONTACTS. //
    for (int i = 0; i < pol_length; i+=reduction_factor) { //Ring ring contacts
        int red_i = i/reduction_factor;
        if (locations[thread_num].find({ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }) != locations[thread_num].end()) {
            for (auto elem : locations[thread_num][{ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }]) {
                contacts[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
            }
        }
        locations[thread_num][{polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2]}].push_back(red_i);
    }
    if (lin_length[stage] != 0) {
        if (lin_length[stage] < pol_length/2) {
            if (oriC + lin_length[stage] >= pol_length) { //need periodic bcs
                for (int i = 0; i < pol_length; i+=reduction_factor) {
                    int red_i=i/reduction_factor;
                    if (i > oriC - lin_length[stage] || i < (oriC + lin_length[stage]) % pol_length) {
                        //Lin lin
                        if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}) != locations_lin[thread_num].end()) {
                            for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}]) {
                                contacts_lin[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                            }
                        }
                        //Ring lin
                        if (locations[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}) != locations[thread_num].end()) {
                            for (auto elem : locations[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}]) {
                                contacts_inter[thread_num][{elem, red_i}] = 0;
                            }
                        }
                        locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}].push_back(red_i);
                    }
                }
            }
            else { //replication not past zero; no periodicity needed
                for (int i = 0; i < pol_length; i+=reduction_factor) {
                    int red_i=i/reduction_factor;

                    if (i > oriC - lin_length[stage] && i < oriC + lin_length[stage]) {
                        //lin lin
                        if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}) !=locations_lin[thread_num].end()) {
                            for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}]) {
                                contacts_lin[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                            }
                        }
                        //ring lin
                        if (locations[thread_num].find({ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }) != locations[thread_num].end()) {
                            for (auto elem : locations[thread_num][{ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }]) {
                                contacts_inter[thread_num][{elem, red_i}] = 0;
                            }
                        }
                        locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}].push_back(red_i);
                    }
                }
            }
        }
        else { //fully replicated
            for (int i = 0; i < pol_length; i+=reduction_factor) {
                int red_i=i/reduction_factor;

                if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}) !=locations_lin[thread_num].end()) {
                    for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}]) {
                        contacts_lin[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                    }
                }
                if (locations[thread_num].find({ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }) != locations[thread_num].end()) {
                    for (auto elem : locations[thread_num][{ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }]) {
                        contacts_inter[thread_num][{elem, red_i}] = 0;
                    }
                }
                locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}].push_back(red_i);
            }
        }
    }
}


void burn_in(int thread_num, int n_steps) {
    for (int m = 0; m < n_steps ; m++) {   //burn-in
        move(thread_num, m);
    }
//    int stage = thread_num % number_of_stages;
    int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)

    std::vector<double> zeroVec(bin_num, 0);
    std::fill(total_contacts[thread_num].begin(), total_contacts[thread_num].end(), zeroVec);
    z_close[thread_num] = 0;
    z_far[thread_num] = 0;
    z_close_squared[thread_num] = 0;
    z_far_squared[thread_num] = 0;
    for(int i=0; i<sites_constrained_mean.size(); i++){
        z_mean_data[thread_num][i] = 0;
    }
    for(int i=0; i<sites_constrained_separation.size(); i++){
        z_separation_data[thread_num][i] = 0;
    }
    //z[thread_num] = 0;
    //z_lin[thread_num] = 0;
    locations[thread_num].clear();
    locations_lin[thread_num].clear();
    contacts[thread_num].clear();
    contacts_lin[thread_num].clear();
    contacts_inter[thread_num].clear();

    //Fill in contacts and locations
    for (int i = 0; i < pol_length; i+=reduction_factor) { //Ring ring contacts
        int red_i = i/reduction_factor;
        if (locations[thread_num].find({ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }) != locations[thread_num].end()) {
            for (auto elem : locations[thread_num][{ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }]) {
                contacts[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
            }
        }
        locations[thread_num][{polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2]}].push_back(red_i);
    }
    if (lin_length[stage] != 0) {
        if (lin_length[stage] < pol_length/2) {
            if (oriC + lin_length[stage] >= pol_length) { //need periodic bcs
                for (int i = 0; i < pol_length; i+=reduction_factor) {
                    int red_i=i/reduction_factor;
                    if (i > oriC - lin_length[stage] || i < (oriC + lin_length[stage]) % pol_length) {
                        //Lin lin
                        if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}) != locations_lin[thread_num].end()) {
                            for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}]) {
                                contacts_lin[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                            }
                        }
                        //Ring lin
                        if (locations[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}) != locations[thread_num].end()) {
                            for (auto elem : locations[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}]) {
                                contacts_inter[thread_num][{elem, red_i}] = 0;
                            }
                        }
                        locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}].push_back(red_i);
                    }
                }
            }
            else { //replication not past zero; no periodicity needed
                for (int i = 0; i < pol_length; i+=reduction_factor) {
                    int red_i=i/reduction_factor;
                    if (i > oriC - lin_length[stage] && i < oriC + lin_length[stage]) {
                        //lin lin
                        if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}) !=locations_lin[thread_num].end()) {
                            for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}]) {
                                contacts_lin[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                            }
                        }
                        //ring lin
                        if (locations[thread_num].find({ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }) != locations[thread_num].end()) {
                            for (auto elem : locations[thread_num][{ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }]) {
                                contacts_inter[thread_num][{elem, red_i}] = 0;
                            }
                        }
                        locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}].push_back(red_i);
                    }
                }
            }
        }
        else { //fully replicated
            for (int i = 0; i < pol_length; i+=reduction_factor) {
                int red_i=i/reduction_factor;
                if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}) !=locations_lin[thread_num].end()) {
                    for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}]) {
                        contacts_lin[thread_num][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                    }
                }
                if (locations[thread_num].find({ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }) != locations[thread_num].end()) {
                    for (auto elem : locations[thread_num][{ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }]) {
                        contacts_inter[thread_num][{elem, red_i}] = 0;
                    }
                }
                locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}].push_back(red_i);
            }
        }
    }
}

#endif
