#ifndef ENERGY_CHANGES_H
#define ENERGY_CHANGES_H

#include <cmath>
#include "global.h"
using namespace Eigen;

bool accept_move(int thread_num, const double& delta_E) {
    return generators[thread_num].disReal() <= exp(-delta_E);
}

double delta_E_other(std::vector<Vector3i>& polym, std::vector<Vector3i>& other_polym, int site, Vector3i prop_move1, int thread_num, int lin_pol) {
    double energy_change = 0;
    int stage = thread_num % number_of_stages;
 //   int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)

    //Energies from contacts if (site+1)%pol_length is a physical monomer
    if((site+1)%reduction_factor==0){
        int red_i_plus_one=((site + 1) % pol_length)/reduction_factor;
        if (locations[thread_num].find({ prop_move1[0],prop_move1[1],prop_move1[2] }) != locations[thread_num].end()) {
            for (auto elem : locations[thread_num][{ prop_move1[0],prop_move1[1],prop_move1[2] }]) {
                energy_change += Interaction_E[elem][red_i_plus_one];
            }
        }
        if (locations_lin[thread_num].find({ prop_move1[0],prop_move1[1],prop_move1[2] }) != locations[thread_num].end()) {
            for (auto elem : locations_lin[thread_num][{ prop_move1[0],prop_move1[1],prop_move1[2] }]) {
                energy_change += Interaction_E[elem][red_i_plus_one];
            }
        }
        for (auto elem : locations[thread_num][{ polym[(site + 1) % pol_length][0],polym[(site + 1) % pol_length][1],polym[(site + 1) % pol_length][2] }]) {
            if (elem != red_i_plus_one) {
                energy_change -= Interaction_E[elem][red_i_plus_one];
            }
        }
        for (auto elem : locations_lin[thread_num][{ polym[(site + 1) % pol_length][0],polym[(site + 1) % pol_length][1],polym[(site + 1) % pol_length][2] }]) {
            if (elem != red_i_plus_one) {
                energy_change -= Interaction_E[elem][red_i_plus_one];
            }
        }
    }


    //Non-contact energy contributions
    if (lin_length[stage] !=0) {
        if (site + 1 == oriC && lin_pol) { // lin_polymer oriC
//            if (pol_close_pole(site+1,thread_num)) {
//                energy_change += beta[thread_num%number_of_stages] * (std::abs(prop_move1[2]-pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
//            else {
//                energy_change += alpha[thread_num%number_of_stages] * (std::abs(prop_move1[2]-pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
            energy_change += beta[thread_num%number_of_stages] * (prop_move1[2]  - polym[site + 1][2] );
            energy_change += beta2[thread_num%number_of_stages] *(pow(prop_move1[2] -xp_z_far_simunits[thread_num%number_of_stages],2) - pow(polym[site + 1][2] -xp_z_far_simunits[thread_num%number_of_stages],2));
        }
        else if (site + 1 == oriC && !lin_pol) {  // polymer oriC
//            if (!pol_close_pole(site+1,thread_num)) {
//                energy_change += beta[thread_num%number_of_stages] * (std::abs(prop_move1[2]-pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
//            else {
//                energy_change += alpha[thread_num%number_of_stages] * (std::abs(prop_move1[2]-pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
            energy_change += alpha[thread_num%number_of_stages] * (prop_move1[2] - polym[site + 1][2] );
            energy_change += alpha2[thread_num%number_of_stages] * (pow(prop_move1[2] -xp_z_close_simunits[thread_num%number_of_stages],2) - pow(polym[site + 1][2] -xp_z_close_simunits[thread_num%number_of_stages],2));
        }
    }
    else { //single chromosome
        if (site + 1 == oriC) {
//            if (pol_close_pole(site+1,thread_num)) {
//                energy_change += alpha[thread_num % number_of_stages] * (std::abs(prop_move1[2] - pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
            energy_change += alpha[thread_num%number_of_stages] * (prop_move1[2] - polym[site + 1][2]);
            energy_change += alpha2[thread_num%number_of_stages] * (pow(prop_move1[2] -xp_z_close_simunits[thread_num%number_of_stages],2) - pow(polym[site + 1][2] -xp_z_close_simunits[thread_num%number_of_stages],2));
        }
    }

    int moved_site1 = (site + 1) % pol_length;
    stage = thread_num % number_of_stages;

    //energy differences from means
    if(energy_mean_map[stage].find(moved_site1) != energy_mean_map[stage].end()){ //if site moved_site1 is constrained
        if(not is_replicated[thread_num][moved_site1]) {
            energy_change += energy_mean_map[stage][moved_site1] * (prop_move1[2] - polym[moved_site1][2]);
        }
        else if(include_replicated_in_mean){
            energy_change += energy_mean_map[stage][moved_site1] * (prop_move1[2] - polym[moved_site1][2]) / 2.;
        }
    }
    //energy differences from separations
    if(is_replicated[thread_num][moved_site1] && energy_separation_map[stage].find(moved_site1) != energy_separation_map[stage].end()){
        // keys of energy_separation_map should only contain replicated sites (two needed to calculate separation)
        // should add a checker for this here DONE
        energy_change += energy_separation_map[stage][moved_site1] * (std::abs(other_polym[moved_site1][2] - prop_move1[2]) - std::abs(other_polym[moved_site1][2] - polym[moved_site1][2]));
    }

    return energy_change;
}

double delta_E_crankshaft(std::vector<Vector3i>& polym, std::vector<Vector3i>& other_polym, int site, Vector3i prop_move1, Vector3i prop_move2, int thread_num, int lin_pol) {
    double energy_change1 = delta_E_other(polym, other_polym, site, prop_move1, thread_num, lin_pol);
    double energy_change2 = delta_E_other(polym, other_polym, site+1, prop_move2, thread_num, lin_pol);
    return energy_change1 + energy_change2;
}

#endif
