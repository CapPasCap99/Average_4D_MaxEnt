#ifndef MV_H
#define MV_H

#include <algorithm>
#include "global.h"
#include "energy_changes.h"

//VERIFY THAT THE RIGHT STAGES OF REPLICATION ARE BEING CONSIDERED FOR  THE MOVES WITH COUTS

//Note locations tracks physical monomers, ie pol index/reduction_factor. Same as indeces of total_contacts and contacts
void update_contacts(int thread_num, int moved_site, Eigen::Vector3i prop_move, int m){
    if (moved_site%reduction_factor != 0){ //unphysical; not tracked
        return;
    }
    int red_site=moved_site/reduction_factor;
    //Throw away old contacts_lin[thread_num], at the same time update contact frequency map
    for (auto elem : locations[thread_num][{polymer[thread_num][moved_site][0], polymer[thread_num][moved_site][1], polymer[thread_num][moved_site][2]}]) {
        if (elem != red_site) {
            int i = std::min(elem, red_site);
            int j = std::max(elem, red_site);
            total_contacts[thread_num][i][j] += double(m) - contacts[thread_num][{i,j}];
            // GG: here one can save other total contacts
            contacts[thread_num].erase({i, j});
        }
    }
    for (auto elem : locations_lin[thread_num][{polymer[thread_num][moved_site][0], polymer[thread_num][moved_site][1], polymer[thread_num][moved_site][2]}]) {
        if (elem != red_site) {
            int i = std::min(elem, red_site);
            int j = std::max(elem, red_site);
            total_contacts[thread_num][i][j] += double(m) - contacts_inter[thread_num][{red_site, elem}]; // [{index within circular, index within linear}]
            // GG: here one can save other total contacts
            contacts_inter[thread_num].erase({red_site, elem});
        }
    }
    //put in new contacts[thread_num]
    if (locations[thread_num].find({prop_move[0], prop_move[1], prop_move[2]}) != locations[thread_num].end()) {
        for (auto elem : locations[thread_num][{prop_move[0], prop_move[1], prop_move[2]}]) {
            contacts[thread_num][{std::min(elem, red_site),std::max(elem, red_site)}] = m;
        }
    }
    if (locations_lin[thread_num].find({prop_move[0], prop_move[1], prop_move[2]}) != locations_lin[thread_num].end()) {
        for (auto elem : locations_lin[thread_num][{prop_move[0], prop_move[1], prop_move[2]}]) {
            contacts_inter[thread_num][{red_site, elem}] = m; // [{index within circular, index within linear}]
        }
    }
}


void update_contacts_lin(int thread_num, int moved_site, Eigen::Vector3i prop_move, int m){
    if (moved_site%reduction_factor != 0){ //unphysical; not tracked
        return;
    }
    int red_site=moved_site/reduction_factor;
    //Throw away old contacts_lin[thread_num], at the same time update contact frequency map
    for (auto elem : locations_lin[thread_num].find({ lin_polymer[thread_num][moved_site][0],lin_polymer[thread_num][moved_site][1],lin_polymer[thread_num][moved_site][2] })->second) {
        if (elem != red_site) {
            int i = std::min(elem, red_site);
            int j = std::max(elem, red_site);
            total_contacts[thread_num][i][j] += double(m) - contacts_lin[thread_num][{i, j}];
            // GG: here one can save other total contacts
            contacts_lin[thread_num].erase({i, j});
        }
    }
    for (auto elem : locations[thread_num][{ lin_polymer[thread_num][moved_site][0],lin_polymer[thread_num][moved_site][1],lin_polymer[thread_num][moved_site][2] }]) {
        if (elem != red_site) {
            int i = std::min(elem, red_site);
            int j = std::max(elem, red_site);
            total_contacts[thread_num][i][j] += double(m) - contacts_inter[thread_num][{elem, red_site}];
            // GG: here one can save other total contacts
            contacts_inter[thread_num].erase({ elem, red_site });
        }
    }
    //put in new contacts[thread_num]
    if (locations_lin[thread_num].find({ prop_move[0],prop_move[1],prop_move[2] }) != locations_lin[thread_num].end()) {
        for (auto elem : locations_lin[thread_num].find({ prop_move[0],prop_move[1],prop_move[2] })->second) {
            contacts_lin[thread_num][{std::min(elem, red_site), std::max(elem, red_site)}] = m;
        }
    }
    if (locations[thread_num].find({ prop_move[0],prop_move[1],prop_move[2] }) != locations[thread_num].end()) {
        for (auto elem : locations[thread_num][{ prop_move[0],prop_move[1],prop_move[2] }]) {
            contacts_inter[thread_num][{elem, red_site}] = m;
        }
    }
}


void update_locations(int thread_num, int moved_site, Eigen::Vector3i prop_move){
    if (moved_site%reduction_factor != 0){ //unphysical; not tracked
        return;
    }
    int red_site=moved_site/reduction_factor;

    //update hash map locations[thread_num]
    std::vector<int> site_coords = { polymer[thread_num][moved_site][0], polymer[thread_num][moved_site][1], polymer[thread_num][moved_site][2] };
    std::vector<int> prop_move_vec = {prop_move[0], prop_move[1], prop_move[2]};
    // remove old locations
    if (locations[thread_num][site_coords].size() == 1) {
        locations[thread_num].erase(site_coords);
    }
    else {
        locations[thread_num][site_coords].erase(std::find(locations[thread_num][site_coords].begin(), locations[thread_num][site_coords].end(), red_site));
    }
    // save new locations
    if (locations[thread_num].find(prop_move_vec) != locations[thread_num].end()) {
        locations[thread_num][prop_move_vec].push_back(red_site);
    }
    else {
        locations[thread_num][prop_move_vec] = { red_site };
    }
}


void update_locations_lin(int thread_num, int moved_site, Eigen::Vector3i prop_move){
    if (moved_site%reduction_factor != 0){ //unphysical; not tracked
        return;
    }
    int red_site=moved_site/reduction_factor;
    //update hash map locations[thread_num]
    std::vector<int> site_coords = { lin_polymer[thread_num][moved_site][0],lin_polymer[thread_num][moved_site][1],lin_polymer[thread_num][moved_site][2] };
    std::vector<int> prop_move_vec = {prop_move[0], prop_move[1], prop_move[2]};
    // remove old locations
    if (locations_lin[thread_num][site_coords].size() == 1) {
        locations_lin[thread_num].erase(site_coords);
    }
    else {
        locations_lin[thread_num][site_coords].erase(std::find(locations_lin[thread_num][site_coords].begin(), locations_lin[thread_num][site_coords].end(), red_site));
    }
    // save new locations
    if (locations_lin[thread_num].find(prop_move_vec) != locations_lin[thread_num].end()) {
        locations_lin[thread_num][prop_move_vec].push_back(red_site);
    }
    else {
        locations_lin[thread_num][prop_move_vec] = { red_site };
    }
}


bool constr_move(int thread_num, Eigen::Vector3i prop_move, int lin_pol, int site) {
    int stage = thread_num % number_of_stages;
//    int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)
    if (constrain_pol && site == oriC) {
        if (lin_pol) {
            if (std::abs(prop_move[2] - offset_z[stage]) > std::abs(polymer[thread_num][site][2] - offset_z[stage])) {
                return false; //If the proposed move in the linear polymer brings it further away from the cell center than the circular polymer ori, deny move
            }
            else {
                return true;
            }
        }
        else {
            if(prop_move[2] > offset_z[stage]){ //GG: to keep the "focus closest to a pole" in the left half of the cell. (that's how the experimental data for ecoli is processed)
                return false; //GG: ok, I see now that this condition was already included below
            }
            if (std::abs(prop_move[2] - offset_z[stage]) < std::abs(lin_polymer[thread_num][site][2] - offset_z[stage]) || prop_move[2] > 0) {
                return false; //If the proposed move on the circular polymer brings it closer to the cell center than the linear polymer ori, or it moves past the midpoint, deny move
            }
            else {
                return true;
            }
        }
    }
    else {return true;} //JM: if false, the proposed move will be rejected
}


bool check_boundary_rest(Eigen::Vector3i prop_move1, int thread_num) {
    int stage = thread_num % number_of_stages;
//    int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)
    if (boundary_cond == 1) {
        bool accept_1;
        // GG: check for prop_move1
        if (std::abs(prop_move1[2] - offset_z[stage]) <= length[stage] / 2) { //GG: is inside the cylinder, not in a cap (z-coord)
            accept_1 = (pow(prop_move1[0]+offset[0], 2) + pow(prop_move1[1]+offset[1], 2) <= pow(radius, 2));
        }
        else if (std::abs(prop_move1[2] - offset_z[stage]) <= length[stage] / 2 + radius) { // GG: is in one of the caps (z-coord)
            accept_1 = (pow(prop_move1[0]+offset[0], 2) + pow(prop_move1[1]+offset[1], 2) + pow(std::abs(prop_move1[2] - offset_z[stage]) - length[stage] / 2, 2) <= pow(radius, 2));
        }
        else {
            accept_1 = 0;
        }

        return accept_1;
    }
    else { return 1; }
}

void junction_move(int thread_num, int site, int &m, int lin_pol) {
    int lin_site;
    int stage = thread_num % number_of_stages;
//    int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)

    //JM: Note: scheme below only works if the origin position is past the midway point of the polymer
    if (oriC + lin_length[stage] >= pol_length) { //JM: Replication has passed the 0-coordinate
        if (site == oriC - lin_length[stage] - 1) {//JM: Site is behind Ori
            lin_site = site + 2; //JM:I think this is one site past the junction point of the linear polymer
        }
        else if (site == (oriC + lin_length[stage] - 1) % pol_length) { //JM Site is ahead of ori. This also assigns a linear site if replication is 100%. How is this later used?
            lin_site = site; //JM: I think this is one site past the junction point of the linear polymer
        }
    }
    else {
        lin_site = site + int((oriC + lin_length[stage] - 1 - site)/lin_length[stage]); //JM: equals site+1 if it site is before ori, site-1 if it is after ori
    }
    //JM: If-statement below: if
    // - Polymer forms a two-monomer loop at site and site+2
    // - The linear polymer doesn't also form part of the loop
    // - The linear polymer doesn't make an 180 degree angle with the loop
    if (polymer[thread_num][site] == polymer[thread_num][(site + 2)%pol_length] && lin_polymer[thread_num][lin_site] != polymer[thread_num][site] && lin_polymer[thread_num][lin_site] != 2 * polymer[thread_num][(site + 1)%pol_length] - polymer[thread_num][site]) {
        Eigen::Vector3i prop_move1;
        prop_move1 = polymer[thread_num][site] + lin_polymer[thread_num][lin_site] - polymer[thread_num][(site + 1)%pol_length];
        if (accept_move(thread_num, delta_E_other(polymer[thread_num], lin_polymer[thread_num], site, prop_move1, thread_num, lin_pol)) && check_boundary_rest(prop_move1, thread_num)) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts(thread_num, (site + 1)%pol_length, prop_move1, m);
            //update hash map locations[thread_num]
            update_locations(thread_num, (site + 1)%pol_length, prop_move1);

            polymer[thread_num][(site + 1)%pol_length] = prop_move1;
            lin_polymer[thread_num][(site + 1)%pol_length] = prop_move1;
        }
        m++; //monte carlo step counted (whether move accepted or not)
    }
        //JM: If-statement below: if
        //- Linear polymer and polymer form loop
        //- Other polymer site is at a 90 degree angle
    else if ((lin_polymer[thread_num][lin_site] == polymer[thread_num][(site + 2)%pol_length] || lin_polymer[thread_num][lin_site] == polymer[thread_num][site]) && polymer[thread_num][(site + 2)%pol_length] != polymer[thread_num][site] && polymer[thread_num][(site + 2)%pol_length] != 2 * polymer[thread_num][(site + 1)%pol_length] - polymer[thread_num][site]) {
        Eigen::Vector3i prop_move1;
        prop_move1 = polymer[thread_num][site] + polymer[thread_num][(site + 2)%pol_length] - polymer[thread_num][(site + 1)%pol_length];
        if (accept_move(thread_num, delta_E_other(polymer[thread_num], lin_polymer[thread_num], site, prop_move1, thread_num, lin_pol)) && check_boundary_rest(prop_move1, thread_num)) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts(thread_num, (site + 1)%pol_length, prop_move1, m);
            //update hash map locations[thread_num]
            update_locations(thread_num, (site + 1)%pol_length, prop_move1);

            polymer[thread_num][(site + 1)%pol_length] = prop_move1;
            lin_polymer[thread_num][(site + 1)%pol_length] = prop_move1;
        }
        m++;
    }
}

void kink_move(int thread_num, int site, int &m, int lin_pol) {
    if ((polymer[thread_num][(site + 2) % pol_length] != polymer[thread_num][site]) && (polymer[thread_num][(site + 2) % pol_length] != 2 * polymer[thread_num][(site + 1) % pol_length] - polymer[thread_num][site])) {
        Eigen::Vector3i prop_move1;
        prop_move1 = polymer[thread_num][site] + polymer[thread_num][(site + 2) % pol_length] - polymer[thread_num][(site + 1) % pol_length];
        if (accept_move(thread_num, delta_E_other(polymer[thread_num], lin_polymer[thread_num], site, prop_move1, thread_num, lin_pol)) == 1 && check_boundary_rest(prop_move1, thread_num) == 1 && constr_move(thread_num,prop_move1,lin_pol,site+1)) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts(thread_num, (site + 1)%pol_length, prop_move1, m);
            //update hash map locations[thread_num]
            update_locations(thread_num, (site + 1)%pol_length, prop_move1);

            //update polymer[thread_num]
            polymer[thread_num][(site + 1) % pol_length] = prop_move1;
        }
        m++;
    }
}

void kink_move_lin(int thread_num, int site, int &m, int lin_pol) {
    if ((lin_polymer[thread_num][(site + 2) % pol_length] != lin_polymer[thread_num][site%pol_length]) && (lin_polymer[thread_num][(site + 2) % pol_length] != 2 * lin_polymer[thread_num][(site + 1) % pol_length] - lin_polymer[thread_num][site%pol_length])) {
        Eigen::Vector3i prop_move1;
        prop_move1 = lin_polymer[thread_num][site%pol_length] + lin_polymer[thread_num][(site + 2) % pol_length] - lin_polymer[thread_num][(site + 1) % pol_length];
        if (accept_move(thread_num, delta_E_other(lin_polymer[thread_num], polymer[thread_num], site, prop_move1, thread_num, lin_pol)) == 1 && check_boundary_rest(prop_move1, thread_num) == 1 && constr_move(thread_num,prop_move1,lin_pol,site+1)) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts_lin(thread_num, (site + 1)%pol_length, prop_move1, m);
            //update hash map locations[thread_num]
            update_locations_lin(thread_num, (site + 1)%pol_length, prop_move1);

            //update lin_polymer[thread_num]
            lin_polymer[thread_num][(site + 1) % pol_length] = prop_move1;
        }
        m++;
    }
}

void crankshaft_move(int thread_num, int site, int &m, int lin_pol) {
    if ((polymer[thread_num][(site + 2) % pol_length] != polymer[thread_num][site]) && (polymer[thread_num][(site + 2) % pol_length] != 2 * polymer[thread_num][(site + 1) % pol_length] - polymer[thread_num][site]) && (polymer[thread_num][(site + 3) % pol_length] - polymer[thread_num][(site + 2) % pol_length] == polymer[thread_num][site] - polymer[thread_num][(site + 1) % pol_length])) {
        int direction = generators[thread_num].unidir();
        Eigen::Vector3i prop_move1; Eigen::Vector3i prop_move2;
        if (direction == 1) { //180 degree flip
            prop_move1 = 2 * polymer[thread_num][site] - polymer[thread_num][(site + 1) % pol_length];
            prop_move2 = 2 * polymer[thread_num][(site + 3) % pol_length] - polymer[thread_num][(site + 2) % pol_length];
        }
        else { //90 degree flip
            prop_move1 = polymer[thread_num][site] + (direction - 1) * (polymer[thread_num][(site + 1) % pol_length] - polymer[thread_num][site]).cross(polymer[thread_num][(site + 3) % pol_length] - polymer[thread_num][site]);
            prop_move2 = polymer[thread_num][(site + 3) % pol_length] + (direction - 1) * (polymer[thread_num][(site + 2) % pol_length] - polymer[thread_num][(site + 3) % pol_length]).cross(polymer[thread_num][(site + 3) % pol_length] - polymer[thread_num][site]);;
        }
        if (accept_move(thread_num, delta_E_crankshaft(polymer[thread_num], lin_polymer[thread_num], site, prop_move1, prop_move2, thread_num, lin_pol)) == 1 && check_boundary_rest(prop_move1, thread_num) == 1 && check_boundary_rest(prop_move2, thread_num) == 1 && constr_move(thread_num,prop_move1,lin_pol,site+1) && constr_move(thread_num,prop_move2,lin_pol,site+2)) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts(thread_num, (site + 1)%pol_length, prop_move1, m);
            update_contacts(thread_num, (site + 2)%pol_length, prop_move2, m);
            //update hash map locations[thread_num]
            update_locations(thread_num, (site + 1)%pol_length, prop_move1);
            update_locations(thread_num, (site + 2)%pol_length, prop_move2);

            //update polymer[thread_num]
            polymer[thread_num][(site + 1) % pol_length] = prop_move1;
            polymer[thread_num][(site + 2) % pol_length] = prop_move2;
        }
        m++;
    }
}

void crankshaft_move_lin(int thread_num, int site, int &m, int lin_pol) {
    if ((lin_polymer[thread_num][(site + 2) % pol_length] != lin_polymer[thread_num][site%pol_length]) && (lin_polymer[thread_num][(site + 2) % pol_length] != 2 * lin_polymer[thread_num][(site + 1) % pol_length] - lin_polymer[thread_num][site%pol_length]) && (lin_polymer[thread_num][(site + 3) % pol_length] - lin_polymer[thread_num][(site + 2) % pol_length] == lin_polymer[thread_num][site%pol_length] - lin_polymer[thread_num][(site + 1) % pol_length])) {
        int direction = generators[thread_num].unidir();
        Eigen::Vector3i prop_move1; Eigen::Vector3i prop_move2;
        if (direction == 1) { //180 degree flip
            prop_move1 = 2 * lin_polymer[thread_num][site%pol_length] - lin_polymer[thread_num][(site + 1) % pol_length];
            prop_move2 = 2 * lin_polymer[thread_num][(site + 3) % pol_length] - lin_polymer[thread_num][(site + 2) % pol_length];
        }
        else { //90 degree flip
            prop_move1 = lin_polymer[thread_num][site%pol_length] + (direction - 1) * (lin_polymer[thread_num][(site + 1) % pol_length] - lin_polymer[thread_num][site%pol_length]).cross(lin_polymer[thread_num][(site + 3) % pol_length] - lin_polymer[thread_num][site%pol_length]);
            prop_move2 = lin_polymer[thread_num][(site + 3) % pol_length] + (direction - 1) * (lin_polymer[thread_num][(site + 2) % pol_length] - lin_polymer[thread_num][(site + 3) % pol_length]).cross(lin_polymer[thread_num][(site + 3) % pol_length] - lin_polymer[thread_num][site%pol_length]);;
        }
        if (accept_move(thread_num, delta_E_crankshaft(lin_polymer[thread_num], polymer[thread_num], site, prop_move1, prop_move2, thread_num, lin_pol)) == 1 && check_boundary_rest(prop_move1, thread_num) == 1 && check_boundary_rest(prop_move2, thread_num) == 1  && constr_move(thread_num,prop_move1,lin_pol,site+1) && constr_move(thread_num,prop_move2,lin_pol,site+2)) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts_lin(thread_num, (site + 1)%pol_length, prop_move1, m);
            update_contacts_lin(thread_num, (site + 2)%pol_length, prop_move2, m);
            //update hash map locations[thread_num]
            update_locations_lin(thread_num, (site + 1)%pol_length, prop_move1);
            update_locations_lin(thread_num, (site + 2)%pol_length, prop_move2);

            //update lin_polymer[thread_num]
            lin_polymer[thread_num][(site + 1) % pol_length] = prop_move1;
            lin_polymer[thread_num][(site + 2) % pol_length] = prop_move2;
        }
        m++;
    }
}

void loop_move(int thread_num, int site, int &m, int lin_pol) {
    if (polymer[thread_num][site] == polymer[thread_num][(site + 2) % pol_length]) {
        int direction = generators[thread_num].unidir_loop();
        Eigen::Vector3i rotated_vector(3); Eigen::Vector3i prop_move1;
        for (int i = 0; i < 3; i++) { // rotate loop in one of 5 possible new directions
            rotated_vector[(i + direction / 2) % 3] = (-2 * (direction % 2) + 1) * (polymer[thread_num][(site + 1) % pol_length] - polymer[thread_num][site])[i];
        }
        prop_move1 = polymer[thread_num][site] + rotated_vector;
        if (accept_move(thread_num, delta_E_other(polymer[thread_num], lin_polymer[thread_num], site, prop_move1, thread_num, lin_pol)) == 1 && check_boundary_rest(prop_move1, thread_num) == 1 && constr_move(thread_num,prop_move1,lin_pol,site+1)) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts(thread_num, (site + 1)%pol_length, prop_move1, m);
            //update hash map locations[thread_num]
            update_locations(thread_num, (site + 1)%pol_length, prop_move1);

            //update polymer[thread_num]
            polymer[thread_num][(site + 1) % pol_length] = prop_move1;
        }
        m++;
    }
}

void loop_move_lin(int thread_num, int site, int &m, int lin_pol) {
    if (lin_polymer[thread_num][site% pol_length] == lin_polymer[thread_num][(site + 2) % pol_length]) {
        int direction = generators[thread_num].unidir_loop();
        Eigen::Vector3i rotated_vector(3); Eigen::Vector3i prop_move1;
        for (int i = 0; i < 3; i++) { // rotate loop in one of 5 possible new directions
            rotated_vector[(i + direction / 2) % 3] = (-2 * (direction % 2) + 1) * (lin_polymer[thread_num][(site + 1) % pol_length] - lin_polymer[thread_num][site% pol_length])[i];
        }
        prop_move1 = lin_polymer[thread_num][site% pol_length] + rotated_vector;
        if (accept_move(thread_num, delta_E_other(lin_polymer[thread_num], polymer[thread_num], site, prop_move1, thread_num, lin_pol)) == 1 && check_boundary_rest(prop_move1, thread_num) == 1 && constr_move(thread_num,prop_move1,lin_pol,site+1)) {

            //Throw away old contacts[thread_num], at the same time update contact frequency map
            update_contacts_lin(thread_num, (site + 1)%pol_length, prop_move1, m);
            //update hash map locations[thread_num]
            update_locations_lin(thread_num, (site + 1)%pol_length, prop_move1);

            //update lin_polymer[thread_num]
            lin_polymer[thread_num][(site + 1) % pol_length] = prop_move1;
        }
        m++;
    }
}

void move(int thread_num, int &m) {
    int lin_pol = generators[thread_num].weightedpol(); //JH: prob proportional to replicated length
    //Pick site on linear or full chromosome:
    int site = (lin_pol==1) ? generators[thread_num].unisitelin() : generators[thread_num].unisitering();

    int stage = thread_num % number_of_stages;
//    int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)

    int old_m = m; //to make sure that for m%100==0 the z-position data is not saved again when no step was made
    //JM it seems like the relative position of the site is calculated twice: once to decide the move, and once while performing the move. This can be done more efficiently...
    if (lin_length[stage] != 0) {
        if (lin_pol) { //picks linear polymer
            if (lin_length[stage] < pol_length/2) { //replicating chromosome
                if (oriC + lin_length[stage] >= pol_length + 3) { //JM: replication has passed zero coordinate, plus three?!
                    if (site >= oriC - lin_length[stage] ||
                        site <= (oriC + lin_length[stage] - 3) % pol_length) { //JM: Site is on the replicated portion of the linear polymer, and not at the junction site
                        int action = generators[thread_num].unimove();
                        if (action == 0) {
                            kink_move_lin(thread_num, site, m, lin_pol);
                        } else if (action == 1) {
                            crankshaft_move_lin(thread_num, site, m, lin_pol);
                        } else {
                            loop_move_lin(thread_num, site, m, lin_pol);
                        }
                    } else if (site == (oriC + lin_length[stage] - 2) % pol_length) {
                        int action2 = generators[thread_num].unimove2();
                        if (action2 == 0) {
                            kink_move_lin(thread_num, site, m, lin_pol);
                        } else {
                            loop_move_lin(thread_num, site, m, lin_pol);
                        }
                    }
                } else { //JM: replication has not passed the zero coordinate, plus three
                    if (site >= oriC - lin_length[stage] && site <= oriC + lin_length[stage] - 3) {//JM: site on the replicated portion of the linear polymer, and not at the junction site
                        int action = generators[thread_num].unimove();
                        if (action == 0) {
                            kink_move_lin(thread_num, site, m, lin_pol);
                        } else if (action == 1) {
                            crankshaft_move_lin(thread_num, site, m, lin_pol);
                        } else {
                            loop_move_lin(thread_num, site, m, lin_pol);
                        }
                    } else if (site == (oriC + lin_length[stage] - 2) % pol_length) {//JM: if the site is one closer to the replication point: no crankshaft allowed
                        int action2 = generators[thread_num].unimove2();
                        if (action2 == 0) {
                            kink_move_lin(thread_num, site, m, lin_pol);
                        } else {
                            loop_move_lin(thread_num, site, m, lin_pol);
                        }
                    }
                }
            }
            else {  //fully replicated chromosome
                int action = generators[thread_num].unimove();
                if (action == 0) {
                    kink_move_lin(thread_num, site, m, lin_pol);
                } else if (action == 1) {
                    crankshaft_move_lin(thread_num, site, m, lin_pol);
                } else {
                    loop_move_lin(thread_num, site, m, lin_pol);
                }
            }
        } else {    //picks ring polymer
            if (lin_length[stage] < pol_length/2) { //replicating chromosome
                if (oriC + lin_length[stage] >= pol_length + 3) {
                    if (site == oriC - lin_length[stage] - 1 ||
                        site == (oriC + lin_length[stage] - 1) % pol_length) {
                        junction_move(thread_num, site, m, lin_pol);
                    } else if ((site <= oriC - lin_length[stage] - 3 &&
                                site >= (oriC + lin_length[stage]) % pol_length) ||
                               site >= oriC - lin_length[stage] ||
                               site <= (oriC + lin_length[stage] - 3) % pol_length) {
                        int action = generators[thread_num].unimove();
                        if (action == 0) {
                            kink_move(thread_num, site, m, lin_pol);
                        } else if (action == 1) {
                            crankshaft_move(thread_num, site, m, lin_pol);
                        } else {
                            loop_move(thread_num, site, m, lin_pol);
                        }
                    } else if (site == oriC - lin_length[stage] - 2 ||
                               site == (oriC + lin_length[stage] - 2) % pol_length) {
                        int action2 = generators[thread_num].unimove2();
                        if (action2 == 0) {
                            kink_move(thread_num, site, m, lin_pol);
                        } else {
                            loop_move(thread_num, site, m, lin_pol);
                        }
                    }
                } else {
                    if (site == oriC - lin_length[stage] - 1 || site == oriC + lin_length[stage] - 1) {
                        junction_move(thread_num, site, m, lin_pol);
                    } else if (site <= oriC - lin_length[stage] - 3 ||
                               (site >= oriC - lin_length[stage] && site <= oriC + lin_length[stage] - 3) ||
                               site >= oriC + lin_length[stage]) {
                        int action = generators[thread_num].unimove();
                        if (action == 0) {
                            kink_move(thread_num, site, m, lin_pol);
                        } else if (action == 1) {
                            crankshaft_move(thread_num, site, m, lin_pol);
                        } else {
                            loop_move(thread_num, site, m, lin_pol);
                        }
                    } else if (site == oriC - lin_length[stage] - 2 ||
                               site == (oriC + lin_length[stage] - 2) % pol_length) {
                        int action2 = generators[thread_num].unimove2();
                        if (action2 == 0) {
                            kink_move(thread_num, site, m, lin_pol);
                        } else  {
                            loop_move(thread_num, site, m, lin_pol);
                        }
                    }
                }
            }
            else {  //fully replicated chromosome
                int action = generators[thread_num].unimove();
                if (action == 0) {
                    kink_move(thread_num, site, m, lin_pol);
                } else if (action == 1) {
                    crankshaft_move(thread_num, site, m, lin_pol);
                } else {
                    loop_move(thread_num, site, m, lin_pol);
                }
            }
        }
        if (m%res==0 && m>old_m) { ///WHAT IS RES ??????
            z_close[thread_num] += polymer[thread_num][oriC][2];
            z_far[thread_num] += lin_polymer[thread_num][oriC][2];
            z_close_squared[thread_num] += pow(polymer[thread_num][oriC][2], 2);
            z_far_squared[thread_num] += pow(lin_polymer[thread_num][oriC][2], 2);

            for(int i=0; i<sites_constrained_mean.size(); i++){
                if(not is_replicated[thread_num][sites_constrained_mean[i]]){
                    z_mean_data[thread_num][i] += polymer[thread_num][ sites_constrained_mean[i] ][2];
                }
                else if(include_replicated_in_mean){
                    z_mean_data[thread_num][i] += (polymer[thread_num][ sites_constrained_mean[i] ][2] + lin_polymer[thread_num][ sites_constrained_mean[i] ][2])/2.;
                }
            }
            for(int i=0; i<sites_constrained_separation.size(); i++) {
                if (is_replicated[thread_num][sites_constrained_separation[i]]) {
                    z_separation_data[thread_num][i] += std::abs(polymer[thread_num][sites_constrained_separation[i]][2] - lin_polymer[thread_num][sites_constrained_separation[i]][2]);
                }
            }
        }
    }
    else if(lin_pol==0) {  //simulate only one polymer // GG: with a non_replicated polymer, it still does something if lin_pol=1... (corrected now)
        int action = generators[thread_num].unimove();
        if (action == 0) {
            kink_move(thread_num, site, m, lin_pol);
        } else if (action == 1) {
            crankshaft_move(thread_num, site, m, lin_pol);
        } else {
            loop_move(thread_num, site, m, lin_pol);
        }
        if (m%res==0 && m>old_m) {
            z_close[thread_num] += polymer[thread_num][oriC][2];
            z_close_squared[thread_num] += pow(polymer[thread_num][oriC][2], 2);

            for(int i=0; i<sites_constrained_mean.size(); i++){
                z_mean_data[thread_num][i] += polymer[thread_num][ sites_constrained_mean[i] ][2];
            }
        }
    }
}

#endif
