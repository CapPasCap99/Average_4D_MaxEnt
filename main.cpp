/*
ORIGINALLY
Code accompanying the paper 'Learning the dynamic organization of a replicating bacterial chromosome from time-course Hi-C data' by Harju et al. (2025)
Here, the energies of the MaxEnt model are obtained via iterative Monte Carlo simulations.

The code is structured as follows:

global.h
- Initiation of global variables

files.h
- Reading of input files, with checks on correct data structures. 
- Writing of output files.

Initialization.h
- Define starting set of polymer configurations, either via reading in saved configurations, or by default construction.
- Configure polymer replication stage.
- Burn in polymer configurations.

RandomGenerator.h
- Definition of the RandomGenerator class, used for random number generation during Monte Carlo simulations.

moves.h
- All polymer moves used during the Monte Carlo simulation, with supporting functions.

energy_changes.h
- Computation of energy changes associated with a Monte Carlo move.

main.cpp
- Setting of simulation parameters
- Orchestration of iterative Monte Carlo simulation

NOW MODIFIED BY C.BERAUD
to do an averaged over replication stages (0,10,30,45,60,75min) MaxEnt model.

//IN 4D CODE EACH RUN IS ONE STAGE OF REPLICATION. EACH RUN IS ACTUALLY 3D.(SHOULD MAKE SURE OF THIS) I WANT TO HAVE ALL STAGES SIMULTADE IN ONE RUN. SO I NEED TO THINK HOW TO IMPLEMENT THIS.
*/
#include <thread>
#include <chrono>
#include <iomanip>
#include "RandomGenerator.h"
#include "global.h"
#include "files.h"
#include "initialization.h"

////// Simulation properties //////
const int number_of_threads = 30; //Number run in parallel
const std::vector<int> stages = {0, 10, 30, 45, 60, 75}; // Crescentus stages: {0, 10, 30, 45, 60, 75}
const std::vector<int> reps = {0, 72, 292, 456, 624, 788};
const std::vector<double> cell_len = {1.8, 1.9, 2.2, 2.4, 2.8, 3.1};
const int number_of_stages = stages.size(); //number of replication stages
const int total_simulations = 144; // Total number of simulations
const int replicates_per_stage = total_simulations / stages.size();
const int stage_count = number_of_stages;


const std::string bacteria_name = "crescentus_separations_3"; //Determines which input files will be used
//const int bin_num = 405; //Crescentus
//const int reduction_factor = 4; //Crescentus
const int bin_num = 405; //Crescentus zoomed
const int reduction_factor = 4; //Crescentus zoomed
const int oriC = 0; // Crescentus   //position of the origin of replication. Determines where replication is initiated.
double radius{ 3.2 }; // Crescentus

const int pol_length = reduction_factor*bin_num;
// const int mc_moves = 4e6;
// const int burn_in_steps = 3e7;
// const int it_steps = 80; //JM: number of steps inverse algorithm
const int mc_moves = 4e2;
const int burn_in_steps = 3e2;
const int it_steps = 2; //JM: number of steps inverse algorithm
bool boundary_cond = true;
bool constrain_pol = true; // constrains one origin to always be the closer one to a cell pole

///////// initial configurations ////////////
bool initConfig = false;
std::string configuration_data_folder = "configurations/";
int init_config_number = 0; //iteration number of the initial configuration used
const int res{1000}; //JM: how often the mean positions are sampled

///////// initial energies //////////////
bool initEnerg = false; //load initial energies?
std::string energy_data_folder = configuration_data_folder; //only used if initErg=True
std::string energy_data_iteration = std::to_string(init_config_number);//only used if initErg=True

///////// input Hi-C data ////////////
// std::string dir = "/";
// std::string HiC_file = dir + "myfile.txt";
std::string dir = "/home/capucine/Documents/test/4D/HiC_data/averaged_normalized/";
std::string HiC_file = dir + "t_averaged_normalized_array_bio_weights.txt";

////// Learning rates //////
// Should be a dynamical parameter! //
double learning_rate{ 0.5 }; // contacts update

//Unused learning rates 
double learning_rate_close { 0 };
double learning_rate_far { 0 };
double learning_rate_close_var { 0 };
double learning_rate_far_var { 0 };
double learning_rate_means { 0.5 };
double learning_rate_separations { 0.5 };
double update_cap_factor {0.1};
int update_cap_onset {5};

////// positional constraints ///////

//const int onset_ori_constraint = 0; //Stage from which the ori constraints are imposed
std::vector<int> sites_constrained_mean = {oriC}; //sites for which constraints on mean are imposed
std::vector<int> sites_constrained_separation = {oriC}; //sites for which constraints on separation are imposed
const int n_constrained_mean = sites_constrained_mean.size();
const int n_constrained_separation = sites_constrained_separation.size();
std::unordered_map<int, int> sites_constrained_mean_map; // GG: keeps index of a site in the "sites_constrained_mean" list
std::unordered_map<int, int> sites_constrained_separation_map;
bool include_replicated_in_mean = false; //GG: include (z1+z2)/2 as contribution to mean z for replicated sites?

////// to know if we have data for a given site at a given stage (to avoid updating corresponding energies if we don't) ///////
std::vector<std::vector<bool>> is_constrained_ori(number_of_stages, std::vector<bool>(4, false));
std::vector<std::vector<bool>> is_constrained_mean(number_of_stages, std::vector<bool>(n_constrained_mean, false));
std::vector<std::vector<bool>> is_constrained_separation(number_of_stages, std::vector<bool>(n_constrained_separation, false));

std::vector<double> offset {0.5,0.5}; //JM: offsets in x and y directions, determines the center of the cylinder
std::vector<double> offset_z(total_simulations, 0);

/////// Set properties of randomly generated numbers ////////
std::vector<RandomGenerator> generators;

/////// variables to be used later ///////
std::string output_folder;

std::vector<int> pole (total_simulations,0);  // JM: stores the z value of the close cell pole, redundant in the current implementation because always the same
/////// Experimental constraints ////////
std::vector<double> xp_z_close(number_of_stages);// Experimental constraints [units of cell length]
std::vector<double> xp_z_far(number_of_stages);
std::vector<double> xp_z_close_var(number_of_stages);
std::vector<double> xp_z_far_var(number_of_stages);

std::vector<double> xp_z_close_simunits(number_of_stages); // Constraints in the units of simulation lattice (needed for calculating variance energy)
std::vector<double> xp_z_far_simunits(number_of_stages);

///////// Storing ori z-coordinate statistics  ///////////////
std::vector<double> z_close(total_simulations,0);
std::vector<double> z_far(total_simulations,0);
std::vector<double> z_close_squared(total_simulations,0);
std::vector<double> z_far_squared(total_simulations,0);
std::vector<double> z_close_var(total_simulations,0);
std::vector<double> z_far_var(total_simulations,0);

std::vector<double> z_close_tot(number_of_stages,0);
std::vector<double> z_far_tot(number_of_stages,0);
std::vector<double> z_close_squared_tot(number_of_stages,0);
std::vector<double> z_far_squared_tot(number_of_stages,0);

////////// Interaction, position, separation ENERGIES ///////////
std::vector<std::vector<double>> Interaction_E(pol_length, std::vector<double>(bin_num, 0));

std::vector<double> alpha(number_of_stages,0); // Energy close origin, for each stage of replication
std::vector<double> beta(number_of_stages,0); // Energy far origin
std::vector<double> alpha2(number_of_stages,0);// Energy close origin (variance), for each stage of replication
std::vector<double> beta2(number_of_stages,0);// Energy far origin (variance)

std::vector<std::vector<double>> energ_coeff_mean(number_of_stages, std::vector<double>(sites_constrained_mean.size(),0));// = {{}}; //[stage] equivalent of Lucas's alpha/beta, multiplies (z_close + z_far)
std::vector<std::vector<double>> energ_coeff_separation(number_of_stages, std::vector<double>(sites_constrained_separation.size(),0));// = {{}}; //[stage] separation energy - multiplies (|z_far - z_close|)
std::vector<std::unordered_map<int, double>> energy_mean_map(number_of_stages, std::unordered_map<int, double>()); //[stage] site_index -> energy coefficient
std::vector<std::unordered_map<int, double>> energy_separation_map(number_of_stages, std::unordered_map<int, double>());

/////////// Storing z-coordinate statistics of constrained sites /////////////////
std::vector<std::vector<double>> z_mean_data(total_simulations, std::vector<double>(sites_constrained_mean.size(),0));
std::vector<std::vector<double>> z_separation_data(total_simulations, std::vector<double>(sites_constrained_separation.size(),0));
std::vector<std::vector<double>> z_mean_data_tot(number_of_stages, std::vector<double>(sites_constrained_mean.size(),0));
std::vector<std::vector<double>> z_separation_data_tot(number_of_stages, std::vector<double>(sites_constrained_separation.size(),0));

//Important vectors//
std::vector<std::vector<Eigen::Vector3i>> polymer(total_simulations);
std::vector<std::vector<Eigen::Vector3i>> lin_polymer(total_simulations);

std::vector<std::vector<bool>> is_replicated(total_simulations, std::vector<bool>(pol_length, true)); // to have an easy check for which sites are replicated. Unreplicated are set to "false" in "initialization.h"

std::vector<std::vector<double>> final_contacts(pol_length, std::vector<double>(bin_num, 0));
std::vector<std::vector<std::vector<double>>> total_contacts(total_simulations, std::vector< std::vector<double>>(bin_num, std::vector<double>(bin_num, 0))); // This stores the contact frequencies during the simulation
std::vector<std::vector<double>> xp_contacts(bin_num, std::vector<double>(bin_num, 0));


//What are those?// 
std::vector<std::vector<double>> target_means;
std::vector<std::vector<double>> target_separations;
bool use_fork_distribution = true; // or true, depending on your needs

/////// cellular properties ///////
std::vector<double> length; // cellular length
std::vector<int> lin_length; // number of replicated sites (physical ones)

/////// functions to be used later ///////
void run(int thread_num, int move_num);
void normalize();
void update_E(int step);
void normalize_measured_z();
void update_alpha_beta(int step);
void update_energ_coeff_mean (int step);
void update_energ_coeff_separation (int step);
void clean_up();
void set_unconstrained_energies_to_zero();

/////// script ///////

int main() {

    std::cout << "bin number: " << bin_num << '\n' << "polymer length: " << pol_length << '\n' << "mc moves: " << mc_moves <<'\n' << "boundaries: " << boundary_cond << '\n' << "learning rate: " << learning_rate << '\n' << "saved in " << dir <<'\n';

    length.resize(total_simulations);
    lin_length.resize(total_simulations);
    for (int i = 0; i < stage_count; i++) {
        for (int r = 0; r < replicates_per_stage; r++) {
            lin_length[i * replicates_per_stage + r] = reps[i];
            length[i * replicates_per_stage + r] = cell_len[i];
        }
    }
    polymer.resize(total_simulations);
    lin_polymer.resize(total_simulations);
    is_replicated.resize(total_simulations);
    for (int i = 0; i < total_simulations; ++i) {
        is_replicated[i].resize(pol_length, true);
    }

    locations.resize(total_simulations);
    locations_lin.resize(total_simulations);
    contacts.resize(total_simulations);
    contacts_lin.resize(total_simulations);
    contacts_inter.resize(total_simulations);
    z_mean_data.resize(total_simulations, std::vector<double>(sites_constrained_mean.size(), 0));
    z_separation_data.resize(total_simulations, std::vector<double>(sites_constrained_separation.size(), 0));
    offset_z.resize(total_simulations, 0.0);


    //Read in constraints//
    //read_input_data(); //Not needed for now
    read_file(xp_contacts, HiC_file);
    std::cout << "HiC read " ;

    //Read in energies//
    if(initEnerg) {
        check_input_energies_compatibility(energy_data_folder, energy_data_iteration);
        // GG: throws an error if the stages in the input file are not the same as in the simulation. Could be modified to automatically find needed stages in the input folder...

        read_interaction_energies(energy_data_folder, energy_data_iteration);
        read_ori_energies(energy_data_folder, energy_data_iteration);
        read_position_energies(energy_data_folder, energy_data_iteration);
        read_separation_energies(energy_data_folder, energy_data_iteration);

        set_unconstrained_energies_to_zero();
    }
    std::cout << "energy read (if exist). Spoiler alert: " << initEnerg;

    //Create output directory (with date and time in the name)//
    time_t time_now = time(0);   // get time now
    struct tm * now = localtime( & time_now );
    char buffer [20];
    strftime(buffer, 20, "%Y-%m-%d_%H%M", now);
    std::string buffer_str = buffer; // converts to string
    output_folder = buffer_str + "_" + std::to_string(stages[0]);

    std::string command = "mkdir " + dir + output_folder;
    std::system(command.c_str());

    command = "mkdir " + dir + output_folder + "/Configurations";
    std::system(command.c_str());
    command = "mkdir " + dir + output_folder + "/Contacts";
    std::system(command.c_str());
    command = "mkdir " + dir + output_folder + "/Energies";
    std::system(command.c_str());
    command = "mkdir " + dir + output_folder + "/Positions";
    std::system(command.c_str());

    // create and seed random generators //
for (int i = 0; i < total_simulations; i++) {
    generators.push_back(RandomGenerator(int(time_now) + i, pol_length, lin_length[i], oriC));
}

    // set values offset_z to calculate the center of the cell while accounting for odd cell-lengths??//
for (int i = 0; i < total_simulations; i++) {
    offset_z[i] = (int(length[i] )% 2) / 2.0;
}

    //Save the input parameters of the simulation in "sim_params.txt"//
    get_sim_params();

    //Initialize configurations//
    std::cout<<"Initialising..."<<std::endl;
    std::vector<std::thread> iniThreads;
    for (int i = 0; i < total_simulations; i += number_of_threads) {
        iniThreads.clear();
        for (int j = i; j < std::min(i + number_of_threads, total_simulations); j++) {
            std::cout << "About to initialize simulation " << j << std::endl;
            iniThreads.emplace_back(initialize, j, init_config_number);
        }
        for (auto& t : iniThreads) t.join();
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::cout<<"Burning in..."<<std::endl;
    

    //Burn in  configurations//
    std::vector<std::thread> burnThreads;
    for (int i = 0; i < total_simulations; i += number_of_threads) {
        burnThreads.clear();
        for (int j = i; j < std::min(i + number_of_threads, total_simulations); j++) {
            burnThreads.emplace_back(burn_in, j, burn_in_steps);
        }
        for (auto& t : burnThreads) t.join();
    }

    std::cout<<"Start gradient descent..."<<std::endl;

    // Perform inverse algorithm //
    for (auto step = 0; step <= it_steps; step++) {
        
        
        //Forward run//
        for (int sim = 0; sim < total_simulations; sim += 30) {
            std::vector<std::thread> threads;
            for (int l = sim; l < std::min(sim + 30, total_simulations); l++) {
                threads.emplace_back(run, l, std::round(mc_moves * sqrt(step+1)));
            }
            for (auto& t : threads) t.join();
        }

        //MAKE SURE ALL STAGES ARE BEING RUN IN THE FORWARD RUN WITH COUTS IN MOVE FUNCTION

        //After all threads finish, their contact maps (total_contacts) are summed into a global matrix final_contacts//
        std::vector<double> w = {1., 1., 1., 1., 0.9270073, 0.35766423}; 
        double w_sum = std::accumulate(w.begin(), w.end(), 0.0);
        for (auto& wi : w) wi /= w_sum;  // Normalize to sum to 1


        for (int s = 0; s < stage_count; s++) {
            double weight = w[s] / replicates_per_stage;
            for (int r = 0; r < replicates_per_stage; r++) {
                int idx = s * replicates_per_stage + r;
                for (int i = 0; i < bin_num; i++) {
                    for (int j = 0; j < bin_num; j++) {
                        final_contacts[i][j] += weight * total_contacts[idx][i][j];
                    }
                }
            }
        }  
        
        //Ensures final_contacts is symmetric//
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j++) {
                final_contacts[j][i] = final_contacts[i][j];
            }
        }

        //Seems to save the configurations at the end of each iteration, ring= unreplicated/original part, strand= replicated part
        for (int l = 0; l < total_simulations; l++) {
            get_configuration(step, "strand", l);
            get_configuration(step, "ring", l);
        }

        //calculation of positional constraint values after forward run//FOR NOW I DONT CARE ABOUT THE ORI CONSTRAINT

        for (int l = 0; l < total_simulations; l++) { // JM: calculation of average z coordinates & squares of both oris
            double norm_factor = total_simulations / number_of_stages; // replicates_per_stage
            double normalization = 1.0 / (norm_factor * std::round(mc_moves * sqrt(step + 1)) * res);

            int stage_idx = l % number_of_stages;

            if (lin_length[stage_idx] == 0) {
                z_close_tot[stage_idx] += z_close[l] * normalization;
                z_close_squared_tot[stage_idx] += z_close_squared[l] * normalization;
            } else {
                z_close_tot[stage_idx] += z_close[l] * normalization;
                z_far_tot[stage_idx] += z_far[l] * normalization;
                z_close_squared_tot[stage_idx] += z_close_squared[l] * normalization;
                z_far_squared_tot[stage_idx] += z_far_squared[l] * normalization;
            }
        }

// Compute variance after summing across all simulations
for (int s = 0; s < number_of_stages; s++) {
    z_close_var[s] = z_close_squared_tot[s] - pow(z_close_tot[s], 2);
    z_far_var[s] = z_far_squared_tot[s] - pow(z_far_tot[s], 2);
}

        //update simulation parameters//
        normalize();
        update_E(step); //JM: Updates pairwise interaction energies

        normalize_measured_z(); //GG: express measured values of z in the units of cell length, and shift such that z=0 at the left pole
        update_energ_coeff_mean(step);
        update_energ_coeff_separation(step);
        //update_alpha_beta(step);

        //save simulation parameters//
        get_final_contacts(step);
        get_energ_coeff_mean(step);
        get_energ_coeff_separation(step);
        //get_alpha_beta(step); //JM: Saves new values alpha and beta
        get_energies_plot(step); //JM: Saves pairwise interaction energies

        clean_up(); //JM: sets everything to zero again before next round
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;
}


void run(int thread_num, int move_num) {
    int m = 0;
    while (m < move_num) {    //perform polymer simulation
        move(thread_num, m); //GG: "m" is passed by reference now
    }
    for (const auto& elem : contacts[thread_num]) { //add the contacts remaining at the end of the simulation (during a the simulation, a contact is only added to 'total_contacts' when a contact is removed)
        total_contacts[thread_num][std::min(elem.first.first,elem.first.second)][std::max(elem.first.first,elem.first.second)] += double(move_num) - elem.second;
    }
    for (const auto& elem : contacts_lin[thread_num]) { //add the contacts remaining at the end of the simulation (during a the simulation, a contact is only added to 'total_contacts' when a contact is removed)
        total_contacts[thread_num][std::min(elem.first.first,elem.first.second)][std::max(elem.first.first,elem.first.second)] += double(move_num) - elem.second;
    }
    for (const auto& elem : contacts_inter[thread_num]) { //add the contacts remaining at the end of the simulation (during a the simulation, a contact is only added to 'total_contacts' when a contact is removed)
        total_contacts[thread_num][std::min(elem.first.first,elem.first.second)][std::max(elem.first.first,elem.first.second)] += double(move_num) - elem.second;
    }
}

void normalize() {
    double sum = 0;
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j==bin_num - 1) && !(i==bin_num - 1 && j==0)) {
                sum += final_contacts[i][j];
            }
        }
    }
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            final_contacts[i][j] *= float(bin_num) / sum;
        }
    }
}

void update_E(int step) {
    for (int i = 0; i < bin_num; i++) {                 //update energy map
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j==bin_num - 1) && !(i==bin_num - 1 && j==0)) {
                auto comp = xp_contacts[i][j];
                if (comp != 0) {
                    if (step>update_cap_onset) { //JM we impose an upper limit on the update after the first few stages, to facilitate convergence
                        double prop_update = learning_rate * (final_contacts[i][j] - comp) / pow(std::abs(comp), 0.5);
                        if (abs(prop_update)<update_cap_factor * learning_rate){
                            Interaction_E[i][j] +=prop_update;
                        }
                        else{
                            if (prop_update<0){
                                Interaction_E[i][j] -=update_cap_factor * learning_rate;
                            }
                            else{
                                Interaction_E[i][j] +=update_cap_factor * learning_rate;
                            }
                        }
                    }
                    else{
                        Interaction_E[i][j] +=
                                learning_rate * (final_contacts[i][j] - comp) / pow(std::abs(comp), 0.5);
                    }
                }
            }
        }
    }
    double meanE = 0;                                   //mean energy
    double sum = 0;
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j==bin_num - 1) && !(i==bin_num - 1 && j==0)) {
                auto comp = xp_contacts[i][j];
                meanE += Interaction_E[i][j] * comp;
                sum += comp;
            }
        }
    }
    if (sum != 0) {
        meanE /= sum;
    }
    std::cout << std::setw(4) << std::left << step << std::setw(10)  << std::setprecision(3) << std::right << meanE << "\n";
    for (int i = 0; i < bin_num; i++) {                 //set average energy to zero
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j==bin_num - 1) && !(i==bin_num - 1 && j==0)) {
                Interaction_E[i][j] -= meanE;
            }
        }
    }
}

void normalize_measured_z(){ // express measured values of z in the units of cell length, also shifts such that z=0 at the left pole
    for (int s=0; s<number_of_stages; s++){
        double total_cell_length = length[s] + 2*radius;
        double left_pole_z = -(length[s]/2 + radius - offset_z[s]);

        z_close_tot[s] = (z_close_tot[s] - left_pole_z)/total_cell_length;
        z_far_tot[s] = (z_far_tot[s] - left_pole_z)/total_cell_length;
        z_close_var[s] = z_close_var[s] / pow(total_cell_length, 2);
        z_far_var[s] = z_far_var[s] / pow(total_cell_length, 2);
        for(int i=0; i<sites_constrained_mean.size(); i++){
            z_mean_data_tot[s][i] = (z_mean_data_tot[s][i] - left_pole_z) / total_cell_length;
        }
        for(int i=0; i<sites_constrained_separation.size(); i++){
            z_separation_data_tot[s][i] = z_separation_data_tot[s][i] / total_cell_length;
        }
    }
}

void update_alpha_beta (int step) {
    for (int s = 0; s < number_of_stages; s++) {
        if(is_constrained_ori[s][0]){
            alpha[s] += learning_rate_close*(z_close_tot[s]-xp_z_close[s]);
        }
        if(is_constrained_ori[s][2]){
            alpha2[s] += learning_rate_close_var * (z_close_var[s] - xp_z_close_var[s]);
        }
        if(is_replicated[s][oriC]){
            if(is_constrained_ori[s][1]){
                beta[s] += learning_rate_far * (z_far_tot[s]-xp_z_far[s]);
            }
            if(is_constrained_ori[s][3]){
                beta2[s] += learning_rate_far_var * (z_far_var[s]- xp_z_far_var[s]);
            }
        }
    }
}

void update_energ_coeff_mean (int step) {
    for(int stage=0; stage<number_of_stages; stage++) {
        for (int i = 0; i < sites_constrained_mean.size(); i++) {
            if(is_constrained_mean[stage][i]){
                energ_coeff_mean[stage][i] += learning_rate_means * (z_mean_data_tot[stage][i] - target_means[stage][i]);
            }
        }
    }
}

void update_energ_coeff_separation (int step) {
    for(int stage=0; stage<number_of_stages; stage++) {
        for (int i = 0; i < sites_constrained_separation.size(); i++) {
            if(is_replicated[stage][sites_constrained_separation[i]] && is_constrained_separation[stage][i]) {
                energ_coeff_separation[stage][i] += learning_rate_separations * (z_separation_data_tot[stage][i] - target_separations[stage][i]);
            }
        }
    }
}

void set_unconstrained_energies_to_zero() { //for safety: makes sure that the initial energies with no corresponding experimental constraints are 0
    for (int s=0; s<number_of_stages; s++){
        if(!is_constrained_ori[s][0]){
            alpha[s] = 0;
        }
        if(!is_constrained_ori[s][1]){
            beta[s] = 0;
        }
        if(!is_constrained_ori[s][2]){
            alpha2[s] = 0;
        }
        if(!is_constrained_ori[s][3]){
            beta2[s] = 0;
        }
        for(int i=0; i<n_constrained_mean; i++) {
            if (!is_constrained_mean[s][i]) {
                energ_coeff_mean[s][i] = 0;
            }
        }
        for(int i=0; i<n_constrained_separation; i++) {
            if(!is_constrained_separation[s][i]){
                energ_coeff_separation[s][i] = 0;
            }
        }
    }
}

void clean_up() {
    std::vector<double> zeroVec(bin_num, 0);
    std::fill(final_contacts.begin(), final_contacts.end(), zeroVec);
    for (int l = 0; l < total_simulations; l++) {
    //CB change 
    //for (int l = 0; l < total_simulations; l++) {
        std::fill(total_contacts[l].begin(), total_contacts[l].end(), zeroVec);
        z_close[l] = 0;
        z_far[l] = 0;
        z_close_tot[l%number_of_stages]=0;
        z_far_tot[l%number_of_stages]=0;
        z_close_squared[l] = 0;
        z_far_squared[l] = 0;
        z_close_var[l%number_of_stages] = 0;
        z_far_var[l%number_of_stages] = 0;
        z_close_squared_tot[l%number_of_stages]=0;
        z_far_squared_tot[l%number_of_stages]=0;

        for(int i=0; i<sites_constrained_mean.size(); i++){
            z_mean_data[l][i] = 0;
            z_mean_data_tot[l%number_of_stages][i] = 0;
        }
        for(int i=0; i<sites_constrained_separation.size(); i++){
            z_separation_data[l][i] = 0;
            z_separation_data_tot[l%number_of_stages][i] = 0;
        }

        contacts[l].clear();
        contacts_lin[l].clear();
        contacts_inter[l].clear();
        for (int i = 0; i < pol_length; i+=reduction_factor) {
            int red_i=i/reduction_factor;
            if (locations[l].find({polymer[l][i][0], polymer[l][i][1], polymer[l][i][2]}) != locations[l].end()) {
                for (auto elem : locations[l][{polymer[l][i][0], polymer[l][i][1],polymer[l][i][2]}]) {
                    contacts[l][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                }
            }
        }
        for (int i = 0; i < pol_length; i+=reduction_factor) {
            int red_i=i/reduction_factor;
            if (oriC + lin_length[l%number_of_stages] >= pol_length) {
                if (i > oriC - lin_length[l%number_of_stages] || i < (oriC + lin_length[l%number_of_stages]) % pol_length) {
                    if (locations_lin[l].find({lin_polymer[l][i][0], lin_polymer[l][i][1],lin_polymer[l][i][2]}) != locations_lin[l].end()) {
                        for (auto elem : locations_lin[l][{lin_polymer[l][i][0],lin_polymer[l][i][1],lin_polymer[l][i][2]}]) {
                            contacts_lin[l][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                        }
                    }
                    if (locations[l].find({lin_polymer[l][i][0], lin_polymer[l][i][1],lin_polymer[l][i][2]}) !=
                        locations[l].end()) {
                        for (auto elem : locations[l][{lin_polymer[l][i][0],lin_polymer[l][i][1],lin_polymer[l][i][2]}]) {
                            contacts_inter[l][{elem, red_i}] = 0;
                        }
                    }
                }
            } else {
                if (i > oriC - lin_length[l%number_of_stages] && i < oriC + lin_length[l%number_of_stages]) {
                    if (locations_lin[l].find({lin_polymer[l][i][0], lin_polymer[l][i][1],lin_polymer[l][i][2]}) != locations_lin[l].end()) {
                        for (auto elem : locations_lin[l][{lin_polymer[l][i][0],lin_polymer[l][i][1],lin_polymer[l][i][2]}]) {
                            contacts_lin[l][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                        }
                    }
                    if (locations[l].find({lin_polymer[l][i][0], lin_polymer[l][i][1],lin_polymer[l][i][2]}) !=
                        locations[l].end()) {
                        for (auto elem : locations[l][{lin_polymer[l][i][0],lin_polymer[l][i][1],lin_polymer[l][i][2]}]) {
                            contacts_inter[l][{elem, red_i}] = 0;
                        }
                    }
                }
            }
        }
    }
}
