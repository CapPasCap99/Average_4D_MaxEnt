#include <fstream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include "global.h"
using namespace Eigen;

void read_file(std::vector<std::vector<double>>& map, std::string& filename) {

    std::cout<<"Reading in "<<filename<<std::endl;
    std::ifstream couplings;
    couplings.open(filename);
    if(couplings.fail()){throw std::invalid_argument( "Hi-C input file doesnt exist." );}
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            couplings >> map[i][j];
        }
    }
    couplings.close();
}

void read_constraints(std::vector<double> &target_vector, std::string& filename){
    std::ifstream values;
    std::cout<<"Reading in "<<filename<<std::endl;
    values.open(filename);
    for (int i =0;i<target_vector.size();i++){
        values>>target_vector[i];
    }
    values.close();
}

void read_interaction_energies(std::string& old_simulation_folder, std::string& iteration) {
    std::ifstream couplings;
    
    std::cout<<"Reading in "<<dir + old_simulation_folder + "Energies/energies_" + iteration + ".txt"<<std::endl;
    couplings.open(dir + old_simulation_folder + "Energies/energies_" + iteration + ".txt");  //read in energies
    if(couplings.fail()){throw std::invalid_argument( "Interaction energies input file doesnt exist." );}
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            couplings >> Interaction_E[i][j];
        }
    }
    couplings.close();
}

void check_input_energies_compatibility(std::string& old_simulation_folder, std::string& iteration){
    std::ifstream input_file;
    std::string word = " ";
    std::cout<<"Reading in "<<dir + old_simulation_folder + "sim_params.txt"<<std::endl;
    input_file.open(dir + old_simulation_folder + "sim_params.txt");
    if(input_file.fail()){throw std::invalid_argument( "Input energy file doesnt exist (or sim_params.txt missing)." );}
    while(word != "Stages:"){
        input_file >> word;
    }
    for (int s : stages){
        input_file >> word;
        if(word.find_first_not_of( "0123456789" ) != std::string::npos || std::stoi(word) != s ){ // check if "word" is an int
            throw std::invalid_argument( "Incompatible energy input stages." );
        }
    }
}

void check_input_configuration_compatibility(std::string& old_simulation_folder){
    std::ifstream input_file;
    std::string word = " ";
    input_file.open(dir + old_simulation_folder + "sim_params.txt");
    if(input_file.fail()){throw std::invalid_argument( "Input configuration file doesnt exist (or sim_params.txt missing)." );}
    while(word != "Thread"){
        input_file >> word;
    }
    input_file >> word;
    input_file >> word; // word = number of threads
    if(number_of_threads > std::stoi(word)){
        throw std::invalid_argument( "Incompatible configuration input (not enough configurations)." );
    }

    while(word != "Stages:"){
        input_file >> word;
    }
    for (int s : stages){
        input_file >> word;
        if(word.find_first_not_of( "0123456789" ) != std::string::npos || std::stoi(word) != s ){ // check if "word" is an int
            throw std::invalid_argument( "Incompatible configuration input stages." );
        }
    }
}

void read_ori_energies(std::string& old_simulation_folder, std::string& iteration) {
    std::ifstream energies_file;
    energies_file.open(dir + old_simulation_folder + "Energies/alpha_" + iteration + ".txt");
    if(energies_file.fail()){throw std::invalid_argument( "Alpha energies file missing." );}
    for (int s=0; s<number_of_stages; s++){
        energies_file >> alpha[s];
    }
    energies_file.close();

    energies_file.open(dir + old_simulation_folder + "Energies/beta_" + iteration + ".txt");
    if(energies_file.fail()){throw std::invalid_argument( "Beta energies file missing." );}
    for (int s=0; s<number_of_stages; s++){
        energies_file >> beta[s];
    }
    energies_file.close();

    energies_file.open(dir + old_simulation_folder + "Energies/alpha2_" + iteration + ".txt");
    if(energies_file.fail()){throw std::invalid_argument( "Alpha2 energies file missing." );}
    for (int s=0; s<number_of_stages; s++){
        energies_file >> alpha2[s];
    }
    energies_file.close();

    energies_file.open(dir + old_simulation_folder + "Energies/beta2_" + iteration + ".txt");
    if(energies_file.fail()){throw std::invalid_argument( "Beta2 energies file missing." );}
    for (int s=0; s<number_of_stages; s++){
        energies_file >> beta2[s];
    }
    energies_file.close();
}

void read_position_energies(std::string& old_simulation_folder, std::string& iteration) {
    for (int i=0; i<sites_constrained_mean.size(); i++){
        std::ifstream energies_file;
        energies_file.open(dir + old_simulation_folder + "Energies/pos_energ_" + iteration + "_site" + std::to_string(sites_constrained_mean[i]) + ".txt");
        if(energies_file.fail()){throw std::invalid_argument( "Missing position energy file for site " + std::to_string(sites_constrained_mean[i]) + "." );}
        for (int s=0; s<number_of_stages; s++){
            energies_file >> energ_coeff_mean[s][i];
        }
        energies_file.close();
    }
}

void read_separation_energies(std::string& old_simulation_folder, std::string& iteration) {
    for (int i=0; i<sites_constrained_separation.size(); i++){
        std::ifstream energies_file;
        energies_file.open(dir + old_simulation_folder + "Energies/sep_energ_" + iteration + "_site" + std::to_string(sites_constrained_separation[i]) + ".txt");
        if(energies_file.fail()){throw std::invalid_argument( "Missing separation energy file for site " + std::to_string(sites_constrained_separation[i]) + "." );}
        for (int s=0; s<number_of_stages; s++){
            energies_file >> energ_coeff_separation[s][i];
        }
        energies_file.close();
    }
}

void read_configuration(std::string name, int thread_num, int step) {
    std::ostringstream filename; filename << dir << configuration_data_folder << "Configurations/configuration_" << name << "_" << step << "_" << thread_num << ".txt";
    std::ifstream monomers;
    monomers.open(filename.str());  //read in monomer positions
    if(monomers.fail()){throw std::invalid_argument( "configuration input file doesnt exist" );}

    if (name == "ring") {
        for (int i = 0; i < pol_length; i++) {
            Vector3i location = {0,0,0};
            for (int j = 0; j < 3; j++) {
                monomers >> location[j];
            }
            polymer[thread_num].push_back(location);
        }
    }
    else if (name == "strand") {
        for (int i = 0; i < pol_length; i++) {
            Vector3i location = {0,0,0};
            for (int j = 0; j < 3; j++) {
                monomers >> location[j];
            }
            lin_polymer[thread_num].push_back(location);
        }
    }
    monomers.close();
}

void read_input_data(){
    std::ifstream input_file;
    std::string word;
    int site;
    float exp_data;
    for(int s=0; s<number_of_stages; s++) {
        std::ostringstream filename;
        filename << dir << "Input/" << bacteria_name << "/" << bacteria_name << "_stage_" << stages[s] << ".txt";
        input_file.open(filename.str());
        if(input_file.fail()){throw std::invalid_argument( "Missing input data for "+ bacteria_name + " at stage " + std::to_string(stages[s]) +"." );}
        while (!input_file.eof()) {
            input_file >> word;
            if (word == "lin_length") {
                input_file >> lin_length[s];
            }
            else if(word == "length"){
                input_file >> length[s];
            }
            else if(word == "ori_pos_near"){
                input_file >> xp_z_close[s];
                is_constrained_ori[s][0] = true;
            }
            else if(word == "ori_pos_far"){
                input_file >> xp_z_far[s];
                is_constrained_ori[s][1] = true;
            }
            else if(word == "ori_var_near"){
                input_file >> xp_z_close_var[s];
                is_constrained_ori[s][2] = true;
            }
            else if(word == "ori_var_far"){
                input_file >> xp_z_far_var[s];
                is_constrained_ori[s][3] = true;
            }
            else if(word == "mean_positions"){
                input_file >> word;
                while(word != "end"){
                    input_file >> site;
                    input_file >> exp_data;
                    if(sites_constrained_mean_map.find(site) != sites_constrained_mean_map.end() ){ //check if site constrained
                        target_means[s][ sites_constrained_mean_map[site] ] = exp_data;
                        is_constrained_mean[s][ sites_constrained_mean_map[site] ] = true;
                    }
                    input_file >> word;
                }
            }
            else if(word == "mean_separations"){
                input_file >> word;
                while(word != "end"){
                    input_file >> site;
                    input_file >> exp_data;
                    if(sites_constrained_separation_map.find(site) != sites_constrained_separation_map.end() ){ //check if site constrained
                        target_separations[s][ sites_constrained_separation_map[site] ] = exp_data;
                        is_constrained_separation[s][ sites_constrained_separation_map[site] ] = true;
                    }
                    input_file >> word;
                }
            }
        }
        input_file.close();

        // GG: extend vectors to length "number_of_threads" (needed for the distribution functions and maybe something else from the old forward code)
        for(int i=0; i<number_of_threads; i++){
            length[i] = length[i % number_of_stages];
            lin_length[i] = lin_length[i % number_of_stages];
        }
    }
}


void read_fork_distribution(std::string fork_distribution_file){
    std::ifstream input_file;
    int multiplicity;
    int fork_pos;
    int counter = 0;

    input_file.open( dir + "Input/fork_distributions/" + fork_distribution_file + "_stage_" + std::to_string(stages[0]) + ".txt" );
    if(input_file.fail()){throw std::invalid_argument( "Missing fork distribution data for stage " + std::to_string(stages[0]) +"." );}
    while (!input_file.eof()) {
        input_file >> multiplicity;
        input_file >> fork_pos;
        if(multiplicity<0){
            break;
        }
        for(int i=0; i<multiplicity; i++){
            if(counter==number_of_threads){ throw std::invalid_argument( "Fork distribution doesn't match the thread number!"); }
            lin_length[counter] = std::max(fork_pos, 0);
            counter++;
        }
    }
    if(counter<number_of_threads){ throw std::invalid_argument( "Fork distribution doesn't match the thread number!"); }
}


void get_energies_plot(const int& steps) { //saves energies as matrix
    std::ostringstream fn;
    fn << dir << output_folder << "/" << "Energies/energies_" << steps << ".txt";
    std::ofstream final_energies;
    final_energies.open(fn.str().c_str(), std::ios_base::binary);

    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            double contact;
            contact = Interaction_E[i][j];
            final_energies << contact << ' ';
        }
        final_energies << '\n';
    }
    final_energies.close();
}

void get_final_contacts(const int& step) {
    std::ostringstream fn;
    fn << dir << output_folder << "/" << "Contacts/contacts_" << step << ".txt";
    std::ofstream final_cont;
    final_cont.open(fn.str().c_str(), std::ios_base::binary); //write contact frequencies
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            double contact;
            if (std::abs(i-j) <= 1 || (i==0 && j == bin_num - 1) || (j==0 && i == bin_num - 1)) { contact = 0; }
            else {
                contact = final_contacts[i][j];
            }
            final_cont << contact << ' ';
        }
        final_cont << '\n';
    }
}

void get_contacts_xp() {
    std::ostringstream fn;
    fn << dir << "contacts_xp.txt";
    std::ofstream final_cont;
    final_cont.open(fn.str().c_str(), std::ios_base::binary); //write contact frequencies
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            double contact;
            if (std::abs(i-j) <= 1 || (i==0 && j == bin_num - 1) || (j==0 && i == bin_num - 1)) { contact = 0; }
            else {
                contact = xp_contacts[i][j];
            }
            final_cont << contact << ' ';
        }
        final_cont << '\n';
    }
}

void get_configuration(int step, std::string name, int thread_num) {
    std::ostringstream fn;
    fn << dir << output_folder << "/" << "Configurations/configuration_" << name << "_" << step << "_" << thread_num << ".txt";
    std::ofstream final_conf;
    final_conf.open(fn.str().c_str(), std::ios_base::binary);
    if (name == "ring") {
        for (int i = 0; i < pol_length; i++) {
            final_conf << polymer[thread_num][i].transpose() << '\n';
        }
    }
    else if (name == "strand"){
        for (int i = 0; i < pol_length; i++) {
            final_conf << lin_polymer[thread_num][i].transpose() << '\n';
        }
    }
    final_conf.close();
    std::ostringstream fn_old;
    fn_old << dir << output_folder << "/" << "Configurations/configuration_" << name << "_" << step - 1 << "_" << thread_num << ".txt";
    std::remove(fn_old.str().c_str());
}

void get_alpha_beta(int step) {
    std::ostringstream fn;
    fn << dir << output_folder << "/" << "Energies/alpha_" << step << ".txt";
    std::ofstream values;
    values.open(fn.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        values << std::setw(7) << std::setprecision(4) << alpha[s] << " ";
    }
    values << '\n';
    values.close();
    std::ostringstream fn2;
    fn2 << dir << output_folder << "/" << "Energies/beta_" << step << ".txt";
    std::ofstream values2;
    values2.open(fn2.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        values2 << std::setw(7) << std::setprecision(4) << beta[s] << " ";
    }
    values2 << '\n';
    values2.close();
    std::ostringstream fn3;
    fn3 << dir << output_folder << "/" << "Energies/beta2_" << step << ".txt";
    std::ofstream values3;
    values3.open(fn3.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        values3 << std::setw(7) << std::setprecision(4) << beta2[s] << " ";
    }
    values3 << '\n';
    values3.close();
    std::ostringstream fn4;
    fn4 << dir << output_folder << "/" << "Energies/alpha2_" << step << ".txt";
    std::ofstream values4;
    values4.open(fn4.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        values4 << std::setw(7) << std::setprecision(4) << alpha2[s] << " ";
    }
    values4 << '\n';
    values4.close();
}

void get_z_lin_far_close(int step) {
    std::ostringstream fn;
    fn << dir << output_folder << "/" << "Positions/close_" << step << ".txt";
    std::ofstream values;
    values.open(fn.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        values << std::setw(5) << std::setprecision(3) << z_close_tot[s] << " ";
    }
    values << '\n';
    values.close();
    std::ostringstream fn2;
    fn2 << dir << output_folder << "/" << "Positions/far_" << step << ".txt";
    std::ofstream values2;
    values2.open(fn2.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        values2 << std::setw(5) << std::setprecision(3) << z_far_tot[s] << " ";
    }
    values2 << '\n';
    values2.close();
    std::ostringstream fn3;
    fn3 << dir << output_folder << "/" << "Positions/close_var_" << step << ".txt";
    std::ofstream values3;
    values3.open(fn3.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        values3 << std::setw(5) << std::setprecision(3) << z_close_var[s] << " ";
    }
    values3 << '\n';
    values3.close();
    std::ostringstream fn4;
    fn4 << dir << output_folder << "/" << "Positions/far_var_" << step << ".txt";
    std::ofstream values4;
    values4.open(fn4.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        values4 << std::setw(5) << std::setprecision(3) << z_far_var[s] << " ";
    }
    values4 << '\n';
    values4.close();
}

void get_z_mean_rest(int site_index, int step) {
    std::ostringstream fn;
    fn << dir << output_folder << "/" << "Positions/means_" << step << "_site" + std::to_string(sites_constrained_mean[site_index])+ ".txt";
    std::ofstream values;
    values.open(fn.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        if (s<(number_of_stages-1)) {
            values << std::setprecision(3) << z_mean_data_tot[s][site_index] << " ";
        }
        else{
            values << std::setprecision(3) << z_mean_data_tot[s][site_index];
        }
    }
    values.close();
}

void get_z_separation_rest(int site_index, int step) {
    std::ostringstream fn;
    fn << dir << output_folder << "/" << "Positions/separations_" << step << "_site" + std::to_string(sites_constrained_separation[site_index])+ ".txt";
    std::ofstream values;
    values.open(fn.str().c_str(), std::ios_base::app);
    for (int s = 0; s < number_of_stages; s++) {
        if (s<(number_of_stages-1)) {
            values << std::setprecision(3) << z_separation_data_tot[s][site_index] << " ";
        }
        else{
            values << std::setprecision(3) << z_separation_data_tot[s][site_index];
        }
    }
    values.close();
}


void get_energ_coeff_mean(int step){
    for (int i = 0; i < sites_constrained_mean.size(); i++) {
        std::ostringstream fn;
        fn << dir << output_folder << "/" << "Energies/pos_energ_" << step << "_site" << std::to_string(sites_constrained_mean[i]) << ".txt";
        std::ofstream values;
        values.open(fn.str().c_str(), std::ios_base::app);
        for (int s = 0; s < number_of_stages; s++) {
            values << std::setprecision(4) << energ_coeff_mean[s][i] << " ";
        }
        values.close();
    }
}

void get_energ_coeff_separation(int step){
    for (int i = 0; i < sites_constrained_separation.size(); i++) {
        std::ostringstream fn;
        fn << dir << output_folder << "/" << "Energies/sep_energ_" << step << "_site" << std::to_string(sites_constrained_separation[i]) << ".txt";
        std::ofstream values;
        values.open(fn.str().c_str(), std::ios_base::app);
        for (int s = 0; s < number_of_stages; s++) {
            if (s<(number_of_stages-1)) {
                values << std::setprecision(4) << energ_coeff_separation[s][i] << " ";
            }
            else{
                values << std::setprecision(4) << energ_coeff_separation[s][i];
            }
        }
        values.close();
    }
}

void get_sim_params() {
    std::ostringstream fn;
    fn << dir << output_folder << "/" << "sim_params.txt";
    std::ofstream params;
    params.open(fn.str().c_str(), std::ios_base::binary);
    params << "Thread number: " << number_of_threads << '\n';
    params << "Bin number: " << bin_num << '\n';
    params << "Polymer length: " << pol_length << '\n';
    params << "MC moves: " << mc_moves << '\n';
    params << "Bacteria: " << bacteria_name << '\n';

    params << "Stages: ";
    for (auto el: stages) {
        params << el << " ";
    }
    params << '\n';
    if(initConfig) {
        params << "Input configurations folder: " << configuration_data_folder << '\n';
    }
    params << "Experimental HiC data file: " << HiC_file << '\n';

    params << '\n';
    params << "Learning rates:" << '\n';
    params << "interactions: " << learning_rate << '\n';
    params << "near ori mean: " << learning_rate_close << '\n';
    params << "far ori mean: " << learning_rate_far << '\n';
    params << "near ori variance: " << learning_rate_close_var << '\n';
    params << "far ori variance: " << learning_rate_far_var << '\n';
    params << "other sites means: " << learning_rate_means << '\n';
    params << "other sites separations: " << learning_rate_separations << '\n';

    params << '\n';
    params << "Confinement/cell shape:" << '\n';
    params << "cell radius: " << radius << '\n';
    params << "cell lengths: ";
    // for (int s = 0; s < number_of_stages; s++) {
    //     params << std::setw(3) << length[s] << " ";
    // }
    for (int s = 0; s < number_of_stages; s++) {
        for (int r = 0; r < replicates_per_stage; r++) {
            int idx = s * replicates_per_stage + r;
            params << std::setw(3) << length[idx] << " ";
        }
        params << '\n';
        if (!use_fork_distribution) break;
    }
    params << '\n';
    params << "offset: " << offset[0] << " "<< offset[1] << '\n';
    params << "oriC: " << oriC << '\n';
    params << "lin_lengths: ";
    // for (int s = 0; s < number_of_threads; s++) {
    //     params << std::setw(3) << lin_length[s] << " ";
    //     if(not use_fork_distribution and s==number_of_stages-1){break;}
    // }
    for (int s = 0; s < number_of_stages; s++) {
        for (int r = 0; r < replicates_per_stage; r++) {
            int idx = s * replicates_per_stage + r;
            params << std::setw(3) << lin_length[idx] << " ";
        }
        params << '\n';
        if (!use_fork_distribution) break;
    }

    params << '\n';

    params <<'\n';
    params <<"Constraints on the ori:";
    params << '\n';
    params << std::setw(19) << "Near ori mean: ";
    for (int s=0; s<number_of_stages; s++){
        if(is_constrained_ori[s][0]){
            params << std::setw(5) << xp_z_close[s] << " ";
        }
        else{
            params << std::setw(5) << "*" << " ";
        }
    }
    params << '\n';
    params << std::setw(19) << "Far ori mean: ";
    for (int s=0; s<number_of_stages; s++){
        if(is_constrained_ori[s][1]){
            params << std::setw(5) << xp_z_far[s] << " ";
        }
        else{
            params << std::setw(5) << "*" << " ";
        }
    }
    params << '\n';
    params << std::setw(19) << "Near ori variance: ";
    for (int s=0; s<number_of_stages; s++){
        if(is_constrained_ori[s][2]){
            params << std::setw(5) << xp_z_close_var[s] << " ";
        }
        else{
            params << std::setw(5) << "*" << " ";
        }
    }
    params << '\n';
    params << std::setw(19) << "Far ori variance: ";
    for (int s=0; s<number_of_stages; s++){
        if(is_constrained_ori[s][3]){
            params << std::setw(5) << xp_z_far_var[s] << " ";
        }
        else{
            params << std::setw(5) << "*" << " ";
        }
    }
    params << '\n';

    params << '\n';
    params << "Constraints on the mean:";
    params << '\n';
    for (int i=0; i<sites_constrained_mean.size(); i++){
        params << "site " << std::setw(3) << sites_constrained_mean[i] << ": ";
        for (int s=0; s<number_of_stages; s++) {
            if(is_constrained_mean[s][i]) {
                params << std::setw(5) << target_means[s][i] << " ";
            }
            else{
                params << std::setw(5) << "*" << " ";
            }
        }
        params << '\n';
    }

    params << '\n';
    params << "Constraints on the separation:";
    params << '\n';
    for (int i=0; i<sites_constrained_separation.size(); i++){
        params << "site " << std::setw(3) << sites_constrained_separation[i] << ": ";
        for (int s=0; s<number_of_stages; s++) {
            if(is_constrained_separation[s][i]) {
                params << std::setw(5) << target_separations[s][i] << " ";
            }
            else{
                params << std::setw(5) << "*" << " ";
            }
        }
        params << '\n';
    }

}
