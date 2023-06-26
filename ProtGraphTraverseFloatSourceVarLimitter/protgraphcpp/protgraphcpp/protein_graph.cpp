#include "protein_graph.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
//DEBUG and time measurement
#include <inttypes.h>
#include <chrono>

#include <ulimit.h>

// Constructor
ProteinGraph::ProteinGraph(std::uint32_t num_acc, std::ifstream &input, std::unordered_map<std::string, std::vector<uint8_t>> max_vars) {
    uint32_t num_n, num_e, num_pdbs, cur_32bit;
    uint8_t cur_8bit;
    uint16_t cur_16bit;
    uint64_t cur_64bit;
    std::string cur_string;

    // Get next entry and its global information about the entry
    num_acc = be32toh(num_acc);
    input.read((char*) &num_n, 4); num_n = be32toh(num_n);
    input.read((char*) &num_e, 4); num_e = be32toh(num_e);
    input.read((char*) &num_pdbs, 4); num_pdbs = be32toh(num_pdbs);

    // Set the most important Protein-Graph-Parameters
    this->N = num_n;
    this->E = num_e;
    this->PDB = num_pdbs;

    // Read AC
    for (int i = 0; i < num_acc; i++) {
        std::getline(input, cur_string, '\0');
        this->accessions.push_back(cur_string);
    }

    // Read NO
    this->nodes = new uint32_t[this->N];
    for (int i = 0; i < num_n; i++) {
        input.read((char*) &cur_32bit, 4); cur_32bit = be32toh(cur_32bit);
        this->nodes[i] = cur_32bit;
    }

    // Read ED
    this->edges = new uint32_t[this->E];
    for (int i = 0; i < num_e; i++) {
        input.read((char*) &cur_32bit, 4); cur_32bit = be32toh(cur_32bit);
        this->edges[i] = cur_32bit;
    }

    // Read SQ
    // Get mapping first
    std::unordered_map<std::string, std::vector<uint32_t>> seq_map;
    for (int i = 0; i < num_n; i++) {
        std::getline(input, cur_string, '\0');
        seq_map[cur_string].push_back(i);
    }
    
    cur_32bit = 0;
    for (auto const& [key, val] : seq_map) {
        cur_32bit += key.length() + 1; //? Chararray length
    }

    // Initialize via map
    this->sequence_str = new char[cur_32bit];
    cur_32bit = 0;
    this->sequence_str_index = new std::uint32_t[this->N];
    for (auto const& [key, val] : seq_map) {
        std::strcpy(&this->sequence_str[cur_32bit], key.c_str());
        // Now go through each entry in vector
        for (uint32_t node_idx: val) {
            this->sequence_str_index[node_idx] = cur_32bit;
        }
        cur_32bit += key.length() + 1;
        this->sequence_str[cur_32bit - 1] = 0;        
    }

    // Read PO
    this->position = new uint16_t[this->N];
    for (int i = 0; i < num_n; i++) {
        input.read((char*) &cur_16bit, 2); cur_16bit = be16toh(cur_16bit);
        this->position[i] = cur_16bit;
    }
    
    // Read IS
    this->iso_index = new uint8_t[this->N];
    for (int i = 0; i < num_n; i++) {
        input.read((char*) &cur_8bit, 1);
        this->iso_index[i] = cur_8bit;
    }

    // Read IP
    this->iso_position = new uint16_t[this->N];
    for (int i = 0; i < num_n; i++) {
        input.read((char*) &cur_16bit, 2); cur_16bit = be16toh(cur_16bit);
        this->iso_position[i] = cur_16bit;
    }

    // Read MW
    this->mono_weight = new double[this->N];
    for (int i = 0; i < num_n; i++) {
        input.read((char*) &cur_64bit, 8); cur_64bit = be64toh(cur_64bit);
        this->mono_weight[i] = (double)cur_64bit;
    }

    // Read CL
    this->cleaved = std::vector<bool>(E, false);
    for (int i = 0; i < num_e; i++) {
        input.read((char*) &cur_8bit, 1);
        if (cur_8bit != 0) {
            this->cleaved.at(i).flip();
        }
    }

    // Read QU
    // Get mapping first
    std::unordered_map<std::string, std::vector<uint32_t>> qu_map;
    for (int i = 0; i < num_e; i++) {
        std::getline(input, cur_string, '\0');
        qu_map[cur_string].push_back(i);
    }

    cur_32bit = 0;
    for (auto const& [key, val] : qu_map) {
        cur_32bit += key.length() + 1; //? Chararray length
    }

    // Initialize via map
    this->qualifiers_str = new char[cur_32bit];
    cur_32bit = 0;
    this->qualifiers_str_index = new std::uint32_t[this->E];
    for (auto const& [key, val] : qu_map) {
        std::strcpy(&this->qualifiers_str[cur_32bit], key.c_str());
        // Now go through each entry in vector
        for (uint32_t edge_idx: val) {
            this->qualifiers_str_index[edge_idx] = cur_32bit;
        }
        cur_32bit += key.length() + 1;
        this->qualifiers_str[cur_32bit - 1] = 0;        
    }

    // Read VC
    this->variant_count = new std::uint8_t[this->E];
    for (int i = 0; i < num_e; i++) {
        input.read((char*) &cur_8bit, 1);
        this->variant_count[i] = cur_8bit;
    }

    // Read PDBs (sequentially in RAM)
    this->pdbs = new double[this->N*2*this->PDB];
    for (int i=0; i<num_n; i++) {
        for (int j=0; j<2*num_pdbs; j++) {
            input.read((char*) &cur_64bit, 8); cur_64bit = be64toh(cur_64bit);
            if (cur_64bit != uint64_t(-1)) {
                this->pdbs[i*2*num_pdbs + j] = (double)cur_64bit;
            } else {
                this->pdbs[i*2*num_pdbs + j] = (double)INT64_MAX;
            }
        }
    }
    
    // Lastly Read max_vars into the proteingraph
    this->max_vars_bins = new std::uint8_t[max_vars[this->accessions[0]].size()];
    this->num_bins = max_vars[this->accessions[0]].size();
    for (int i=0; i < this->num_bins; i++) {
        this->max_vars_bins[i] = max_vars[this->accessions[0]].at(i);
    }

}



// Gets the edge from the source and target node. NOTE: the edge should exist, otherwise it returns 0!
uint32_t ProteinGraph::get_edge_index(uint32_t source_node, uint32_t target_node) {
    // if we are at the start node, we set it to 0
    uint32_t start_idx = (source_node != 0) ? this->nodes[source_node-1] : 0;
    for (uint32_t idx = start_idx; idx <= this->nodes[source_node]; idx++) {
        if (this->edges[idx] == target_node) {
            return idx;
        }
    }
    return 0;
}



std::string ProteinGraph::convert_paths_to_fasta(std::vector<std::vector<uint32_t>> paths) {
    std::string sequence, sequence_intermediate, sequence_new_lined;       // by nodes
    std::string qualifiers;   // by edges
    uint8_t iso_idx;  // by nodes
    uint32_t mssclvg, cur_edge_id;
    std::string spos, epos;
    int64_t tot_weight;
    bool spos_retrieved;

    std::string fasta = "";


    for (std::vector<uint32_t> path: paths) {
        qualifiers = "";
        sequence = "";
        sequence_new_lined = "";
        iso_idx = 0;
        mssclvg = 0;
        tot_weight = 0;
        spos_retrieved = false;

        // Get sequence, spos, mssclvg, iso_idx and qualifiers
        for (uint32_t idx = 1; idx < path.size()-1; idx++) {
            sequence_intermediate = std::string(&this->sequence_str[this->sequence_str_index[path.at(idx)]]);
            sequence.append(sequence_intermediate);  // concat sequences
            iso_idx = std::max(this->iso_index[path.at(idx)], iso_idx); // get the accession (maybe iso accession)

            if (!spos_retrieved && 0 != sequence_intermediate.compare("")) {
                // Considering n term modifications (if applicable)
                // Get starting pos
                spos_retrieved = true;
                if (this->iso_position[path.at(idx)] != UINT16_MAX){
                    spos = std::to_string(this->iso_position[path.at(idx)]);
                } else if (this->position[path.at(idx)] != UINT16_MAX) {
                    spos = std::to_string(this->position[path.at(idx)]);
                } else {
                    spos = '?';
                }
            }
            
            // We need to look up via the specific edge index!!!
            cur_edge_id = get_edge_index(path.at(idx-1), path.at(idx));
            if (this->cleaved.at(cur_edge_id)) {
                mssclvg++; // count misscleavages
            }
            // We need to look up via the specific edge index!!!
            if (std::string(&this->qualifiers_str[this->qualifiers_str_index[cur_edge_id]]).length() > 0) {
                qualifiers.append(std::string(&this->qualifiers_str[this->qualifiers_str_index[cur_edge_id]]));
                qualifiers.append(",");
            }
            tot_weight += this->mono_weight[path.at(idx)]; // Only a Sanity Check!!!
        }
        // Edge Case, there might by a qualifier to the end node
        cur_edge_id = get_edge_index(path.at(path.size()-2), path.at(path.size()-1));
        if (std::string(&this->qualifiers_str[this->qualifiers_str_index[cur_edge_id]]).length() > 0) {
                qualifiers.append(std::string(&this->qualifiers_str[this->qualifiers_str_index[cur_edge_id]]));
                qualifiers.append(",");
            }

        // Get epos
        for (uint32_t idx = path.size()-2; idx > 0; idx--) {
            sequence_intermediate = std::string(&this->sequence_str[this->sequence_str_index[path.at(idx)]]);
            if (0 != sequence_intermediate.compare("")) {
                // Retrieve ending pos
                if (this->iso_position[path.at(idx)] != UINT16_MAX){
                    epos = std::to_string(this->iso_position[path.at(idx)] + std::string(&this->sequence_str[this->sequence_str_index[path.at(idx)]]).length() -1);
                } else if (this->position[path.at(idx)] != UINT16_MAX) {
                    epos = std::to_string(this->position[path.at(idx)] + std::string(&this->sequence_str[this->sequence_str_index[path.at(idx)]]).length() -1);
                } else {
                    epos = '?';
                }
                break;
            }
        }

        // Make sequence FASTA-conform by adding "\n" every 60 characters
        for (int i = 0; i < sequence.size(); i+=60) {
            if (i+60 > sequence.size()) {
                sequence_new_lined.append(sequence.substr(i, sequence.size()-i) + "\n");
            } else {
                sequence_new_lined.append(sequence.substr(i, 60) + "\n");
            }
        }

        fasta += ">pg|TODO|" + this->accessions[iso_idx] + "(" + spos + ":" + epos
        + ",mssclvg:" +  std::to_string(mssclvg) + "," + qualifiers.substr(0, qualifiers.length()-1) + ")\n"
        + sequence_new_lined; //+ "W=" + std::to_string(tot_weight) + "\n";


    }
    return fasta;
}


bool ProteinGraph::overlapping_interval(uint32_t node_num, double lower, double upper) {
    uint32_t pdb_n = 2*this->PDB;
    double pdb_lower, pdb_upper;

    // Iterate over all pdb entries 
    for (uint32_t i = 0; i < pdb_n; i+=2) {
        pdb_lower = pdbs[node_num*pdb_n  + i];
        // As long as an interval is specified
        if ((double)INT64_MAX == pdb_lower) { break; }  // No more intervals to look into!

        // Check if it is overlapping!
        pdb_upper = this->pdbs[node_num*pdb_n  + i + 1];
        if (
            pdb_lower <= (double)upper
            &&
            pdb_upper >= (double)lower
        ) {
            return true;
        }
    }

    return false;
};

/*------------------------------------------------------------------------------------------------------*/
/*----------------------------------Float Implementation----------------------------------------------*/
/*------------------------------------------------------------------------------------------------------*/
std::string ProteinGraph::tvs_traverse_naive(int64_t lower, int64_t upper) {
    
    // State information
    std::unordered_map<uint32_t, std::vector<double>> tv_vals;  // tv vals achieved and currently achieved by expanding
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> paths;  // Paths which was taken to achieve the corresponding tv_val
    double f_lower = (double)lower, f_upper = (double)upper;  // Convert query to doubles

    // Variables during traversal
    uint32_t e_b, e_e, target_node;
    double new_lower, new_upper, achieved;

    // Initial values for traversal
    tv_vals[0] = {0};
    paths[0] = {{0}};

    // For every node (in top order)
    for (uint32_t i = 0; i < this->N-1; i++) {
        // skip if tv_vals is already empty (no paths!)
        if (tv_vals[i].size() == 0) {continue;}

        // Get beginning and ending of edge-ids
        if (i == 0) {
            e_b = 0; e_e = this->nodes[i];
        } else {
            e_b = this->nodes[i - 1]; e_e = this->nodes[i];
        }
        
        // For every possible path
        for (uint32_t j = 0; j < tv_vals[i].size(); j++) {   
            // For every outgoing edge of the node
            for (uint32_t k = e_b; k < e_e; k++) {
                // Calculated the achieved weight, lower, upper and target_node
                achieved = (tv_vals[i][j] + this->mono_weight[this->edges[k]]);
                new_lower = f_lower - achieved;  // New lower 
                new_upper = f_upper - achieved;  // New upper
                target_node = this->edges[k];  // Target of Edge

                // Check if we expand on this node
                if (this->overlapping_interval(
                        target_node, 
                        new_lower, new_upper
                        )
                    ) {
                    // CASE: Expanding
                    // Add new tv_val
                    tv_vals[target_node].push_back(achieved);

                    // Add new path how we achieved it
                    paths[target_node].push_back(paths[i][j]);
                    paths[target_node].back().push_back(target_node);
                } 
                // CASE: No Exanding --> Skip entry
            }
        }
        
        // Free memory during traversal, since older results can be removed (-> dag)!
        tv_vals.erase(i);
        paths.erase(i);
    };

    // Return results
    return convert_paths_to_fasta(paths[this->N-1]);
};


std::string ProteinGraph::tvs_traverse_varcount_naive(int64_t lower, int64_t upper, uint8_t max_vars) {
    
    // State information
    std::unordered_map<uint32_t, std::vector<double>> tv_vals;  // tv vals achieved and currently achieved by expanding
    std::unordered_map<uint32_t, std::vector<uint8_t>> var_count;  // tv vals achieved and currently achieved by expanding
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> paths; // Paths which was taken to achieve the corresponding tv_val

    double f_lower = (double)lower, f_upper = (double)upper;  // Convert query to doubles

    // Variables during traversal
    uint32_t e_b, e_e, target_node;
    uint16_t current_var_count;
    double new_lower, new_upper, achieved;

    // Initial values for traversal
    tv_vals[0] = {0};
    var_count[0] = {0};
    paths[0] = {{0}};

    // For every node (in top order)
    for (uint32_t i = 0; i < this->N-1; i++) {
        // skip if tv_vals is already empty (no paths!)
        if (tv_vals[i].size() == 0) {continue;}

        // Get beginning and ending of edge-ids
        if (i == 0) {
            e_b = 0; e_e = this->nodes[i];
        } else {
            e_b = this->nodes[i - 1]; e_e = this->nodes[i];
        }
        
        // For every possible path
        for (uint32_t j = 0; j < tv_vals[i].size(); j++) {   
            // For every outgoing edge of the node
            for (uint32_t k = e_b; k < e_e; k++) {
                // Calculated the achieved weight, lower, upper and target_node
                achieved = (tv_vals[i][j] + this->mono_weight[this->edges[k]]);
                new_lower = f_lower - achieved;  // New lower 
                new_upper = f_upper - achieved;  // New upper
                target_node = this->edges[k];  // Target of Edge

                // Additionally count the variants
                current_var_count = var_count[i][j] + this->variant_count[k];

                // Check if we expand on this node
                if (
                    (current_var_count <= max_vars)
                    &&
                    this->overlapping_interval(
                        target_node, 
                        new_lower, new_upper
                        )
                    ) {
                    // CASE: Expanding
                    // Add new tv_val
                    tv_vals[target_node].push_back(achieved);

                    // Count up used variants
                    var_count[target_node].push_back(current_var_count);

                    // Add new path how we achieved it
                    paths[target_node].push_back(paths[i][j]);
                    paths[target_node].back().push_back(target_node);
                } 
                // CASE: No Exanding --> Skip entry
            }
        }
        
        // Free memory during traversal, since older results can be removed (-> dag)!
        tv_vals.erase(i);
        paths.erase(i);
    };

    // Return results
    return convert_paths_to_fasta(paths[this->N-1]);
};
