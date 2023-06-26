#ifndef PROTEINGRAPH_H
#define PROTEINGRAPH_H
#include <fstream>

#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <atomic>


class ProteinGraph {
    public:
        ProteinGraph(std::uint32_t num_acc, std::ifstream &input);
        ~ProteinGraph() = default;

        uint32_t N;
        uint32_t E;
        uint32_t PDB;
        std::vector<std::string> accessions;
        std::uint32_t* nodes; // Node
        std::uint8_t* iso_index;
        std::int64_t* mono_weight;

        std::uint32_t* edges; // Edge
        std::vector<bool> cleaved; // Edge Attrs  <-- as bool vector TODO DL


        // Other specific information
        std::uint8_t* variant_count;  // On Edges  <-- maybe parsing?
        std::int64_t* pdbs;  // On Nodes
        
        std::uint16_t* position;  // On Nodes TODO DL
        std::uint16_t* iso_position;  // On Nodes
        char* sequence_str; // Node/Edge Attrs (compacted, this could probably be also made seperately)
        std::uint32_t* sequence_str_index; // Node/Edge Attrs (compacted, this could probably be also made seperately)

        char* qualifiers_str; // Node/Edge Attrs (compacted, this could probably be also made seperately)
        std::uint32_t* qualifiers_str_index; // Node/Edge Attrs (compacted, this could probably be also made seperately)
        // char* qualifiers_str[]; 

        // TODO what methods to include?
        bool overlapping_interval(uint32_t node_num, int64_t lower, int64_t upper);
        std::string tvs_traverse_naive(int64_t lower, int64_t upper, std::atomic<bool>*);
        std::string tvs_traverse_varcount_naive(int64_t lower, int64_t upper, uint8_t max_vars, std::atomic<bool>*);

        std::uint32_t get_edge_index(uint32_t source_node, uint32_t target_node);
        std::string convert_paths_to_fasta(std::vector<std::vector<uint32_t>> paths);
};


#endif