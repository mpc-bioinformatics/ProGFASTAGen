#include <cstdint>
#include <cstring>
#include <string>
#include <fstream>
#include <endian.h>
#include <iostream>
#include "graph_loader.hpp"
#include "protein_graph.hpp"
#include <ulimit.h>
#include <unordered_map>
#include <vector>

        
std::vector<ProteinGraph>* GraphLoaderBinary::loadGraphs(std::string fileLoc, std::unordered_map<std::string, std::vector<uint8_t>> max_vars) {

    // malloc and then return reference?
    std::vector<ProteinGraph>* pgs = new std::vector<ProteinGraph>;  // Ptr to vector

    std::ifstream input(fileLoc, std::ios::binary);

    uint32_t num_acc, num_n, num_e, num_pdbs, cur_32bit;
    uint8_t cur_8bit;
    uint16_t cur_16bit;
    uint64_t cur_64bit;
    std::string cur_string;

    // Debug Counter 
    uint32_t counter = 0;

    // // Get next entry and its global information about the entry
    while (input.read((char*) &num_acc, 4)) {
        pgs->push_back(ProteinGraph(num_acc, input, max_vars));
        // Print Progress (TODO this might slow us down!!)
        // printf("\rParsed %i ProteinGraphs", counter);
        std::fflush(stdout);
        counter++;
        
    }

    // Closing the binary file, since eof is reached
    input.close();

    return pgs; // TODO DL return pointer
};
