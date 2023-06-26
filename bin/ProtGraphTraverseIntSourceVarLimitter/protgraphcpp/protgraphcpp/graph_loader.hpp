#ifndef GRAPHLOADER_H
#define GRAPHLOADER_H

#include <string>

#include "protein_graph.hpp"
#include <unordered_map>


class GraphLoader{
    public:
        virtual ~GraphLoader() = default;
        virtual std::vector<ProteinGraph>* loadGraphs(std::string fileLoc, std::unordered_map<std::string, std::vector<uint8_t>> max_vars) = 0;
};


class GraphLoaderBinary: public GraphLoader {
    public:
        ~GraphLoaderBinary() = default;
        
        std::vector<ProteinGraph>* loadGraphs(std::string fileLoc, std::unordered_map<std::string, std::vector<uint8_t>> max_vars);
};

#endif