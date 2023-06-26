#ifndef GRAPHLOADER_H
#define GRAPHLOADER_H

#include <string>

#include "protein_graph.hpp"


class GraphLoader{
    public:
        virtual ~GraphLoader() = default;
        virtual std::vector<ProteinGraph>* loadGraphs(std::string fileLoc) = 0;
};


class GraphLoaderBinary: public GraphLoader {
    public:
        ~GraphLoaderBinary() = default;
        
        std::vector<ProteinGraph>* loadGraphs(std::string fileLoc);
};

#endif