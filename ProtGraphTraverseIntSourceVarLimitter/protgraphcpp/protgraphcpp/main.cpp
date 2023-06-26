#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mutex>
#include <ostream>
#include <queue>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <vector>
#include <array>
#include <memory>
#include <thread>
#include <tuple>
#include <chrono>
#include <unordered_map>
#include <cmath>



#include "protein_graph.hpp"
#include "graph_loader.hpp"


#define QUEUE_SIZE 10000

template<class T, size_t MaxQueueSize>
class Queue
{
    std::condition_variable consumer_, producer_;
    std::mutex mutex_;
    using unique_lock = std::unique_lock<std::mutex>;

    std::queue<T> queue_;

public:
    template<class U>
    void push_back(U&& item) {
        unique_lock lock(mutex_);
        while(MaxQueueSize == queue_.size())
            producer_.wait(lock);
        queue_.push(std::forward<U>(item));
        consumer_.notify_one();
    }

    T pop_front() {
        unique_lock lock(mutex_);
        while(queue_.empty())
            consumer_.wait(lock);
        auto full = MaxQueueSize == queue_.size();
        auto item = queue_.front();
        queue_.pop();
        if(full)
            producer_.notify_all();
        return item;
    }
};


// TODO output queue
void thread_lifecycle(
    Queue<std::tuple<int64_t, int64_t, uint8_t>, QUEUE_SIZE>& query, 
    Queue<std::string, QUEUE_SIZE>& output_queue,
    std::vector<ProteinGraph>& pgs,
    std::atomic<bool>* pgs_executed,
    std::atomic<uint32_t>& atomic_pgs_finished,
    int64_t max_query,
    uint32_t num_bins
    ){

        bool iter_bool = false, should_execute = iter_bool;

        std::tuple<int64_t, int64_t, uint8_t> next_query;
        std::string output_benchmark = "";

        int64_t lower, upper;
        uint32_t used_bin;
        uint8_t num_vars_max;
        while(true) {
            // Get next query
            next_query = query.pop_front();

            // Stop Condition reached terminate thread
            if (std::get<2>(next_query) == uint8_t(-1)) break;

            //Set query params
            lower = std::get<0>(next_query);
            upper = std::get<1>(next_query);

            // Get the bin to use
            used_bin = (uint32_t) std::ceil( (upper / (max_query / num_bins))) - 1;
            if (used_bin >= num_bins) {
                used_bin = num_bins - 1;
            }


            // Scan bool vector
            for (uint32_t i = 0; i < pgs.size(); i++){
                // Check if it was already executed by another thread
                pgs_executed[i].compare_exchange_strong(should_execute, !iter_bool);
                if (should_execute == iter_bool) {
                    // It was not executed by another thread, execute now!
                    num_vars_max = pgs.at(i).max_vars_bins[used_bin];
                    if (num_vars_max == 255) {
                        output_benchmark += pgs.at(i).tvs_traverse_naive(lower,  upper);
                    } else {
                        output_benchmark += pgs.at(i).tvs_traverse_varcount_naive(lower,  upper, num_vars_max);
                    }


                    if (output_benchmark.length() != 0) {
                        output_queue.push_back(output_benchmark);
                        output_benchmark = "";
                    }
                    int finished_count = atomic_pgs_finished.fetch_add(1, std::memory_order_acq_rel);
                    if (finished_count >= pgs.size() - 1) {
                        //  Only 1 Thread should be active here!!!
                        // std::cerr << "Thread pushing end signal" << std::endl;
                        output_queue.push_back("TODO FINISHED CALCULATING");
                        atomic_pgs_finished.exchange(0);
                    }
                } 
                // Update should_execute due to the atomic replacing its value
                should_execute = iter_bool;
            }
            // Update bools
            iter_bool = !iter_bool;
            should_execute = iter_bool;
            
            
        }

}

int main(int argc, char *argv[]) {
    // printf("ProtGraphCpp, a (hopefully fast) Peptide Query Engine\n");


    // Retrieve parameters from argparse (database )
    std::string FILENAME = argv[1];
    std::string::size_type idx = FILENAME.rfind('.');
    std::string extension = FILENAME.substr(idx+1);

    // Limits of each protein
    std::string var_limit_csv = argv[5];
    std::ifstream var_limits_if(var_limit_csv);
    std::string var_line;
    std::stringstream var_ss_line;
    std::string var_entry;
    std::string protein_entry;

    uint32_t num_bins;
    std::vector<uint64_t> bins;
    std::unordered_map<std::string, std::vector<uint8_t>> max_vars;

    std::cout << "Loading Variants" << std::endl;
    while (std::getline(var_limits_if, var_line)) {
        // Parse Query
        var_ss_line.clear();
        var_ss_line.str(var_line);
        std::getline(var_ss_line, var_entry, ',');

        if (var_entry.compare("#bins") == 0) {
            std::getline(var_ss_line, var_entry, '\n');
            num_bins = (uint32_t) std::stod(var_entry);
        } else if (var_entry.compare("bins") == 0) {
            for (int i=0; i<num_bins-1; i++) {
                std::getline(var_ss_line, var_entry, ',');
                bins.push_back((int64_t)(std::stod(var_entry) * 1000000000));
            }
            std::getline(var_ss_line, var_entry, '\n');
            bins.push_back((int64_t)(std::stod(var_entry) * 1000000000));
        } else /* It is a protein */  {
            protein_entry = var_entry;
            for (int i=0; i<num_bins-1; i++) {
                std::getline(var_ss_line, var_entry, ',');
                max_vars[protein_entry].push_back((uint8_t)(std::stoi(var_entry)));
            }
            std::getline(var_ss_line, var_entry, '\n');
            max_vars[protein_entry].push_back((uint8_t)(std::stoi(var_entry)));
        }
    }

    std::cout << "Loading Graphs" << std::endl;
    GraphLoader* gl;
    if (extension.compare("bpcsr") == 0) {
        gl = new GraphLoaderBinary();
    }

    // Get Protein Graphs
    std::vector<ProteinGraph>* pgs = gl->loadGraphs(FILENAME, max_vars);
    printf("\n");

    // Get the number of available threads
    int num_threads = 0;
    if (atoi(argv[3]) == -1) {
        num_threads = std::thread::hardware_concurrency();
    } else {
        num_threads = atoi(argv[3]);
    }    

    std::ofstream output_file(argv[4]);

    

    // Set Queues
    Queue<std::tuple<int64_t, int64_t, uint8_t>, QUEUE_SIZE> query;
    Queue<std::string, QUEUE_SIZE> output_queue;

    // Set bool vector if executed
    // std::vector<bool> pgs_executed(pgs->size(), false);
    // bool* pgs_executed = new bool[pgs->size()];
    // for (int i = 0; i < pgs->size(); i++){
    //     pgs_executed[i] = false;
    // }
    std::atomic<bool>* pgs_executed = new std::atomic<bool>[pgs->size()];
    for (int i = 0; i < pgs->size(); i++){
        pgs_executed[i].exchange(false);
    }

    // Create atomic counter of executed entries
    std::atomic<uint32_t> atomic_pgs_finished{0};


    // Create Thread pool after retrieving all needed information
    std::cout << "Starting Threads" << std::endl;
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++) {
        threads.push_back(std::thread(
            thread_lifecycle, // Method
            std::ref(query), std::ref(output_queue), //Params
            std::ref(*pgs), std::ref(pgs_executed),
            std::ref(atomic_pgs_finished),
            bins.back(),
            num_bins
            
        ));
    }


    // Now serve queries
    // TODO Currently we simply read a csv
    std::string QUERY_FILE = argv[2];
    std::ifstream query_file(QUERY_FILE);
    std::string line;
    std::stringstream ss_line;
    std::string entry;
    int64_t lower;
    int64_t upper;
    int query_counter = 1;

    while (std::getline(query_file, line)) {
        // Parse Query
        ss_line.clear();
        ss_line.str(line);
        std::getline(ss_line, entry, ',');
        lower = (int64_t)(std::stod(entry) * 1000000000);
        std::getline(ss_line, entry, '\n');
        upper = (int64_t)(std::stod(entry) * 1000000000);
        std::tuple<int64_t, int64_t, uint8_t> query_tuple(lower, upper, 1);

        //Submit Query
        for (int i = 0; i < num_threads; i++) {
            query.push_back(query_tuple);
        }

        // Wait for the results!
        while (true) {
            std::string output = output_queue.pop_front();
            if (output.starts_with("TODO FINISHED CALCULATING")){
                break;
            }
            // std::cout << output; // E.G.: here we could simply pass it through the socket
            output_file << output;

        }

        std::cerr << "Processed Query " << query_counter << " with: " << lower << ":" << upper << std::endl;
        query_counter++;
        // 18 446 744 073.709553
        //  9 223 372 036.854776
    }
    // Repeat next Query TODO 



    // Spin down
    // std::cout << "Spinning Down..." << std::endl;
    std::tuple<int64_t, int64_t, uint8_t> stop_tuple(0, 0, uint8_t(-1));
    for (int i = 0; i < num_threads*2; i++)
    {
        query.push_back(stop_tuple);
    }
    for (int i = 0; i < num_threads; i++)
    {
        threads.at(i).join();
    }

    // printf("Completely finished!\n");
    return  0;
}