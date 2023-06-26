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


void throw_event_function(
    std::atomic<bool>* event,
    std::atomic<bool>* timout,
    int seconds
) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(); // Start measuring time
    // poll time until 
    
    while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - begin).count() < seconds) {
        // Constantly poll to exit directly (happens on most proteins, which only need a few microseconds)
        // For other proteins we simply poll until the time is passed
        if (timout->load() == true) {
            return;
        }
    }
    event->store(true);
}

// TODO output queue
void thread_lifecycle(
    Queue<std::tuple<int64_t, int64_t, uint8_t>, QUEUE_SIZE>& query, 
    Queue<std::string, QUEUE_SIZE>& output_queue,
    std::vector<ProteinGraph>& pgs,
    std::atomic<bool>* pgs_executed,
    std::atomic<bool>* pgs_timeout,
    std::atomic<uint32_t>& atomic_pgs_finished,
    uint32_t limit_query_in_seconds
    ){

        bool iter_bool = false, should_execute = iter_bool;

        std::tuple<int64_t, int64_t, uint8_t> next_query;
        std::string output_benchmark = "";

        int64_t lower, upper;
        while(true) {
            // Get next query
            next_query = query.pop_front();

            // Stop Condition reached terminate thread
            if (std::get<2>(next_query) == uint8_t(-1)) break;

            //Set query params
            lower = std::get<0>(next_query);
            upper = std::get<1>(next_query);

            // Scan bool vector
            for (uint32_t i = 0; i < pgs.size(); i++){
                // Check if it was already executed by another thread
                pgs_executed[i].compare_exchange_strong(should_execute, !iter_bool);
                if (should_execute == iter_bool) {

                    std::atomic<bool> timeout_timeout;
                    timeout_timeout.store(false);
                    std::thread throw_event = std::thread(
                        throw_event_function,
                        &pgs_timeout[i],
                        &timeout_timeout,
                        limit_query_in_seconds
                    );

                    // It was not executed by another thread, execute now!
                    output_benchmark += pgs.at(i).tvs_traverse_naive(lower,  upper, &pgs_timeout[i]);

                    timeout_timeout.store(true);
                    throw_event.join();
                    pgs_timeout[i].store(false);


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

// TODO output queue
void thread_lifecycle_var_count(
    Queue<std::tuple<int64_t, int64_t, uint8_t>, QUEUE_SIZE>& query, 
    Queue<std::string, QUEUE_SIZE>& output_queue,
    std::vector<ProteinGraph>& pgs,
    std::atomic<bool>* pgs_executed,
    std::atomic<bool>* pgs_timeout,
    std::atomic<uint32_t>& atomic_pgs_finished,
    uint8_t varcount,
    uint32_t limit_query_in_seconds
    ){

        bool iter_bool = false, should_execute = iter_bool;

        std::tuple<int64_t, int64_t, uint8_t> next_query;
        std::string output_benchmark = "";

        int64_t lower, upper;
        while(true) {
            // Get next query
            next_query = query.pop_front();

            // Stop Condition reached terminate thread
            if (std::get<2>(next_query) == uint8_t(-1)) break;

            //Set query params
            lower = std::get<0>(next_query);
            upper = std::get<1>(next_query);

            // Scan bool vector
            for (uint32_t i = 0; i < pgs.size(); i++){
                // Check if it was already executed by another thread
                pgs_executed[i].compare_exchange_strong(should_execute, !iter_bool);
                if (should_execute == iter_bool) {

                    std::atomic<bool> timeout_timeout;
                    timeout_timeout.store(false);
                    std::thread throw_event = std::thread(
                        throw_event_function,
                        &pgs_timeout[i],
                        &timeout_timeout,
                        limit_query_in_seconds
                    );

                    // It was not executed by another thread, execute now!
                    output_benchmark += pgs.at(i).tvs_traverse_varcount_naive(lower,  upper, varcount, &pgs_timeout[i]);

                    timeout_timeout.store(true);
                    throw_event.join();
                    pgs_timeout[i].store(false);


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

    GraphLoader* gl;
    if (extension.compare("bpcsr") == 0) {
        gl = new GraphLoaderBinary();
    }

    // Get Protein Graphs
    std::vector<ProteinGraph>* pgs = gl->loadGraphs(FILENAME);
    printf("\n");

    // Get the number of available threads
    int num_threads = 0;
    if (atoi(argv[3]) == -1) {
        num_threads = std::thread::hardware_concurrency();
    } else {
        num_threads = atoi(argv[3]);
    }    

    // Output CSV File
    std::ofstream output_file(argv[4]);

    // Check if we want to limit by number of variants
    int var_limit = atoi(argv[5]);

    // Time in seconds when to stop a search and return -1 (not in time)
    uint32_t limit_query_in_seconds = atoi(argv[6]);

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
    std::atomic<bool>* pgs_timeout  = new std::atomic<bool>[pgs->size()];

    for (int i = 0; i < pgs->size(); i++){
        pgs_executed[i].exchange(false);
        pgs_timeout[i].exchange(false);
    }

    // Create atomic counter of executed entries
    std::atomic<uint32_t> atomic_pgs_finished{0};


    // Create Thread pool after retrieving all needed information
    // std::cout << "Spinning Up..." << std::endl;
    std::vector<std::thread> threads;
    if (var_limit == -1) {
        for (int i = 0; i < num_threads; i++) {
            threads.push_back(std::thread(
                thread_lifecycle, // Method
                std::ref(query), std::ref(output_queue), //Params
                std::ref(*pgs), std::ref(pgs_executed), std::ref(pgs_timeout),
                std::ref(atomic_pgs_finished),
                limit_query_in_seconds
            ));
        }
    } else {
        for (int i = 0; i < num_threads; i++) {
            threads.push_back(std::thread(
                thread_lifecycle_var_count, //Method
                std::ref(query), std::ref(output_queue), //Params
                std::ref(*pgs), std::ref(pgs_executed), std::ref(pgs_timeout),
                std::ref(atomic_pgs_finished),
                var_limit,
                limit_query_in_seconds
            ));
        }
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
    }
    // Repeat next Query TODO 



    // Spin down
    // std::cout << "Spinning Down..." << std::endl;
    std::tuple<int64_t, int64_t, uint8_t> stop_tuple(0, 0, uint8_t(-1));
    for (int i = 0; i < num_threads; i++)
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