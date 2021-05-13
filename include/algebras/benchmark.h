#pragma once

#include <chrono>
#include <iostream>

class Timer {
public:
    Timer() : elapsed_recent(0), bPrinted(false) { saved_time = std::chrono::system_clock::now(); }
    ~Timer() {
        if (!bPrinted) {
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed = end - saved_time;
            std::cout << "\033[0;32m" << "Elapsed time: " << elapsed.count() << "s\033[0m\n";
        }
    }
    void print(const char* msg = "") {
        bPrinted = true;
        auto end = std::chrono::system_clock::now();
        elapsed_recent = end - saved_time;
        std::cout << "\033[0;32m" << msg << "Elapsed time: " << elapsed_recent.count() << "s\033[0m\n";
        saved_time = std::chrono::system_clock::now();
    }
    double Elapsed() const {
        std::chrono::duration<double> result = std::chrono::system_clock::now() - saved_time;
        return result.count();
    }
    void start() { saved_time = std::chrono::system_clock::now(); }
    std::chrono::duration<double> get_elapsed_recent() const { return elapsed_recent; }
    void SuppressPrint() { bPrinted = true; }
private:
    std::chrono::time_point<std::chrono::system_clock> saved_time;
    std::chrono::duration<double> elapsed_recent;
    bool bPrinted;
};