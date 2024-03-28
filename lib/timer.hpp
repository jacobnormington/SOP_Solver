#ifndef TIMER_H
#define TIMER_H
#include <unordered_map>
#include <chrono>

class timer
{
private:
    std::unordered_map<std::string, std::chrono::time_point<std::chrono::system_clock>> markers;
    std::chrono::time_point<std::chrono::system_clock> epoch;

public:
    timer();                                     // contructor - starts timer
    void add_new_marker(std::string marker);     // adds a marker to the timer
    void restart();                              // restarts the main timer
    double get_time_seconds();                   // returns the main time in seconds
    double get_time_seconds(std::string marker); // returns the specified marker time in seconds
    double get_time_millis();                    // returns the main time in milliseconds
    double get_time_millis(std::string marker);  // returns the specified marker time in milliseconds
};

#endif