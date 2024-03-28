#include "timer.hpp"
#include <chrono>
timer::timer()
{
    epoch = std::chrono::system_clock::now();
}

void timer::restart()
{
    epoch = std::chrono::system_clock::now();
}

void timer::add_new_marker(std::string marker)
{
    markers[marker] = std::chrono::system_clock::now();
}

double timer::get_time_seconds()
{
    return (std::chrono::duration<double>(std::chrono::system_clock::now() - epoch)).count();
}

double timer::get_time_seconds(std::string marker)
{
    if (!markers.count(marker))
        return -1;
    return (std::chrono::duration<double>(std::chrono::system_clock::now() - markers[marker])).count();
}

double timer::get_time_millis()
{
    return (std::chrono::duration<double>(std::chrono::system_clock::now() - epoch)).count() * 1000;
}

double timer::get_time_millis(std::string marker)
{
    if (!markers.count(marker))
        return -1;
    return (std::chrono::duration<double>(std::chrono::system_clock::now() - markers[marker])).count() * 1000;
}
