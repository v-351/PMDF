#pragma once

#include <chrono>

/**
 * @brief Тайминг.
 * Регистрация временных промежутков.
 */
class Timing{
private:
    std::chrono::steady_clock::time_point start_point;
    std::chrono::seconds duration;
    bool is_ticking = false;
public:

    Timing(){
        duration.zero();
        is_ticking = false;
    }

    void start(){
        if(!is_ticking){
            start_point = std::chrono::steady_clock::now();
            is_ticking = true;
        }
    }
    void pause(){
        duration += std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_point);
        is_ticking = false;
    }

    void reset(){
        duration.zero();
        is_ticking = false;
    }
    
    double result(){
        if(is_ticking){
            pause();
        }
        return duration.count();
    }
};