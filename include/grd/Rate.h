#pragma once

#include <chrono>
#include <thread>

namespace grd {

    /** Inspired by ros::Rate.
     */
    class Rate {
    public:
        typedef std::chrono::duration<int, std::milli> milliseconds_type;

        /** Constructor.
         */
        Rate(double frequency)
        : period_((int) (1000.0 / frequency)), start_(std::chrono::system_clock::now()) {
        }

        /**
         */
        void sleep() {
            // Computes the remaining time in the cycle
            //    std::chrono::system_clock::time_point stop = start_ + period_; 
            //    auto remaining = stop - chrono::system_clock::now();
            //    auto remaining_ms = chrono::duration_cast<std::chrono::milliseconds>(diff);

            // Sleep
            std::this_thread::sleep_until(start_ + period_);
            start_ = std::chrono::system_clock::now();
        }

    private:
        std::chrono::system_clock::time_point start_;
        milliseconds_type period_;
    };

} // end of namespace 
