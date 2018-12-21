/**
 * GRD - Geometric Relation Distribution
 * Copyright (C) 2018 Dario Lodi Rizzini.
 *
 * GRD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * GRD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with GRD.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef GRD_RATE_H
#define GRD_RATE_H

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

#endif
