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
#ifndef GRD_PROFILER_H
#define GRD_PROFILER_H

#include <iostream>
#include <map>
#include <chrono>

namespace grd {

    // Declarations 
    class Profiler;
    class ScopedTimer;

    /**
     * Class Profiler collects the statistics about the execution time. 
     * It is implemented as a Singleton. 
     */
    class Profiler {
    public:

        /**
         * Struct for storing time stats. 
         */
        struct MeasureStatistic {
            double timeAvg;
            double timeVar;
            double timeMin;
            double timeMax;
            size_t count;

            MeasureStatistic(double time) : timeAvg(time), timeVar(0.0), timeMin(time), timeMax(time), count(1) {
            }

            double getVariance() const {
                if (count > 1) {
                    return (timeAvg / (count - 1));
                } else {
                    return 0.0;
                }
            }
        };

        /**
         * Returns the only instance of profiler.
         */
        static inline Profiler& getProfiler() {
            static Profiler profiler;
            return profiler;
        }

        /**
         * Updates the statistics associated to a given label. 
         * @param label identifier of a series of measurements
         * @param time the last measurement of the event associated to a label
         */
        void updateStat(std::string label, double time);

        /**
         * Prints the statistics on the given output stream
         * @param out the output stream
         */
        void printStats(std::ostream& out) const;

    protected:
        std::map<std::string, MeasureStatistic> stats_;

        /**
         * Default constructor as private member of the class. 
         */
        Profiler() : stats_() {
        }

        /**
         * Default destructor. 
         */
        ~Profiler() {
        }

        /**
         * Copy constructor as private member of the class. 
         */
        Profiler(const Profiler& p) = delete;

        /**
         * Assignment operator as private member. 
         */
        void operator=(const Profiler& p) = delete;
    };

    /**
     * ScopedTimer measures the time elapsed from constructor to the invocation of destructor. 
     */
    class ScopedTimer {
    public:
        typedef std::chrono::steady_clock timer_type;
        //        typedef std::chrono::high_resolution_clock timer_type;

        /**
         * Constructor of timer inside scope with the given label
         * @param label label associated to a measurement serie
         */
        ScopedTimer(std::string label);

        /**
         * Destructor. It saves stats when destructor is called. 
         */
        ~ScopedTimer();

        /**
         * Returns the elapsed time in milliseconds. 
         */
        double elapsedTimeMs() const;

    protected:
        std::string label_;
        //std::chrono::time_point<std::chrono::high_resolution_clock> timeStart_;
        std::chrono::time_point<timer_type> timeStart_;
    };

} // end of namespace 

#endif /* PROFILER_H */