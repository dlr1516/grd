/**
 * LIRE - LIdar 3D REgistration
 * Copyright (C) 2018 Dario Lodi Rizzini, Fabio Riccardo Turdo
 *
 * LIRE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * LIRE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LIRE.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <grd/Profiler.h>


namespace grd {

    // ----------------------------------------------------
    // PROFILER
    // ----------------------------------------------------

    void Profiler::updateStat(std::string label, double time) {
        auto it = stats_.find(label);
        if (it == stats_.end()) {
            MeasureStatistic ms(time);
            stats_.insert(std::make_pair(label, ms));
        } else {
            // Computes mean value and variance with online Welford’s method 
            double timeAvgOld = it->second.timeAvg;
            it->second.timeAvg += (time - timeAvgOld) / (double) it->second.count;
            it->second.timeVar += (time - it->second.timeAvg) * (time - timeAvgOld);
            if (time < it->second.timeMin) {
                it->second.timeMin = time;
            }
            if (time > it->second.timeMax) {
                it->second.timeMax = time;
            }
            it->second.count++;
        }
    }

    void Profiler::printStats(std::ostream& out) const {
        for (auto s : stats_) {
            out << s.first << " \t"
                    << s.second.timeAvg << " \t"
                    << s.second.getVariance() << " \t"
                    << s.second.timeMin << " \t"
                    << s.second.timeMax << " \t"
                    << s.second.count << " \t"
                    << std::endl;
        }
    }

    // ----------------------------------------------------
    // SCOPED TIMER
    // ----------------------------------------------------

    ScopedTimer::ScopedTimer(std::string label) : label_(label), timeStart_(timer_type::now()) {
        // before: timeStart_(std::chrono::high_resolution_clock::now())
    }

    ScopedTimer::~ScopedTimer() {
        Profiler::getProfiler().updateStat(label_, elapsedTimeMs());
    }

    double ScopedTimer::elapsedTimeMs() const {
        //auto elapsedNanosec = std::chrono::duration_cast<double,std::chrono::milliseconds>(std::chrono::steady_clock::now() - timeStart_).count();
        //double timeElapsed = 1.0 * (double)elapsedNanosec;
        std::chrono::duration<double, std::milli> timeElapsed = timer_type::now() - timeStart_;
        return timeElapsed.count();
    }

} // end of namespace 