/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017 Dario Lodi Rizzini.
 *
 * ARS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ARS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ARS.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <grd/BBOptimizer1d.h>


namespace grd {

    // --------------------------------------------------------
    // OPTIMIZATION
    // --------------------------------------------------------

    BBOptimizer1d::BBOptimizer1d() : xtol_(0.001), ytol_(0.001), xtollOn_(true), ytollOn_(false) {
    }

    BBOptimizer1d::BBOptimizer1d(double xtol, double ytol) : xtol_(xtol), ytol_(ytol), xtollOn_(true), ytollOn_(false) {
    }

    BBOptimizer1d::~BBOptimizer1d() {
    }

    void BBOptimizer1d::findGlobalMax(double xmin, double xmax, double& x, double& ylower, double& yupper) {
        LeastUpperBoundFirstQueue queue;
        IntervalBound curr, left, right, global;

        if (xmax < xmin) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid interval endpoints: xmin " << xmin << " must be less than xmax " << xmax
                    << "\n  swapping values to continue" << std::endl;
            std::swap(xmin, xmax);
        }

        // At least one tolerance must be enabled to guarantee a stop condition! 
        // If not, default choice (xtol_ on) is chosen. 
        if (!xtol_ && !ytol_) {
            xtol_ = true;
        }

        // Initializes the queue with the minimum and maximum
        global.xmin = xmin;
        global.xmax = xmax;
        findLU(global.xmin, global.xmax, global.ylower, global.yupper);
        queue.push(global);
        // Extracts candidate interval from queue and split them to refine lower/upper value estimation
        while (!queue.empty()) {
            curr = queue.top();
            queue.pop();
            // It processes the interval curr only if the further analysis can improve the current solution. 
            // In practice, if curr.upper < global.lower, then curr is ignored.
            if (curr.yupper >= global.ylower) {
                // Updates: the current candidate interval to contain global maximum with the lower and upper bounds
                // of global maximum
                if (global.ylower <= curr.ylower) {
                    global.xmin = curr.xmin;
                    global.xmax = curr.xmax;
                    global.ylower = curr.ylower;
                    global.yupper = curr.yupper;
                }
                //      std::cout << "curr " << STREAM_BOUND_INTERVAL(curr) << std::endl;
                //      std::cout << "  global " << STREAM_BOUND_INTERVAL(global) << std::endl; 
                // Splits curr into intervals left and right, if the interval does not satisfy stop criteria
                if (!checkInterval(curr)) {
                    left.xmin = curr.xmin;
                    left.xmax = 0.5 * (curr.xmin + curr.xmax);
                    findLU(left.xmin, left.xmax, left.ylower, left.yupper);
                    right.xmin = left.xmax;
                    right.xmax = curr.xmax;
                    findLU(right.xmin, right.xmax, right.ylower, right.yupper);
                    queue.push(left);
                    queue.push(right);
                }
            }
            //    else {
            //      std::cout << "  discarding " << STREAM_BOUND_INTERVAL(curr) << std::endl;
            //    }
        }
        x = 0.5 * (global.xmin + global.xmax);
        ylower = global.ylower;
        yupper = global.yupper;
    }

    bool BBOptimizer1d::checkInterval(const IntervalBound& ib) const {
        if (xtollOn_ && (ib.xmax - ib.xmin > xtol_)) {
            return false;
        }
        if (ytollOn_ && (ib.yupper - ib.ylower > ytol_)) {
            return false;
        }
        return true;
    }

    // --------------------------------------------------------

    FourierOptimizerBB1D::FourierOptimizerBB1D()
    : BBOptimizer1d(), coeffs_(), orderMax_(0) {
    }

    FourierOptimizerBB1D::FourierOptimizerBB1D(const std::vector<double>& coeffs)
    : BBOptimizer1d(), coeffs_(), orderMax_(0) {
        setCoefficients(coeffs);
    }

    FourierOptimizerBB1D::~FourierOptimizerBB1D() {
    }

    void FourierOptimizerBB1D::findLU(double xmin, double xmax, double& ylower, double& yupper) {
        findLUFourier(coeffs_, xmin, xmax, ylower, yupper);
    }

    // --------------------------------------------------------

    void findGlobalMaxBBFourier(const std::vector<double>& coeffs, double theta0, double theta1, double thetaTol, double fourierTol, double& thetaMax, double& fourierMax) {
        FourierOptimizerBB1D fopt(coeffs);
        double fourierLower, fourierUpper;
        if (thetaTol > 0.0) {
            fopt.setXTolerance(thetaTol);
            fopt.enableXTolerance(true);
        }
        if (fourierTol > 0.0) {
            fopt.setYTolerance(fourierTol);
            fopt.enableYTolerance(true);
        }
        fopt.findGlobalMax(theta0, theta1, thetaMax, fourierLower, fourierUpper);
        fourierMax = 0.5 * (fourierLower + fourierUpper);
    }

} // end of namespace
