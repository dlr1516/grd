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
#ifndef KERNEL_SIGNATURE_H
#define KERNEL_SIGNATURE_H

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

#include <grd/Point.h>

namespace grd {

    class KernelSignature {
    public:

        enum RangeDistribution {
            BIASED_RAYLEIGH, ERLANG
        };

        struct KernelParam {
            double theta;
            double range;
        };

        /**
         * Default constructor.
         */
        KernelSignature();

        /**
         * Destructor.
         */
        virtual ~KernelSignature();

        /**
         * Clears the kernels to set new points.
         */
        void clear() {
            kernels_.clear();
        }

        /**
         * Sets the type of range distribution. 
         * @param distr
         */
        void setRangeDistribution(RangeDistribution distr) {
            distribution_ = distr;
        }

        /**
         * Sets the "concentration" parameter kappa of von Mises distribution. 
         * @param kappa
         */
        void setVonMisesKappa(double kappa) {
            vonMisesKappa_ = kappa;
            vonMisesNorm_ = 1.0 / (2 * M_PI * boost::math::cyl_bessel_i(0, kappa));
        }

        /**
         * Sets the "width" of Bias Rayleigh mode. 
         * @param sigma
         */
        void setBiasedRayleighSigma(double sigma) {
            biasraySigma_ = sigma;
        }

        /** Inserts features: each feature has a field point. 
         */
        void setPoints(const VectorPoint2d& points);

        void displayKernels(std::ostream& out) const;

        /**
         * Evaluates the signature distribution value.
         * @param theta
         * @param r
         * @return the value of the signature distribution 
         */
        double evaluate(double theta, double r) const;

    private:
        std::vector<KernelParam> kernels_;
        RangeDistribution distribution_;
        double vonMisesKappa_;
        double vonMisesNorm_;
        double biasraySigma_;
    };

} // end of namespace 

#endif /* KERNELSIGNATURE_H */
