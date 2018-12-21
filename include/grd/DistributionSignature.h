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
#ifndef GRD_DISTRIBUTION_SIGNATURE_H
#define GRD_DISTRIBUTION_SIGNATURE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include <grd/Point.h>

namespace grd {

    class DistributionSignature {
    public:

        enum RangeDistribution {
            BIASED_RAYLEIGH, ERLANG
        };

        /**
         * Default constructor.
         */
        DistributionSignature();

        /**
         * Default constructor.
         */
        DistributionSignature(const DistributionSignature& ds);

        /** Constructor with the size of Glare Signature.
         */
        DistributionSignature(int thetaN, int rangeN, RangeDistribution distr = BIASED_RAYLEIGH);

        /**
         * Destructor.
         */
        virtual ~DistributionSignature();

        /** 
         * Resets the content of distribution signature.
         */
        void clear();

        /**
         * Sets the dimension of GLARE histogram. 
         * @param thetaNum
         * @param rhoNum
         * @param rhoStep
         */
        void setSize(int thetaN, int rangeN);

        /**
         * Sets the type of range distribution. 
         * @param distr
         */
        void setRangeDistribution(RangeDistribution distr);

        /**
         * Sets the "concentration" parameter kappa of von Mises distribution. 
         * @param kappa
         */
        void setVonMisesKappa(double kappa) {
            vonMisesKappa_ = kappa;
            init();
        }

        /**
         * Sets the "width" of Bias Rayleigh mode. 
         * @param sigma
         */
        void setBiasedRayleighSigma(double sigma) {
            biasraySigma_ = sigma;
        }

        /**
         * Sets the tolerance on the estimation of the
         */
        void setCorrelationThetaTolerance(double thtol) {
            correlThetaTol_ = thtol;
        }

        /**
         * Sets the tolerance on the estimation of the
         */
        void setCorrelationFourierTolerance(double foutol) {
            correlFourierTol_ = foutol;
        }

        /** Inserts features: each feature has a field point. 
         */
        void setPoints(const VectorPoint2d& points);

        /** Returns the const reference to points.
         */
        const VectorPoint2d& getPoints() const {
            return points_;
        }

        /**
         * Returns a const reference to the matrix of coefficients for the 
         * Fourier-Laguerre functions. 
         * @return a const reference to the matrix of coefficients
         */
        const Eigen::MatrixXd& getSignatureCoeffs() const {
            return signature_;
        }

        /**
         * Returns the auto-correlation of signature. 
         * @return 
         */
        double getAutoCorrel() const {
            return autoCorr_;
        }

        /**
         * Evaluates the signature distribution value using the approximant series 
         * decomposition.
         * @param theta
         * @param r
         * @return the 
         */
        double evaluate(double theta, double r) const;

        /**
         * Returns the correlation value of the two signatures. The correlation is 
         * computed on the maximum angular shift that better aligns the two signatures. 
         * @param dsig the signature to be compared
         * @param correlCoeffs
         * @return 
         */
        double correlate(const DistributionSignature& dsig, double& thetaShift, std::vector<double>& correlCoeffs) const;

        /**
         * Returns the correlation value of the two signatures. The correlation is 
         * computed on the maximum angular shift that better aligns the two signatures. 
         * @param dsig the signature to be compared
         * @param thetaShift the best angular shift
         * @return the correlation value
         */
        double correlate(const DistributionSignature& dsig, double& thetaShift) const;

        /**
         * Returns the correlation value of the two signatures. The correlation is 
         * computed on the maximum angular shift that better aligns the two signatures. 
         * @param dsig the signature to be compared
         * @return the correlation value
         */
        double correlate(const DistributionSignature& dsig) const;

    private:
        Eigen::MatrixXd signature_;
        VectorPoint2d points_;
        int thetaN_; // order of Fourier transform
        int rangeN_; // order of Laguerre polynomial series
        RangeDistribution distribution_;
        double vonMisesKappa_;
        std::vector<double> vonMisesBesselRatios_;
        double biasraySigma_;
        double correlThetaTol_;
        double correlFourierTol_;
        double autoCorr_;

        void init();

        void updateCoeffs(double range, double angle);
    };

} // end of namespace 

#endif /* DISTRIBUTIONSIGNATURE_H */

