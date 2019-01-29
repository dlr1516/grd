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
#ifndef GRD_GLAROT_SIGNATURE_H
#define GRD_GLAROT_SIGNATURE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include <grd/Point.h>

namespace grd {

    /** Computes GLARE and GLAROT.
     */
    class GlarotSignature {
    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        struct PointSignature {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            Eigen::Vector2d point;
            Eigen::MatrixXd signature;
        };

        //typedef std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > VectorPoint2d;
        typedef std::vector<PointSignature, Eigen::aligned_allocator<PointSignature> > VectorPointSignature;

        /** Constructor with default Glare signature size. 
         */
        GlarotSignature();

        /** Constructor with the size of Glare Signature.
         */
        GlarotSignature(int thetaNum, int rhoNum, double rhoStep);

        /** Resets the content of GLARE: 
         * -global signature becomes zero;
         * -point signatures are erased.
         */
        void clear();

        /**
         * Sets the dimension of GLARE histogram. 
         * @param thetaNum
         * @param rhoNum
         * @param rhoStep
         */
        void setSize(int thetaNum, int rhoNum, double rhoStep);

        /** Sets the covariance matrix representing uncertainty of vector [theta,rho].
         * 
         *  covar_ = [ \sigma_{\theta}^2  cov_{\theta,\rho} ]
         *           [ cov_{\theta,\rho}  \sigma_{\rho}^2   ]
         */
        void setCovariance(const Eigen::Matrix2d& cov) {
            covar_ = 0.5 * (cov + cov.transpose());
            computeGaussianKernel();
        }

        /** Sets the covariance matrix representing uncertainty of vector [theta,rho]
         * as a diagonal matrix using standard deviation values. 
         */
        void setCovariance(double sigmaTheta, double sigmaRho) {
            covar_ << sigmaTheta*sigmaTheta, 0.0,
                    0.0, sigmaRho*sigmaRho;
            computeGaussianKernel();
        }

        /** Inserts features: each feature has a field point. 
         */
        void setPoints(const VectorPoint2d& points);

        /** Returns the const reference to points.
         */
        const VectorPoint2d& getPoints() const {
            return points_;
        }

        /** Returns a reference to global signature.
         */
        const Eigen::MatrixXd& getSignatureGlobal() const {
            return signatureGlobal_;
        }

        /** Number of points used in the signature.
         */
        int getPointNum() const {
            return signaturePoints_.size();
        }

        /** Returns a reference to the point i-th.
         */
        const Eigen::Vector2d& getPoint(int i) const {
            return signaturePoints_[i].point;
        }

        /** Returns a reference to the signature of point i. 
         */
        const Eigen::MatrixXd& getSignaturePoint(int i) const {
            return signaturePoints_[i].signature;
        }

        /** Computes the L1 distance on global signature.
         */
        double distanceL1(const GlarotSignature& gs) const;

        /**
         */
        double distanceL1Min(const GlarotSignature& gs) const;

        /** Matches points according to their specific point signature.
         */
        void matchPoints(const GlarotSignature& gs, std::vector<std::pair<int, int> >& associations) const;

        /** Exports for gnuplot surface:
         *     theta0 rho0 signature00
         *     theta0 rho1 signature01
         *     ....
         *     theta0 rhoN signature0N
         *
         *     theta1 rho0 signature10
         *     theta1 rho1 signature11
         *     ...
         *     theta1 rhoN signature1N
         *     ...
         */
        void exportGnuplotSignature(std::ostream& out) const;

    protected:
        Eigen::MatrixXd signatureGlobal_;
        VectorPointSignature signaturePoints_;
        VectorPoint2d points_;
        Eigen::MatrixXd kernel_;
        double thetaStep_;
        double rhoStep_;
        Eigen::Matrix2d covar_;

        void computeParam(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, double& theta, double& rho) {
            rho = (p2 - p1).norm();
            if (p1.y() < p2.y()) {
                theta = atan2(p2.y() - p1.y(), p2.x() - p1.x());
            } else {
                theta = atan2(p1.y() - p2.y(), p1.x() - p2.x());
            }
            //    std::cout << " p1 " << p1.transpose() << " p2 " << p2.transpose() 
            //      << ": theta " << (180.0/M_PI*theta) << "  rho " << rho << std::endl;
        }

        void computeGaussianKernel();

        /** Adds the translated kernel to the GLARE grid. 
         */
        void addTranslatedKernel(Eigen::MatrixXd& grid, int itheta, int irho);

        /** Computes distance with shifted value.
         */
        double distanceL1Shift(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2, int shift) const;

        /**
         */
        double distanceL1Min(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2, int& shiftMin) const;
    };

} // end of namespace 

#endif