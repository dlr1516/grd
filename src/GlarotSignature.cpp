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
#include <grd/GlarotSignature.h>

namespace grd {

    GlarotSignature::GlarotSignature()
    : thetaStep_(M_PI / 8.0), rhoStep_(0.2) {
        // Initializes the grid and the default values of covariance matrix
        signatureGlobal_ = Eigen::MatrixXd::Zero(8, 50);
        setCovariance(thetaStep_, rhoStep_);
    }

    GlarotSignature::GlarotSignature(int thetaNum, int rhoNum, double rhoStep)
    : thetaStep_(M_PI / thetaNum), rhoStep_(rhoStep) {
        // Initializes the grid and the default values of covariance matrix
        signatureGlobal_ = Eigen::MatrixXd::Zero(thetaNum, rhoNum);
        setCovariance(thetaStep_, rhoStep_);
    }

    void GlarotSignature::clear() {
        signatureGlobal_.fill(0.0);
        signaturePoints_.clear();
    }

    void GlarotSignature::setSize(int thetaNum, int rhoNum, double rhoStep) {
        signatureGlobal_ = Eigen::MatrixXd::Zero(thetaNum, rhoNum);
        rhoStep_ = rhoStep;
    }

    void GlarotSignature::setPoints(const VectorPoint2d& points) {
        double rho, theta;
        int irho, itheta;
        // Creates the signature of each point
        PointSignature ps;
        signaturePoints_.reserve(points.size());
        for (int i = 0; i < points.size(); ++i) {
            ps.point = points[i];
            ps.signature = Eigen::MatrixXd::Zero(signatureGlobal_.rows(), signatureGlobal_.cols());
            signaturePoints_.push_back(ps);
        }
        // For each pair of points, it computes their angle (theta) and distance (rho).
        for (int i = 0; i < points.size(); ++i) {
            for (int j = i + 1; j < points.size(); ++j) {
                computeParam(points[i], points[j], theta, rho);
                itheta = (int) round(theta / thetaStep_);
                irho = (int) round(rho / rhoStep_);
                //      std::cout << "  (theta,rho) = (" << theta << "," << rho << ")   itheta " << itheta << " irho " << irho << std::endl;
                if (0 <= itheta && itheta < signatureGlobal_.rows() && 0 <= irho && irho < signatureGlobal_.cols()) {
                    addTranslatedKernel(signaturePoints_[i].signature, itheta, irho);
                    addTranslatedKernel(signaturePoints_[j].signature, itheta, irho);
                    signatureGlobal_ += signaturePoints_[i].signature + signaturePoints_[j].signature;
                }
            }
        }
        // Normalizes the global GLARE signature according to the max value 
        ////  std::cout << __FILE__ << "," << __LINE__ << ": signatureGlobal_ " << signatureGlobal_.rows() 
        //    << "*" << signatureGlobal_.cols() << std::endl;
        double valMax = signatureGlobal_.maxCoeff();
        assert(!std::isnan(valMax));
        if (!std::isnan(valMax) && std::abs(valMax) > 1e-4) {
            signatureGlobal_ = signatureGlobal_ / valMax;
        }
        //  for (int r = 0; r < signatureGlobal_.rows(); ++r) {
        //    for (int c = 0; c < signatureGlobal_.cols(); ++c) {
        //      if (std::isnan(signatureGlobal_(r,c))) {
        //        std::cout << __FILE__ << "," << __LINE__ << ":\n" << signatureGlobal_ << std::endl;
        //        assert(0);
        //      }
        //    }
        //  }
        points_.clear();
        points_.reserve(points.size());
        std::copy(points.begin(), points.end(), std::back_inserter(points_));
    }

    double GlarotSignature::distanceL1(const GlarotSignature& gs) const {
        //  double dist = (signatureGlobal_ - gs.getSignatureGlobal()).lpNorm<1>();
        double dist = distanceL1Shift(signatureGlobal_, gs.getSignatureGlobal(), 0);
        return dist;
    }

    double GlarotSignature::distanceL1Min(const GlarotSignature& gs) const {
        //  double dist, distMin;
        //  distMin = 1e+6;
        //  for (int i = 0; i < signatureGlobal_.rows(); ++i) {
        //    dist = distanceL1Shift(signatureGlobal_,gs.getSignatureGlobal(),i);
        //    if (dist < distMin) distMin = dist;
        //  }
        //  return distMin;
        int imin;
        return distanceL1Min(signatureGlobal_, gs.getSignatureGlobal(), imin);
    }

    void GlarotSignature::matchPoints(const GlarotSignature& gs, std::vector<std::pair<int, int> >& associations) const {
        // Computes matrix of distances
        int n1 = getPointNum();
        int n2 = gs.getPointNum();
        int imin, jbest;
        double distBest;
        Eigen::MatrixXd dist = Eigen::MatrixXd::Zero(n1, n2);
        associations.clear();
        for (int i = 0; i < n1; ++i) {
            jbest = -1;
            for (int j = 0; j < n2; ++j) {
                dist(i, j) = distanceL1Min(getSignaturePoint(i), gs.getSignaturePoint(j), imin);
                //      dist(j,i) = dist(i,j);
                if (jbest < 0 || dist(i, j) < distBest) {
                    jbest = j;
                    distBest = dist(i, j);
                }
            }
            associations.push_back(std::make_pair(i, jbest));
        }
    }

    void GlarotSignature::exportGnuplotSignature(std::ostream& out) const {
        for (int i = 0; i < signatureGlobal_.rows(); ++i) {
            for (int j = 0; j < signatureGlobal_.cols(); ++j) {
                out << (thetaStep_ * i) << " " << (rhoStep_ * j) << " " << signatureGlobal_(i, j) << "\n";
            }
            out << "\n";
        }
    }

    // --------------------------------------------------------
    // PRIVATE METHODS
    // --------------------------------------------------------

    void GlarotSignature::computeGaussianKernel() {
        Eigen::Matrix2d covarInv;
        Eigen::Vector2d v;
        double normConst;
        int winTheta, winRho;
        // Computes the size of Gaussian kernel sample matrix according to its covariance matrix
        winTheta = 3 * (int) ceil(sqrt(covar_(0, 0)) / thetaStep_);
        winRho = 3 * (int) ceil(sqrt(covar_(1, 1)) / rhoStep_);
        kernel_ = Eigen::MatrixXd::Zero(2 * winTheta + 1, 2 * winRho + 1);
        // Sets the values of the mask (normalization of kernel function may not be required...)
        covarInv = covar_.inverse();
        //  normConst = 1.0 / (2 * M_PI * covar_.determinant());
        normConst = 1.0;
        for (int t = -winTheta; t <= winTheta; ++t) {
            for (int r = -winRho; r <= winRho; ++r) {
                v << thetaStep_*t, rhoStep_*r;
                //      std::cout << "t " << t << " r " << r << " winTheta " << winTheta << " winRho " << winRho << std::endl
                //        << "  t+winTheta: " << (t+winTheta) << " max " << kernel_.rows() 
                //        << ", r+winRho: " << (r+winRho) << " max " << kernel_.rows() << std::endl;
                kernel_(t + winTheta, r + winRho) = normConst * exp(-0.5 * (v.transpose() * covarInv * v)(0, 0));
            }
        }
        //  std::cout << __FILE__ << "," << __LINE__ << ": kernel\n" << kernel_ << std::endl;
    }

    void GlarotSignature::addTranslatedKernel(Eigen::MatrixXd& grid, int itheta, int irho) {
        int gt, gr, kt, kr;
        int winTheta = kernel_.rows() / 2;
        int winRho = kernel_.cols() / 2;
        for (int t = 0; t < kernel_.rows(); ++t) {
            for (int r = 0; r < kernel_.cols(); ++r) {
                gt = itheta + (t - winTheta);
                gr = irho + (r - winRho);
                //      std::cout << "  row (" << itheta << ") + (" << (t-winTheta) << ")  col (" << irho << ") + (" << (r-winRho) << ") " << std::endl;
                if (0 <= gt && gt < grid.rows() && 0 <= gr && gr < grid.cols()) {
                    grid(gt, gr) += kernel_(t, r);
                }
            }
        }
    }

    double GlarotSignature::distanceL1Shift(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2, int shift) const {
        int rowNum = m1.rows();
        int colNum = m1.cols();
        int rowShift;
        double dist;
        assert(m2.rows() == rowNum && m2.cols() == colNum);
        dist = 0.0;
        for (int r = 0; r < rowNum; ++r) {
            rowShift = (r + shift) % rowNum;
            for (int c = 0; c < colNum; ++c) {
                dist += std::abs(m1(r, c) - m2(rowShift, c));
            }
        }
        return dist;
    }

    double GlarotSignature::distanceL1Min(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2, int& shiftMin) const {
        double dist, distMin;
        distMin = distanceL1Shift(m1, m2, 0);
        shiftMin = 0;
        for (int i = 1; i < m1.rows(); ++i) {
            dist = distanceL1Shift(m1, m2, i);
            if (dist < distMin) {
                distMin = dist;
                shiftMin = i;
            }
        }
        return distMin;
    }

} // end of namespace