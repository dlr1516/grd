/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DistributionSignature.cpp
 * Author: dario
 * 
 * Created on 28 August 2018, 15:05
 */
#include <boost/math/special_functions/laguerre.hpp>

#include <grd/DistributionSignature.h>
#include <grd/functions.h>
#include <grd/BBOptimizer1d.h>

DistributionSignature::DistributionSignature()
: signature_(10, 20), points_(), thetaN_(10), rangeN_(20), distribution_(BIASED_RAYLEIGH),
vonMisesKappa_(1.0), biasraySigma_(1.0),
correlThetaTol_(0.0017453), correlFourierTol_(1.0), autoCorr_(1.0) {
    init();
}

DistributionSignature::DistributionSignature(const DistributionSignature& ds)
: signature_(ds.signature_), points_(ds.points_), thetaN_(ds.thetaN_), rangeN_(ds.rangeN_),
distribution_(ds.distribution_), vonMisesKappa_(ds.vonMisesKappa_), biasraySigma_(ds.biasraySigma_),
correlThetaTol_(ds.correlThetaTol_), correlFourierTol_(ds.correlFourierTol_), autoCorr_(ds.autoCorr_) {
    init();
}

DistributionSignature::DistributionSignature(int thetaN, int rangeN, RangeDistribution distr)
: signature_(thetaN, rangeN), points_(), thetaN_(thetaN), rangeN_(rangeN), distribution_(distr),
vonMisesKappa_(1.0), biasraySigma_(1.0),
correlThetaTol_(0.0017453), correlFourierTol_(1.0), autoCorr_(1.0) {
    init();
}

DistributionSignature::~DistributionSignature() {
}

void DistributionSignature::clear() {
    signature_.fill(0.0);
}

void DistributionSignature::setSize(int thetaN, int rangeN) {
    signature_.resize(thetaN, rangeN);
    thetaN_ = thetaN;
    rangeN_ = rangeN;

    // Recomputes the ratios of modified Bessel functions of first kind used to 
    // computes the Fourier coefficients of von Mises distribution. 
    //    vonMisesBesselRatios_[0] = 1 / (2.0 * M_PI);
    //    vonMisesBesselRatios_[i] = besseli(n, x) / (M_PI * besseli(0, x); 
    if (thetaN + 1 != vonMisesBesselRatios_.size()) {
        init();
    }
}

void DistributionSignature::setRangeDistribution(RangeDistribution distr) {
    distribution_ = distr;
    // Recomputes the ratios of modified Bessel functions of first kind used to 
    // computes the Fourier coefficients of von Mises distribution. 
    //    vonMisesBesselRatios_[0] = 1 / (2.0 * M_PI);
    //    vonMisesBesselRatios_[i] = besseli(n, x) / (M_PI * besseli(0, x); 
}

void DistributionSignature::setPoints(const VectorPoint2d& points) {
    Point2d diff;
    double range, angle;
    int pairNum = points.size() * (points.size() - 1);

    // Initializes the signature 
    signature_ = Eigen::MatrixXd::Zero(2 * thetaN_ + 2, rangeN_ + 1);

    //    std::cout << __FILE__ << "," << __LINE__ << ": vonMisesBesselRatios_: " << std::endl;
    //    for (int i = 0; i < vonMisesBesselRatios_.size(); ++i) {
    //        std::cout << "vonMisesBesselRatios_ " << i << ": " << vonMisesBesselRatios_[i] << std::endl;
    //    }

    // For each pair of points it computes their geometric relation in the form 
    // of a distribution. 
    // The distribution is represented by coefficients of orthogonal basis function
    // (Fourier for angles, Laguerre polynomials for ranges). 
    //std::cout << __FILE__ << "," << __LINE__ << ": modes/kernels " << std::endl;
    for (int i = 0; i < points.size(); ++i) {
        for (int j = 0; j < points.size(); ++j) {
            if (i != j) {
                diff = points[i] - points[j];
                range = diff.norm();
                angle = atan2(diff.y(), diff.x());
                if (angle < 0.0)
                    angle += M_PI;
                //std::cout << "  angle[deg] " << (180.0 / M_PI * angle) << "  range " << range << std::endl;
                updateCoeffs(range, angle);
            }
        }
    }
    if (pairNum > 0) {
        signature_ *= 1.0 / pairNum;
        autoCorr_ = 1.0;     // autoCorr_ is used by correlate()): we start with 1.0
        autoCorr_ = sqrt(fabs(correlate(*this)));
    }
}

double DistributionSignature::evaluate(double theta, double r) const {
    Eigen::VectorXf sinusoids(2 * thetaN_ + 2);
    Eigen::VectorXf laguerrePolynomials(rangeN_ + 1);
    double val = 0.0;

    if (signature_.rows() != 2 * thetaN_ + 2 || signature_.cols() != rangeN_ + 1) {
        std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of signature "
                << signature_.rows() << " x " << signature_.cols()
                << " instead of " << (2 * thetaN_ + 2) << " x " << (rangeN_ + 1)
                << std::endl;
        return val;
    }

    // Fills sinusoids vector with values
    //   sinusoids(2 * i) = cos(i * theta)
    //   sinusoids(2 * i + 1) = sin(i * theta)
    // using recursive formula of cos(i*theta) and sin(i*theta):
    //   cos(i*theta) = cos(theta) * cos((i-1)*theta) - sin(theta) * sin((i-1)*theta)
    //   sin(i*theta) = sin(theta) * cos((i-1)*theta) + cos(theta) * sin((i-1)*theta)
    double cth = cos(theta);
    double sth = sin(theta);
    sinusoids(0) = 1.0;
    sinusoids(1) = 0.0;
    for (int i = 1; i <= thetaN_; ++i) {
        sinusoids(2 * i) = cth * sinusoids(2 * i - 2) - sth * sinusoids(2 * i - 1);
        sinusoids(2 * i + 1) = sth * sinusoids(2 * i - 2) + cth * sinusoids(2 * i - 1);
    }

    // Fills vector of Laguerre polynomials evaluated in r.
    // One could directly use the function:
    //    laguerrePolynomials(j) = boost::math::laguerre(j, r);
    // or use the (more efficient?) recurrent evaluation
    laguerrePolynomials(0) = boost::math::laguerre(0, r);
    if (rangeN_ > 0)
        laguerrePolynomials(1) = boost::math::laguerre(1, r);
    for (int j = 2; j <= rangeN_; ++j) {
        laguerrePolynomials(j) = boost::math::laguerre_next(j, r, laguerrePolynomials[j - 1], laguerrePolynomials(j - 2));
    }

    // Evaluates in point
    for (int i = 0; i <= thetaN_; ++i) {
        for (int j = 0; j <= rangeN_; ++j) {
            val += signature_(2 * i, j) * sinusoids(2 * i) * laguerrePolynomials(j)
                    + signature_(2 * i + 1, j) * sinusoids(2 * i + 1) * laguerrePolynomials(j);
        }
    }

    return val;
}

double DistributionSignature::correlate(const DistributionSignature& dsig, double& thetaShift, std::vector<double>& correlCoeffs) const {
    int thetaNum = std::min(thetaN_, dsig.thetaN_);
    int rangeNum = std::min(rangeN_, dsig.rangeN_);
    //std::vector<double> correlCoeffs(2 * thetaNum + 2);
    double as, bs, ad, bd, ag, bg, valueMax;

    correlCoeffs.resize(2 * thetaNum + 2);

    // Computes the coefficients of correlation  Fourier series
    Eigen::MatrixXd coeffTerms(2 * thetaNum + 2, rangeNum + 1);
    for (int i = 0; i <= thetaNum; ++i) {
        ag = 0.0;
        bg = 0.0;
        for (int j = 0; j <= rangeNum; ++j) {
            as = signature_(2 * i, j);
            bs = signature_(2 * i + 1, j);
            ad = dsig.signature_(2 * i, j);
            bd = dsig.signature_(2 * i + 1, j);
            ag += 0.5 * (as * ad + bs * bd);
            bg += 0.5 * (as * bd - bs * ad);
            coeffTerms(2 * i, j) = 0.5 * (as * ad + bs * bd);
            coeffTerms(2 * i + 1, j) = 0.5 * (as * bd - bs * ad);
        }
        correlCoeffs[2 * i] = ag;
        correlCoeffs[2 * i + 1] = bg;
        //std::cout << "ag[" << i << "] " << ag << "  bg[" << i << "] " << bg << std::endl;
    }
    //std::cout << "\ncoeffTerms:\n" << coeffTerms << "\n\n";

    findGlobalMaxBBFourier(correlCoeffs, 0.0, 2.0 * M_PI, correlThetaTol_, correlFourierTol_, thetaShift, valueMax);

    return (valueMax / (autoCorr_ * dsig.getAutoCorrel()));
}

double DistributionSignature::correlate(const DistributionSignature& dsig, double& thetaShift) const {
    std::vector<double> correlCoeffs;
    return correlate(dsig, thetaShift, correlCoeffs);
}

double DistributionSignature::correlate(const DistributionSignature& dsig) const {
    double thetaShift;
    return correlate(dsig, thetaShift);
}

// --------------------------------------------------------
// PRIVATE METHODS
// --------------------------------------------------------

void DistributionSignature::init() {
    besselIRatio(thetaN_, vonMisesKappa_, vonMisesBesselRatios_);

    //        std::cout << "init vonMisesBesselRatios_: thetaN_ " << thetaN_ << ", vonMisesKappa_ " << vonMisesKappa_ << std::endl;
    //        for (int i = 0; i < vonMisesBesselRatios_.size(); ++i) {
    //            std::cout << "vonMisesBesselRatios_ " << i << ": " << vonMisesBesselRatios_[i] << std::endl;
    //        }
}

void DistributionSignature::updateCoeffs(double range, double angle) {
    std::vector<double> thetaCoeffs(2 * thetaN_ + 2, 0.0);
    std::vector<double> rangeCoeffs(rangeN_ + 1, 0.0);
    double erlangLambda;
    int erlangD;

    if (vonMisesBesselRatios_.size() < thetaN_ + 1) {
        std::cerr << __FILE__ << "," << __LINE__ << ": vonMisesBesselRatios_.size() " << vonMisesBesselRatios_.size()
                << " instead of " << (thetaN_ + 1) << std::endl;
        return;
    }

    // Computes the Fourier coeffients of von Mises distribution
    for (int i = 0; i <= thetaN_; ++i) {
        thetaCoeffs[2 * i] = vonMisesBesselRatios_[i] * cos(i * angle);
        thetaCoeffs[2 * i + 1] = vonMisesBesselRatios_[i] * sin(i * angle);
    }

    // Computes the Laguerre polynomials coeffients of range distribution.
    // The range distribution is chosen among options.  
    switch (distribution_) {
        case ERLANG:
            erlangLambda = 1.0;
            erlangD = 2;
            coeffLaguerreErlang(rangeN_, erlangLambda, erlangD, rangeCoeffs);
            break;
        case BIASED_RAYLEIGH:
        default:
            coeffLaguerreBiasray(rangeN_, range, biasraySigma_, rangeCoeffs);
    }
    //std::cout << __FILE__ << "," << __LINE__ << ":  rangeCoeffs[0] " << rangeCoeffs[0] << std::endl;

    // Updates the signature coefficients
    for (int i = 0; i < thetaCoeffs.size(); ++i) {
        for (int j = 0; j < rangeCoeffs.size(); ++j) {
            signature_(i, j) += thetaCoeffs[i] * rangeCoeffs[j];
        }
    }
}
