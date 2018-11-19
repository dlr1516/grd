/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   KernelSignature.h
 * Author: dario
 *
 * Created on 01 September 2018, 09:36
 */

#ifndef KERNELSIGNATURE_H
#define KERNELSIGNATURE_H

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

#include "Point.h"

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

#endif /* KERNELSIGNATURE_H */

