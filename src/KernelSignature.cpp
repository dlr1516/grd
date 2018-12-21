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
#include <vector>
#include <algorithm>    // std::min

#include <grd/KernelSignature.h>
#include <grd/functions.h>

namespace grd {

    KernelSignature::KernelSignature() : kernels_() {
    }

    KernelSignature::~KernelSignature() {
    }

    void KernelSignature::setPoints(const VectorPoint2d& points) {
        Point2d diff;
        KernelParam kparam;
        int pairNum = points.size() * (points.size() - 1) / 2;
        int pairIdx;

        kernels_.resize(points.size() * (points.size() - 1) / 2);

        // For each pair of points it computes their geometric relation in the form 
        // of a distribution. 
        // The distribution is represented by coefficients of orthogonal basis function
        // (Fourier for angles, Laguerre polynomials for ranges). 
        pairIdx = 0;
        for (int i = 0; i < points.size(); ++i) {
            for (int j = i + 1; j < points.size(); ++j) {
                diff = points[i] - points[j];
                kernels_[pairIdx].range = diff.norm();
                kernels_[pairIdx].theta = atan2(diff.y(), diff.x());
                if (kernels_[pairIdx].theta < 0.0)
                    kernels_[pairIdx].theta += M_PI;
                pairIdx++;
            }
        }
    }

    void KernelSignature::displayKernels(std::ostream& out) const {
        out << __FILE__ << "," << __LINE__ << ": modes/kernels:\n";
        for (auto& ker : kernels_) {
            out << "  angle[deg] " << (180.0 / M_PI * ker.theta) << "  range " << ker.range << std::endl;
        }
    }

    double KernelSignature::evaluate(double theta, double r) const {
        double val = 0.0;
        double distrTheta, distrRange;
        int erlangD;
        double erlangLambda;
        int kernelNum;

        kernelNum = 0;
        for (auto& ker : kernels_) {
            distrTheta = vonMisesNorm_ * exp(vonMisesKappa_ * cos(theta - ker.theta));
            // The range distribution is chosen among options.  
            switch (distribution_) {
                case ERLANG:
                    erlangD = std::min(50, (int) floor((ker.range / biasraySigma_) * (ker.range / biasraySigma_)));
                    erlangLambda = ker.range / erlangD;
                    distrRange = evaluateErlang(r, erlangLambda, erlangD);
                    break;
                case BIASED_RAYLEIGH:
                default:
                    distrRange = evaluateBiasedRayleigh(r, ker.range, biasraySigma_);
            }
            val += distrTheta * distrRange;
        }
        if (kernelNum > 0) {
            val = val / kernelNum;
        }

        return val;
    }

} // end of namespace