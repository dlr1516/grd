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
#include <grd/functions.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/laguerre.hpp>

namespace grd {

    // --------------------------------------------------------
    // FUNCTIONS
    // --------------------------------------------------------

    void fastCosSin(double x, double& c, double& s) {
        constexpr double factor = 1.0 / (2.0 * M_PI);
        x *= factor;

        c = x - (0.25 + floor(x + 0.25));
        c *= 16.0 * (std::abs(c) - 0.5);
        c += 0.225 * c * (std::abs(c) - 1.0);

        s = x - floor(x + 0.5);
        s *= 16.0 * (0.5 - std::abs(s));
        s += 0.225 * s * (std::abs(s) - 1.0);
    }

    double evaluateFourier(const std::vector<double>& coeffs, double theta) {
        double val, cth2, sth2, cth, sth, ctmp, stmp;
        int n;

        if (coeffs.size() % 2 != 0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of coefficients must be even: found " << coeffs.size() << std::endl;
        }
        n = (coeffs.size() / 2) - 1;

        cth2 = cos(theta);
        sth2 = sin(theta);
        cth = 1.0;
        sth = 0.0;
        val = 0.0;

        for (int k = 0; k <= n; ++k) {
            val += coeffs[2 * k] * cth + coeffs[2 * k + 1] * sth;
            ctmp = cth2 * cth - sth2 * sth;
            stmp = sth2 * cth + cth2 * sth;
            cth = ctmp;
            sth = stmp;
        }
        return val;
    }

    double evaluateVonMises(double t, double mean, double kappa) {
        return exp(kappa * cos(t - mean)) / (2 * M_PI * boost::math::cyl_bessel_i(0, kappa));
    }

    double evaluateErlang(double x, int d, double lambda) {
        double erlangNorm = 1.0 / (pow(lambda, d) * std::tgamma(d));
        return erlangNorm * pow(x, d - 1) * exp(-x / lambda);
    }

    double evaluateBiasedRayleigh(double x, double mu, double sigma) {
        double ratio = mu / sigma;
        double biasrayConst = mu * sigma * sqrt(0.5 * M_PI) * (1.0 + std::erf(ratio / sqrt(2.0)))
                + sigma * sigma * exp(-0.5 * ratio * ratio);
        //    std::cout << __FILE__ << "," << __LINE__ << ": biasrayConst " << biasrayConst
        //            << ", 1 / biasrayConst " << (1.0 / biasrayConst) << std::endl;
        double arg = (x - mu) / sigma;
        return x * exp(-0.5 * arg * arg) / biasrayConst;
    }

    double evaluateLaguerre(const std::vector<double>& coeffs, double r) {
        double val = 0.0;
        for (int k = 0; k <= coeffs.size(); ++k) {
            val += coeffs[k] * boost::math::laguerre(k, r);
        }
        return val;
    }

    // --------------------------------------------------------
    // COEFFICIENTS: FOURIER, LAGUERRE, ETC. 
    // --------------------------------------------------------

    void correlationFourierCoeffs(const std::vector<double>& fourierSrc, const std::vector<double>& fourierDst, std::vector<double>& fourierCor) {
        int n;

        if (fourierSrc.size() % 2 != 0 || fourierDst.size() % 2 != 0 || fourierSrc.size() != fourierDst.size()) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of ARSF coefficients must be even and equal: fourierSrc " << fourierSrc.size()
                    << ", fourierDst " << fourierDst.size() << std::endl;
            return;
        }
        n = (fourierSrc.size() / 2) - 1;

        // Resizes the correlation coefficients, if required
        if (fourierCor.size() != 2 * n + 2) {
            fourierCor.resize(2 * n + 2);
        }

        // Computes the coefficients
        for (int k = 0; k <= n; ++k) {
            fourierCor[2 * k] = 0.5 * (fourierSrc[2 * k] * fourierDst[2 * k] + fourierSrc[2 * k + 1] * fourierDst[2 * k + 1]);
            fourierCor[2 * k + 1] = 0.5 * (fourierSrc[2 * k] * fourierDst[2 * k + 1] - fourierSrc[2 * k + 1] * fourierDst[2 * k]);
        }
    }

    void besselIRatio(int n, double x, std::vector<double>& besratio) {
        double factor, seqPrev, seqCurr, seqNext;
        if (besratio.size() < n + 1) {
            besratio.resize(n + 1);
        }

        if (x < 0.0) x = -x;

        // If x~=0, then BesselI(0,x) = 1.0 and BesselI(k,x) = 0.0 for k > 0.
        // Thus, PNEBI(0,x) = 2.0 and PNEBI(k,x) = 0.0 for k > 0.
        if (x < 1e-6) {
            std::fill(besratio.begin(), besratio.end(), 0.0);
            besratio[0] = 2.0;
            return;
        }

        // Computes bessel function using back recursion
        factor = 2.0 / x;
        seqPrev = 0.0; // bip
        seqCurr = 1.0; // bi
        seqNext = 0.0; // bim
        for (int k = 2 * (n + (int) sqrt(40.0 * n)); k >= 0; --k) {
            seqNext = seqPrev + factor * k * seqCurr;
            seqPrev = seqCurr;
            seqCurr = seqNext;
            if (k <= n) {
                besratio[k] = seqPrev;
            }
            // To avoid overflow!
            if (seqCurr > BIG_NUM) {
                seqPrev *= SMALL_NUM;
                seqCurr *= SMALL_NUM;
                for (int i = 0; i < besratio.size(); ++i) {
                    besratio[i] *= SMALL_NUM;
                }
                //std::cerr << __FILE__ << "," << __LINE__ << ": ANTI-OVERFLOW!" << std::endl;
            }
        }

        double scaleFactor = 1.0 / (M_PI * besratio[0]);
        for (int i = 0; i < besratio.size(); ++i) {
            besratio[i] = scaleFactor * besratio[i];
        }
        besratio[0] *= 0.5;
    }

    void coeffFourierVonMises(int n, double thetaMean, double kappa, std::vector<double>& coeffs) {
        std::vector<double> besratio;

        if (coeffs.size() < 2 * n + 2) {
            coeffs.resize(2 * n + 2);
        }

        besselIRatio(n, kappa, besratio);
        for (int i = 0; i <= n; ++i) {
            coeffs[2 * i] = besratio[i] * cos(i * thetaMean);
            coeffs[2 * i + 1] = besratio[i] * sin(i * thetaMean);
        }
    }

    double coeffLaguerreErlang(int n, double lambda, int d) {
        double val = 0.0;
        double val2 = 0.0;
        double base = lambda / (lambda + 1.0);
        double sign = 1.0;
        double binom1 = 1.0;
        double binom2 = 1.0;
        double basePow = 1.0;


        for (int k = 0; k <= n; ++k) {
            val += sign * binom1 * binom2 * basePow;

            //        std::cout << "n " << n << "  k " << k << ": "
            //                << "binom1 " << binom1 << ", (n, k) " << (std::tgamma(n+1) / (std::tgamma(k+1) * std::tgamma(n - k + 1))) << std::endl;
            //        std::cout << "k " << k << "  d " << d << ": "
            //                << "binom2 " << binom2 << ", (k + d - 1, k) " << (std::tgamma(k + d) / (std::tgamma(k+1) * std::tgamma(d))) << std::endl;

            sign = -sign;
            // recurrent relation: binom(n, k+1) = binom(n, k) * (n - k) / (k+1)
            binom1 = binom1 * (n - k) / (k + 1);
            // recurrent relation: binom(k+d,k+1) = binom(k+d-1,k) * (k+d) / (k+1)
            binom2 = binom2 * (k + d) / (k + 1);
            // basePow is base^k
            basePow = basePow * base;
        }
        val = val / pow(lambda + 1.0, d);
        return val;
    }

    void coeffLaguerreErlang(int n, double lambda, int d, std::vector<double>& coeffs) {
        if (coeffs.size() < n + 1) {
            coeffs.resize(n + 1);
        }
        for (int i = 0; i <= n; ++i) {
            coeffs[i] = coeffLaguerreErlang(i, lambda, d);
        }
    }

    void umomentBiasray(int n, double mu, double sigma, std::vector<double>& moments) {
        double sigma2 = sigma * sigma;
        double ratio = -mu / sigma;
        double g0 = sqrt(0.5 * M_PI) * (1.0 + std::erf(ratio / sqrt(2.0)));
        double g1 = exp(-0.5 * ratio * ratio);
        double g2 = ratio * g1 + g0;

        moments.resize(n + 1);
        if (moments.size() >= 1) {
            moments[0] = mu * sigma * g0 + sigma2 * g1;
        }
        if (moments.size() >= 2) {
            moments[1] = sigma * mu * mu * g0 + 2.0 * sigma2 * mu * g1 + sigma * sigma2 * g2;
            for (int k = 2; k < moments.size(); ++k) {
                moments[k] = mu * moments[k - 1] + k * sigma2 * moments[k - 2];
            }
        }
    }

    void umomentBiasrayFactorial(int n, double mu, double sigma, std::vector<double>& momentsFact) {
        double sigma2 = sigma * sigma;
        double ratio = -mu / sigma;
        double g0 = sqrt(0.5 * M_PI) * (1.0 - std::erf(ratio / sqrt(2.0)));
        double g1 = exp(-0.5 * ratio * ratio);
        double g2 = ratio * g1 + g0;

        momentsFact.resize(n + 1);
        if (momentsFact.size() >= 1) {
            momentsFact[0] = mu * sigma * g0 + sigma2 * g1;
        }
        if (momentsFact.size() >= 2) {
            momentsFact[1] = (sigma * mu * mu * g0 + 2.0 * sigma2 * mu * g1 + sigma * sigma2 * g2);
            for (int k = 2; k < momentsFact.size(); ++k) {
                momentsFact[k] = mu * momentsFact[k - 1] / k + sigma2 * momentsFact[k - 2] / (k - 1);
            }
        }
    }

    //void umomentFactorialBiasray(int n, double mu, double sigma, std::vector<double>& moments) {
    //    double sigma2 = sigma * sigma;
    //    double ratio = -mu / sigma;
    //    double g0 = sqrt(0.5 * M_PI) * (1.0 - std::erf(ratio / sqrt(2.0)));
    //    double g1 = exp(-0.5 * ratio * ratio);
    //    double g2 = ratio * g1 + g0;
    //
    //    moments.resize(n + 1);
    //    if (moments.size() >= 1) {
    //        moments[0] = mu * sigma * g0 + sigma2 * g1;
    //    }
    //    if (moments.size() >= 2) {
    //        moments[1] = sigma * mu * mu * g0 + 2.0 * sigma2 * mu * g1 + sigma * sigma2 * g2;
    //        for (int k = 2; k < moments.size(); ++k) {
    //            moments[k] = mu * moments[k - 1] + k * sigma2 * moments[k - 2];
    //        }
    //    }
    //}
    //
    //void umomentFactorialBiasray(int n, double mu, double sigma, std::vector<double>& moments) {
    //    double sigma2 = sigma * sigma;
    //    double ratio = -mu / sigma;
    //    double g0 = sqrt(0.5 * M_PI) * (1.0 - std::erf(ratio / sqrt(2.0)));
    //    double g1 = exp(-0.5 * ratio * ratio);
    //    double g2 = ratio * g1 + g0;
    //
    //    moments.resize(n + 1);
    //    if (moments.size() >= 1) {
    //        moments[0] = mu * sigma * g0 + sigma2 * g1;
    //    }
    //    if (moments.size() >= 2) {
    //        moments[1] = sigma * mu * mu * g0 + 2.0 * sigma2 * mu * g1 + sigma * sigma2 * g2;
    //        for (int k = 2; k < moments.size(); ++k) {
    //            moments[k] = mu * moments[k - 1] + k * sigma2 * moments[k - 2];
    //        }
    //    }
    //}

    void coeffLaguerreBiasray(int n, double mu, double sigma, std::vector<double>& coeffs) {
        std::vector<double> moments;
        double sigma2 = sigma * sigma;
        double ratio = -mu / sigma;
        double biasConst = mu * sigma * sqrt(0.5 * M_PI) * (1.0 - std::erf(ratio / sqrt(2.0))) + sigma2 * exp(-0.5 * ratio * ratio);
        double ea = exp(0.5 * sigma2 - mu);
        double sign, binom, fact;

        //std::cout << __FILE__ << "," << __LINE__ << ": biasConst " << biasConst << std::endl;

        if (coeffs.size() < n + 1) {
            coeffs.resize(n + 1);
        }
        umomentBiasray(n, mu - sigma*sigma, sigma, moments);

        for (int k = 0; k <= n; ++k) {
            coeffs[k] = 0.0;
            sign = 1.0;
            binom = 1.0;
            fact = 1.0;
            for (int i = 0; i <= k; ++i) {
                coeffs[k] += sign * binom * moments[i] / fact;

                //            std::cout << "k " << k << "  i " << i << ": "
                //                    << "binom " << binom << ", (k, i) " << (std::tgamma(k+1) / (std::tgamma(i+1) * std::tgamma(k - i + 1))) << std::endl;
                //            std::cout << " factorial " << i << "! = " << fact << std::endl;

                sign = -sign;
                binom = binom * (k - i) / (i + 1);
                fact = fact * (i + 1);
            }
            coeffs[k] = ea / biasConst * coeffs[k];
        }
    }

    void coeffLaguerreBiasrayFactorial(int n, double mu, double sigma, std::vector<double>& coeffs) {
        std::vector<double> moments;
        double sigma2 = sigma * sigma;
        double ratio = -mu / sigma;
        double biasConst = mu * sigma * sqrt(0.5 * M_PI) * (1.0 - std::erf(ratio / sqrt(2.0))) + sigma2 * exp(-0.5 * ratio * ratio);
        double ea = exp(0.5 * sigma2 - mu);
        double sign, binom, fact;

        std::cout << "biasConst " << biasConst << std::endl;

        if (coeffs.size() < n + 1) {
            coeffs.resize(n + 1);
        }
        umomentBiasrayFactorial(n, mu - sigma*sigma, sigma, moments);

        for (int k = 0; k <= n; ++k) {
            coeffs[k] = 0.0;
            sign = 1.0;
            binom = 1.0;
            //        fact = 1.0;
            for (int i = 0; i <= k; ++i) {
                coeffs[k] += sign * binom * moments[i];

                //            std::cout << "k " << k << "  i " << i << ": "
                //                    << "binom " << binom << ", (k, i) " << (std::tgamma(k+1) / (std::tgamma(i+1) * std::tgamma(k - i + 1))) << std::endl;
                //            std::cout << " factorial " << i << "! = " << fact << std::endl;

                sign = -sign;
                binom = binom * (k - i) / (i + 1);
                //            fact = fact * (i + 1);
            }
            coeffs[k] = ea / biasConst * coeffs[k];
        }
    }


    // --------------------------------------------------------
    // INTERVAL FUNCTIONS
    // --------------------------------------------------------

    void findLUCos(double a, double b, double& cmin, double& cmax) {
        double amod, bmod;

        if (a > b) std::swap(a, b);

        if (b - a >= 2.0 * M_PI) {
            cmin = -1.0;
            cmax = +1.0;
        } else {
            // Normalizes to circular interval [0, 2*M_PI[
            amod = fmod(a, 2.0 * M_PI);
            if (amod < 0.0) amod += 2.0 * M_PI;
            bmod = fmod(b, 2.0 * M_PI);
            if (bmod < 0.0) bmod += 2.0 * M_PI;
            // Case bmod < amod: for example [300,30[ deg: angle 0 is included.
            if (bmod < amod) {
                cmax = +1.0;
                if (bmod < M_PI && M_PI < amod) {
                    cmin = std::min(cos(amod), cos(bmod));
                } else {
                    cmin = -1.0;
                }
                //      if (M_PI < bmod || amod < M_PI) {
                //        cmin = -1.0;
                //      }
                //      else {
                //        cmin = std::min(cos(amod),cos(bmod));
                //      }
            } else {
                cmax = std::max(cos(amod), cos(bmod));
                if (amod < M_PI && M_PI < bmod) {
                    cmin = -1.0;
                } else {
                    cmin = std::min(cos(amod), cos(bmod));
                }
            }
        }
    }

    void findLUFourier(const std::vector<double>& coeffs, double theta0, double theta1, double& fourierMin, double& fourierMax) {
        double amplitude, phase, sinusoidMin, sinusoidMax;
        int n, i0, i1;

        if (coeffs.size() % 2 != 0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of coefficients must be even: found " << coeffs.size() << std::endl;
        }
        n = (coeffs.size() / 2) - 1;

        if (theta1 < theta0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid interval [" << theta0 << "," << theta1 << "]: swapping endpoints to continue" << std::endl;
            std::swap(theta0, theta1);
        }

        // fourierMin and fourierMax initialized with constant component
        fourierMin = coeffs[0];
        fourierMax = coeffs[0];
        for (int k = 1; k <= n; ++k) {
            // t_k = a_k * cos(2*k*theta) + b_k * sin(2*k*theta) = amplitude * cos(2*k*theta - phase)
            // Period of k-th terms is M_PI / k. 
            amplitude = sqrt(coeffs[2 * k] * coeffs[2 * k] + coeffs[2 * k + 1] * coeffs[2 * k + 1]);
            phase = atan2(coeffs[2 * k + 1], coeffs[2 * k]);
            //std::cout << "k " << k << ", amplitude " << amplitude << ", phase[deg] " << (180.0/M_PI*phase) << std::endl;
            // If the [theta0,theta1] interval is larger than period, then the whole sinusoid amplitude is considered.
            // Otherwise, a more refined evaluation is performed.
            findLUCos(k * theta0 - phase, k * theta1 - phase, sinusoidMin, sinusoidMax);
            fourierMin += amplitude * sinusoidMin;
            fourierMax += amplitude * sinusoidMax;
        }
    }

} // enf og namespace 
