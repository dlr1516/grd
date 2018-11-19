#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <vector>
#include <cmath>


const double BIG_NUM = 1.0e+10;
const double SMALL_NUM = 1.0e-10;


// --------------------------------------------------------
// FUNCTIONS
// --------------------------------------------------------

/** Computes sine and cosine values using a parabolic approximation of sine.
 *    y = 16.0 * xn * (abs(xn) - 0.5)
 * where xn = x/(2*M_PI) is the normalized value of x over the period. 
 *
 * Code suggested here:
 *  http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648/6
 *  http://stackoverflow.com/questions/18662261/fastest-implementation-of-sine-cosine-and-square-root-in-c-doesnt-need-to-b
 */
void fastCosSin(double x, double& c, double& s);

/**
 * Computes the value of the given Fourier series at the given point theta. 
 * The user must provide the vector of serie coefficients:
 *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(i*x) + coeffs[2*i+1] * sin(i*x) )
 * @param coeffs the vector of coefficiens (vector size must be an even number!)
 * @param theta the point where to compute the Fourier serie value
 * @return the value of the function
 */
double evaluateFourier(const std::vector<double>& coeffs, double theta);

double evaluateLaguerre(const std::vector<double>& coeffs, double r);

/**
 * Computes the value of von Mises probability density function in t.
 * @param t value in which to compute the PDF value
 * @param mean mean angle
 * @param kappa concentration parameter
 * @return the PDF value
 */
double evaluateVonMises(double t, double mean, double kappa);

/**
 * Computes the value of Erlang probability density function in x.
 * @param x value in which to compute the PDF value
 * @param d order of Erlang
 * @param lambda rate parmater 
 * @return the PDF value
 */
double evaluateErlang(double x, int d, double lambda);

/**
 * Computes the value of Biased Rayleigh probability density function in x.
 * @param x value in which to compute the PDF value
 * @param mu mode of the distribution
 * @param sigma width of the distribution
 * @return the PDF value
 */
double evaluateBiasedRayleigh(double x, double mu, double sigma);

// --------------------------------------------------------
// COEFFICIENTS: FOURIER, LAGUERRE, ETC. 
// --------------------------------------------------------

/**
 * Computes the coefficients Fourier series resulting from the correlation of two source and 
 * destination Fourier series. 
 * @param fourierSrc coefficients of source Fourier series
 * @param fourierDst coefficients of destination Fourier series
 * @param fourierCor coefficients of correlation Fourier series
 */
void correlationFourierCoeffs(const std::vector<double>& fourierSrc, const std::vector<double>& fourierDst, std::vector<double>& fourierCor);

/**
 * Computes the vector of values besratio[] (modified Bessel function ratio) equal to:
 * 
 *   besratio[0] = 1 / (2.0 * M_PI);
 *   besratio[i] = besseli(n, x) / (M_PI * besseli(0, x); 
 * 
 * for i = 1, ..., n.
 * 
 * @param n
 * @param x
 * @param bratio
 */
void besselIRatio(int n, double x, std::vector<double>& besratio);

/**
 * Computes the Fourier coeffients for given von Mises distribution. 
 * @param n
 * @param thetaMean
 * @param kappa
 * @param coeffs
 */
void coeffFourierVonMises(int n, double thetaMean, double kappa, std::vector<double>& coeffs);

/**
 * Returns the coefficient with order n-th of Laguerre polynomial expansion for
 * an Erlang distribution. 
 * @param n
 * @param lamda
 * @param d
 * @return 
 */
double coeffLaguerreErlang(int n, double lambda, int d);

/**
 * Returns all the coefficients with orders from 0 to n of Laguerre polynomial 
 * expansion for an Erlang distribution. 
 * @param n
 * @param lamda
 * @param d
 * @param coeffs
 */
void coeffLaguerreErlang(int n, double lambda, int d, std::vector<double>& coeffs);

/**
 * Computes the unormalized moments of Biased Rayleigh distribution. 
 * @param n
 * @param mu
 * @param sigma
 * @param moments
 */
void umomentBiasray(int n, double mu, double sigma, std::vector<double>& moments);

/**
 * Computes the unormalized moments (order up to n) of Biased Rayleigh distribution 
 * divided by their order. 
 * 
 * 
 * 
 * @param n
 * @param mu
 * @param sigma
 * @param moments
 */
void umomentBiasrayFactorial(int n, double mu, double sigma, std::vector<double>& momentsFact);

/**
 * Returns all the coefficients with orders from 0 to n of Laguerre polynomial 
 * expansion for a Biased Rayleigh distribution. 
 * @param n
 * @param mu
 * @param sigma
 * @param coeffs
 */
void coeffLaguerreBiasray(int n, double mu, double sigma, std::vector<double>& coeffs);

/**
 * Returns all the coefficients with orders from 0 to n of Laguerre polynomial 
 * expansion for a Biased Rayleigh distribution. 
 * @param n
 * @param mu
 * @param sigma
 * @param coeffs
 */
void coeffLaguerreBiasrayFactorial(int n, double mu, double sigma, std::vector<double>& coeffs);

// --------------------------------------------------------
// INTERVAL FUNCTIONS
// --------------------------------------------------------

/** Computes lower and upper bounds of cosine function on a given interval.
 */
void findLUCos(double a, double b, double& cmin, double& cmax);

/** Computes lower and upper bounds of Fourier Series (represented by its coefficients)
 * on a given interval.
 * The vector of coefficients coeffs[i] are used in Fourier series:
 *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(2*i*x) + coeffs[2*i+1] * sin(2*i*x) )
 */
void findLUFourier(const std::vector<double>& coeffs, double theta0, double theta1, double& fourierMin, double& fourierfMax);

#endif /* FUNCTIONS_H */


