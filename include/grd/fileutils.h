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
#ifndef GRD_FILEUTILS_H
#define GRD_FILEUTILS_H

#include <cinttypes>
#include <cmath>
#include <string>
#include <vector>
#include <boost/filesystem.hpp> 


namespace grd {

    // ----------------------------------------------
    // OPERATIONS ON BITS
    // ----------------------------------------------

    /** Computes Morton index (see Z-order) for two indices.
     * From: http://stackoverflow.com/questions/12157685/z-order-curve-coordinates
     */
    inline std::uint32_t morton2d(std::uint16_t xpos, std::uint16_t ypos) {
        std::uint32_t x = xpos; // Interleave lower 16 bits of x and y, so the bits of x
        std::uint32_t y = ypos;
        std::uint32_t z;

        x = (x | (x << 8)) & 0x00FF00FF;
        x = (x | (x << 4)) & 0x0F0F0F0F;
        x = (x | (x << 2)) & 0x33333333;
        x = (x | (x << 1)) & 0x55555555;

        y = (y | (y << 8)) & 0x00FF00FF;
        y = (y | (y << 4)) & 0x0F0F0F0F;
        y = (y | (y << 2)) & 0x33333333;
        y = (y | (y << 1)) & 0x55555555;

        z = x | (y << 1);
        return z;
    }

    /** Computes Morton index (see Z-order) for two indices.
     * Version for 32 bits indices. 
     */
    //inline std::uint64_t morton2d(std::uint32_t xpos,std::uint32_t ypos)
    //{
    //  std::uint16_t xl = (xpos & 0x0000FFFF);
    //  std::uint16_t xu = (xpos & 0xFFFF0000) >> 16;
    //  std::uint16_t yl = (ypos & 0x0000FFFF);
    //  std::uint16_t yu = (ypos & 0xFFFF0000) >> 16;
    //  std::uint64_t z = morton2d(xl,yl) | (morton2d(xu,yu) << 32);
    //  return z;
    //}

    /** Interleaving of bits for z-order.
     */
    inline std::uint32_t shuffle(std::uint16_t x, std::uint16_t y) {
        std::uint32_t z = (y << 16) | x;
        std::uint32_t t;
        t = (z ^ (z >> 8)) & 0x0000FF00;
        z = z ^ t ^ (t << 8);
        t = (z ^ (z >> 4)) & 0x00F000F0;
        z = z ^ t ^ (t << 4);
        t = (z ^ (z >> 2)) & 0x0C0C0C0C;
        z = z ^ t ^ (t << 2);
        t = (z ^ (z >> 1)) & 0x22222222;
        z = z ^ t ^ (t << 1);
        return z;
    }

    /** Interleaving of bits for z-order.
     */
    //inline std::uint64_t shuffle(std::uint32_t x,std::uint32_t y)
    //{
    //  std::uint16_t xl = (x & 0x0000FFFF);
    //  std::uint16_t xu = (x & 0xFFFF0000) >> 16;
    //  std::uint16_t yl = (y & 0x0000FFFF);
    //  std::uint16_t yu = (y & 0xFFFF0000) >> 16;
    //  std::uint64_t z = shuffle(xl,yl) | (shuffle(xu,yu) << 32);
    //  return z;
    //}

    /** De-Interleaving of bits for z-order.
     */
    inline void unshuffle(std::uint32_t z, std::uint16_t& x, std::uint16_t& y) {
        std::uint32_t t;
        t = (z ^ (z >> 1)) & 0x22222222;
        z = z ^ t ^ (t << 1);
        t = (z ^ (z >> 2)) & 0x0C0C0C0C;
        z = z ^ t ^ (t << 2);
        t = (z ^ (z >> 4)) & 0x00F000F0;
        z = z ^ t ^ (t << 4);
        t = (z ^ (z >> 8)) & 0x0000FF00;
        z = z ^ t ^ (t << 8);
        x = z & 0x0000FFFF;
        y = (z & 0x0000FFFF) >> 16;
    }

    /** De-Interleaving of bits for z-order.
     */
    //inline void unshuffle(std::uint64_t z,std::uint32_t& x,std::uint32_t& y)
    //{
    //  std::uint16_t xl = (x & 0x0000FFFF);
    //  std::uint16_t xu = (x & 0xFFFF0000) >> 16;
    //  std::uint16_t yl = (y & 0x0000FFFF);
    //  std::uint16_t yu = (y & 0xFFFF0000) >> 16;
    //  std::uint64_t z = shuffle(xl,yl) | (shuffle(xu,yu) << 32);
    //}

    // ----------------------------------------------
    // GEOMETRY
    // ----------------------------------------------

    /** Computes the line obtained from linear regression on the given points. 
     */
    template <typename It>
    bool computeLine(It pbeg, It pend, double& theta, double& rho) {
        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;
        double sx = 0.0;
        double sy = 0.0;
        int num = 0;
        for (It pit = pbeg; pit != pend; ++pit) {
            sxx += pit->x() * pit->x();
            syy += pit->y() * pit->y();
            sxy += pit->x() * pit->y();
            sx += pit->x();
            sy += pit->y();
            ++num;
        }
        if (num < 1) return false;

        double msxx = (sxx - sx * sx / num) / num;
        double msyy = (syy - sy * sy / num) / num;
        double msxy = (sxy - sx * sy / num) / num;
        double b = 2.0 * msxy;
        double a = msxx;
        double c = msyy;
        theta = 0.5 * (atan2(b, a - c) + M_PI);
        theta = atan2(sin(theta), cos(theta));
        rho = (sx * cos(theta) + sy * sin(theta)) / num;
        if (rho < 0) {
            theta = atan2(sin(theta + M_PI), cos(theta + M_PI));
            rho = -rho;
        }
        return true;
    }

    // ----------------------------------------------
    // I/O OPERATIONS
    // ----------------------------------------------

    /** Returns a list of files based on Unix-like GLOB. 
     */
    void glob(const std::string globPath, std::vector<std::string>& matchingFiles);

    /** Generates a filename dependent on date and time.
     */
    std::string generateStampedString(const std::string prefix = "", const std::string postfix = "");

    std::string getWithoutExtension(std::string filename);

    std::string getPrefix(std::string filename);

    std::string getShortName(std::string filename);

    std::string getLeafDirectory(std::string filename);

} // end of namespace 

#endif