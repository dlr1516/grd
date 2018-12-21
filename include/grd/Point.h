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
#ifndef GRD_POINT_H
#define GRD_POINT_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <vector>
#include <deque>

namespace grd {


    /**
     * Types for a 2D point in the form of Euclidean vector. 
     */
    typedef Eigen::Vector2d Point2d;
    typedef Eigen::Vector2f Point2f;
    typedef Eigen::Vector2i Point2i;

    /**
     * Types for a 2D matrix in the form of Euclidean vector. 
     */
    typedef Eigen::Matrix2d Matrix2d;
    typedef Eigen::Matrix2f Matrix2f;
    typedef Eigen::Matrix2i Matrix2i;

    /**
     * Types for STL vectors of different point types. 
     */
    typedef std::vector<Point2d, Eigen::aligned_allocator<Point2d> > VectorPoint2d;
    typedef std::vector<Point2f, Eigen::aligned_allocator<Point2f> > VectorPoint2f;
    typedef std::vector<Point2i, Eigen::aligned_allocator<Point2i> > VectorPoint2i;
    typedef std::deque<Point2d, Eigen::aligned_allocator<Point2d> > DequePoint2d;
    typedef std::deque<Point2f, Eigen::aligned_allocator<Point2f> > DequePoint2f;
    typedef std::deque<Point2i, Eigen::aligned_allocator<Point2i> > DequePoint2i;

    /**
     *  Type for reference frames
     */
    typedef Eigen::Isometry2d Transformation2d;
    typedef Eigen::Isometry2f Transformation2f;


} // end of namespace 

#endif
