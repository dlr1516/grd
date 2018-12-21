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
#ifndef GRD_SKIP_EXTRACTOR_H
#define GRD_SKIP_EXTRACTOR_H

#include <iostream>
#include <queue>

#include <grd/LaserScan.h>

namespace grd {

    /**
     * Class CurvaturePointExtractor detects high curvature points using elementary 
     * geometric criteria like Discrete Curve Evolution (DCE) presented in several works like: 
     * 
     * Latecki, L.J., Lakamper, R. "Shape similarity measure based on correspondence 
     * of visual parts". IEEE T-PAMI 22(10), 1185–1190 (2000)
     * 
     */
    class SkipExtractor {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        static const double INFTY;

        enum ScoreType {
            SCORE_TRIANGLE_DIFF, SCORE_DCE
        };
    private:

        /** 
         * Data structure to store points and their previous and next point on a curve.  
         */
        struct Vertex {
            int currVertexId;
            int prevVertexId; // index of previous point on curve
            int nextVertexId; // index of next point on curve
            double score; // score value (high score for meaningful keypoints)
            bool isGap; // gap due to occlusions
            bool changed; // flag used internally by the algorithm 
            bool valid; // removed points have valid==true

            /**
             * Sorts the edges according to their score
             * @param e the other edge to compare
             */
            bool operator<(const Vertex& e) {
                return (score < e.score);
            }
        };
        typedef std::vector<Vertex> VertexVector;

        /**
         * Sorting data structure for pointer to edges according to their scores. 
         * Used in the priority queue to remove points with lower score
         */
        struct VertexPtrSort {

            bool operator()(Vertex* e1, Vertex* e2) {
                if (e1 == 0) {
                    return true;
                } else if (e2 == 0) {
                    return false;
                } else {
                    return (e1->score > e2->score);
                }
            }
        };
        typedef std::priority_queue<Vertex*, std::vector<Vertex*>, VertexPtrSort> VertexPtrPriorityQueue;

    public:

        /**
         * Default constructor.
         */
        SkipExtractor();

        /**
         * Destructor. 
         */
        virtual ~SkipExtractor();

        /**
         * 
         */
        void setScoreType(ScoreType st, double thres) {
            scoreType_ = st;
            scoreThres_ = thres;
        }

        void setGapThres(double q, double m) {
            distThresolhQ_ = q;
            distThresolhM_ = m;
        }


        /**
         * Extract simple corner keypoints from a given scan. 	 
         * @param scan input laser scan
         * @param keypoints corner points
         */
        virtual void extract(const LaserScan& scan, VectorPoint2d& keypoints);

        virtual void extract(const LaserScan& scan, VectorPoint2d& keypoints, VectorPoint2d& gaps);

    private:
        ScoreType scoreType_;
        double distThresolhQ_;
        double distThresolhM_;
        double scoreThres_;

        double evaluateScore(const VectorPoint2d& points, int iprev, int icurr, int inext);

        double evaluateScoreEucl(const VectorPoint2d& points, int iprev, int icurr, int inext);

        double evaluateScoreDCE(const VectorPoint2d& points, int iprev, int icurr, int inext);
    };

} // end of namespace 

#endif 