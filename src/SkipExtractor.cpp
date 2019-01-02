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
#include <queue>
#include <grd/SkipExtractor.h>

namespace grd {

    const double SkipExtractor::INFTY = 1e+6;

    SkipExtractor::SkipExtractor()
    : scoreType_(SCORE_TRIANGLE_DIFF), distThresolhQ_(0.1), distThresolhM_(0.02), scoreThres_(0.4) { //0.4
    }

    SkipExtractor::~SkipExtractor() {
    }

    void SkipExtractor::extract(const LaserScan& scan, VectorPoint2d& keypoints) {
        VectorPoint2d gaps;
        extract(scan, keypoints, gaps);
    }

    void SkipExtractor::extract(const LaserScan& scan, VectorPoint2d& keypoints, VectorPoint2d& gaps) {
        VertexPtrPriorityQueue queue;
        Vertex *curr, *prev, *next;
        int n, iprev, inext;

        // Creates edge vector with indices to prev and next point 
        n = scan.points.size();
        VertexVector vertices(n);
        int countGap = 0;
        for (int i = 0; i < n; ++i) {
            iprev = (i + n - 1) % n;
            inext = (i + 1) % n;
            // Finds if there is a large range diff with prev range
            vertices[i].currVertexId = i;
            vertices[i].isGap = false;
            vertices[i].valid = true;
            if (fabs(scan.ranges[iprev] - scan.ranges[i]) < distThresolhQ_ + 0.5 * distThresolhM_ * (scan.ranges[i] + scan.ranges[iprev])) {
                vertices[i].prevVertexId = iprev;

            } else {
                vertices[i].prevVertexId = i;
                vertices[i].isGap = true;
            }
            // Finds if there is a large range diff with next range
            if (fabs(scan.ranges[inext] - scan.ranges[i]) < distThresolhQ_ + 0.5 * distThresolhM_ * (scan.ranges[i] + scan.ranges[inext])) {
                vertices[i].nextVertexId = inext;
            } else {
                vertices[i].nextVertexId = i;
                vertices[i].isGap = true;
            }
            // If there is no gap on left and on right computes the score
            if (vertices[i].prevVertexId != i && vertices[i].nextVertexId != i) {
                vertices[i].score = evaluateScore(scan.points, iprev, i, inext);
                queue.push(&vertices[i]);
            } else {
                vertices[i].score = INFTY;
                countGap++;
                //gaps.push_back(scan.points[i]);
            }
            vertices[i].changed = false;

            //            std::cout << "  point " << i << ": score " << edges[i].score
            //                    << " prev " << edges[i].prevVertexId << " next " << edges[i].nextVertexId
            //                    << ", changed " << edges[i].changed << std::endl;
        }
        //std::cout << "countGap " << countGap << ", queue: " << queue.size() << std::endl;
        //        if (!queue.empty()) {
        //            std::cout << "queue top " << queue.top()->score << std::endl;
        //        }

        // Dicrete Curve Evolution
        while (!queue.empty() && queue.top()->score < scoreThres_) {
            curr = queue.top();
            queue.pop();
            //assert(!curr->isGap);
            //std::cout << "curr " << curr->currVertexId << " score " << curr->score << std::endl;
            if (curr->changed) {
                // The score of the item top has not been updated: thus, the item 
                // is removed and pushed again in priority queue
                curr->score = evaluateScore(scan.points, curr->prevVertexId, curr->currVertexId, curr->nextVertexId);
                curr->changed = false;
                queue.push(curr);
            } else {
                curr->valid = false;
                // After removing curr from the link, its previous node is updated
                prev = &vertices.at(curr->prevVertexId);
                prev->nextVertexId = curr->nextVertexId;
                prev->changed = true;
                //                if (!prev->isGap) {
                //                    prev->score = evaluateScore(scan.points, prev->prevVertexId, curr->prevVertexId, prev->nextVertexId);
                //                }
                // After removing curr from the link, its next node is updated
                next = &vertices.at(curr->nextVertexId);
                next->prevVertexId = curr->prevVertexId;
                next->changed = true;
                //                if (!next->isGap) {
                //                    next->score = evaluateScore(scan.points, next->prevVertexId, curr->nextVertexId, next->nextVertexId);
                //                }
            }
        }
        //        std::cout << "after queue: " << queue.size() << std::endl;

        // Checks the number of items with flag changed sets to true (we expect them to be very limited in number)


        keypoints.clear();
        for (int i = 0; i < vertices.size(); ++i) {
            if (vertices[i].valid && !vertices[i].isGap) {
                //std::cout << "keypoint " << vertices[i].currVertexId << " " << scan.points[i].transpose() << " score " << vertices[i].score << std::endl;
                keypoints.push_back(scan.points[i]);
            } else if (vertices[i].isGap) {
                gaps.push_back(scan.points[i]);
            }
        }
    }

    double SkipExtractor::evaluateScore(const VectorPoint2d& points, int iprev, int icurr, int inext) {
        if (scoreType_ == SCORE_TRIANGLE_DIFF) {
            return evaluateScoreEucl(points, iprev, icurr, inext);
        } else {
            return evaluateScoreDCE(points, iprev, icurr, inext);
        }
    }

    double SkipExtractor::evaluateScoreEucl(const VectorPoint2d& points, int iprev, int icurr, int inext) {
        double lenPrev = (points[iprev] - points[icurr]).norm();
        double lenNext = (points[inext] - points[icurr]).norm();
        double lenTot = (points[inext] - points[iprev]).norm();
        return fabs(lenPrev + lenNext - lenTot);
    }

    double SkipExtractor::evaluateScoreDCE(const VectorPoint2d& points, int iprev, int icurr, int inext) {
        Point2d vecPrev = (points[iprev] - points[icurr]);
        Point2d vecNext = (points[inext] - points[icurr]);
        double lenPrev = vecPrev.norm();
        double lenNext = vecNext.norm();
        double angle = acos(vecPrev.dot(vecNext) / (vecPrev.norm() * vecNext.norm()));
        if (angle < 0.5 * M_PI) {
            angle = M_PI - angle;
        }
        return (angle * lenPrev * lenNext / (lenPrev + lenNext));
    }

} // end of namespace grd