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
#ifndef GRD_CARMEN_V2_READER_H
#define GRD_CARMEN_V2_READER_H

#include <iostream>
#include <string>
#include <fstream>

#include <grd/LaserScan.h>

namespace grd {

    struct RobotLaserMess {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        LaserScan scan;
        Eigen::Vector3d wTb;
        Eigen::Vector3d wTl;
        int vertexId;
    };

    class CarmenV2Reader {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        std::vector<RobotLaserMess> scans;

        CarmenV2Reader() {

        };

        CarmenV2Reader(const std::string& filename) {
            readCarmen(filename);
        }

        void readCarmen(const std::string& filename) {
            scans.clear();
            using namespace std;
            ifstream in(filename.c_str());
            if (!in.is_open()) {
                cerr << "Could not open file: " << filename << endl;
                exit(0);
            }

            cout << "Reading carmen file '" << filename << "'...." << endl;

            string line, temp;
            int lastVertexId;
            while (!in.eof()) {
                getline(in, line);
                if (line.length() > 0) {
                    stringstream ss(line);
                    ss >> temp;
                    if (temp == "VERTEX2") {
                        ss >> lastVertexId;
                    }
                    if (temp == "ROBOTLASER1") {
                        RobotLaserMess mess;
                        //                    LaserScan scan;
                        //                    scan.rangeMin = 0.1;
                        //                    scan.intensityMin = 0.0;
                        //                    scan.intensityMax = 100.0;
                        double angleMin, angleInc, rangeMax, fov, accuracy;
                        std::string sensorName, robotName;
                        std::vector<double> ranges;
                        std::vector<double> remission;
                        unsigned int remissionNumber = 0;
                        double timestamp;
                        int laserType, remissionMode, numBeams;
                        double tv, rv, forward_safety_dist, side_safety_dist, turn_axis;

                        ss >> laserType >> angleMin >> fov >> angleInc >> rangeMax >> accuracy >> remissionMode; // Laser sensor parameters
                        ss >> numBeams;

                        mess.scan.setAngleMin(angleMin);
                        mess.scan.setAngleInc(angleInc);
                        mess.scan.setLaserFoV(fov);
                        mess.scan.setNumBeams(numBeams);

                        //                    scan.ranges.resize(scan.numBeams);
                        //                    scan.points.resize(scan.numBeams);
                        //                    scan.intensities.resize(scan.numBeams, 0.0);
                        ranges.resize(numBeams);
                        for (int i = 0; i < numBeams; ++i) {
                            ss >> ranges[i];
                            //                        double rho = scan.ranges[i];
                            //                        double alpha = (i * scan.angleInc) + scan.angleMin;
                            //                        scan.points[i] = nav::Point2d(rho * cos(alpha), rho * sin(alpha));
                        }
                        mess.scan.fromRanges(ranges);

                        ss >> remissionNumber;
                        remission.resize(remissionNumber);

                        for (uint i = 0; i < remissionNumber; i++) {
                            ss >> remission[i];
                        }

                        ss >> mess.wTl[0]; // wTl x
                        ss >> mess.wTl[1]; // wTl y
                        ss >> mess.wTl[2]; // wTl theta
                        ss >> mess.wTb[0]; // wTb x
                        ss >> mess.wTb[1]; // wTb y
                        ss >> mess.wTb[2]; // wTb theta

                        ss >> tv >> rv >> forward_safety_dist >> side_safety_dist >> turn_axis;

                        ss >> timestamp >> robotName;

                        mess.scan.setTimestamp(timestamp);
                        mess.vertexId = lastVertexId;

                        scans.push_back(std::move(mess));

                        //                    std::cout << "Read vertex " << mess.vertexId << " in "
                        //                            << "[" << mess.wTb[0] << ", " << mess.wTb[1] << ", " << (180.0 / M_PI * mess.wTb[2]) << " deg]"
                        //                            << " scan with " << mess.scan.getNumBeams() 
                        //                            << std::endl;
                    }
                }
                //			if (scans.size() > 100) break;
            }

            cout << "Finished" << endl;
        }
    };

} // end of namespace

#endif
