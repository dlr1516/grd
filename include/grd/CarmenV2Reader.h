#pragma once


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
