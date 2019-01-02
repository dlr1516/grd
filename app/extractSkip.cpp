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
#define __STDCPP_MATH_SPEC_FUNCS__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <grd/Point.h>
#include <grd/LaserScan.h>
#include <grd/ParamMap.h>
#include <grd/SkipExtractor.h>
#include <grd/CarmenV2Reader.h>
#include <grd/GnuplotVisualizer.h>
#include <grd/Rate.h>
#include <grd/GlarotSignature.h>


int main(int argc, char** argv) {
    // General application parameters
    grd::ParamMap params;
    std::string filenameIn;
    std::string filenameCfg;
    std::string filenameOut;
    grd::CarmenV2Reader reader;
    bool plotOn;
    double rate;
    // Skip feature extraction
    grd::SkipExtractor skip; 
    double skipThres;
    std::string skipScoreType;
    double skipGapQ, skipGapM;
    grd::VectorPoint2d features;
    // GLAROT
    grd::GlarotSignature glarot;
    int glarotThetaNum, glarotRhoNum;
    double glarotRhoStep;
 
    
    // Reads file params from command line
    params.read(argc, argv);
    params.getParam<std::string>("in", filenameIn, std::string(""));
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.getParam<std::string>("out", filenameOut, std::string(""));

    // Reads potential parameter value from file and, then again, from command line
    if (!params.read(filenameCfg)) {
        std::cout << "Cannot open configuration file \"" << filenameCfg << "\": using default values" << std::endl;
    }
    params.read(argc, argv);
    
    // SKIP parameters
    params.getParam<double>("skipThres", skipThres, double(0.15));
    params.getParam<std::string>("skipScoreType", skipScoreType, std::string("triangle_diff"));    
    params.getParam<double>("skipGapQ", skipGapQ, double(0.2));
    params.getParam<double>("skipGapM", skipGapM, double(0.1));
    if (skipScoreType == "triangle_diff") {
        skip.setScoreType(grd::SkipExtractor::SCORE_TRIANGLE_DIFF, skipThres);
    } else if (skipScoreType == "dce") {
        skip.setScoreType(grd::SkipExtractor::SCORE_DCE, skipThres);
    } else {
        skip.setScoreType(grd::SkipExtractor::SCORE_TRIANGLE_DIFF, skipThres);
    }
    skip.setGapThres(skipGapQ, skipGapM);
    
    // GLAROT parameters
    params.getParam<int>("glarotThetaNum", glarotThetaNum, int(8));
    params.getParam<int>("glarotRhoNum", glarotRhoNum, int(50));
    params.getParam<double>("glarotRhoStep", glarotRhoStep, double(0.2));
    
    // General options
    params.getParam<bool>("plot", plotOn, bool(true));
    params.getParam<double>("rate", rate, double(10.0));
    
    std::cout << "\nParams:\n";
    params.write(std::cout);
    std::cout << std::endl;
    
    std::cout << "Reading from file \"" << filenameIn << "\"" << std::endl;
    reader.readCarmen(filenameIn);
    
    grd::Rate looprate(rate);
    for (int i = 0; i < reader.scans.size(); ++i) {
        skip.extract(reader.scans[i].scan, features);
        std::cout << "scan " << i << ": found " << features.size() << " features" << std::endl; 
        
        glarot.clear(); 
        glarot.setPoints(reader.scans[i].scan.points);
        
        if (plotOn) {
            std::stringstream ss;
            ss << "scan-" << std::setfill('0') << std::setw(5) << (i);
            grd::g_gnuplot.setWinParam(ss.str(), std::make_pair(0.0, 20.0), std::make_pair(-10.0, 20.0));
            grd::g_gnuplot.addPointVector("scan", reader.scans[i].scan.points, 13, 1, 2);
            grd::g_gnuplot.addPointVector("features", features, 7, 1, 1);
            grd::g_gnuplot.plot();
        }
        
        looprate.sleep();
    }
    
    //std::cout << "laguerre(3,2.7) " << std::laguerre(3, 2.7) << std::endl;
    return 0;
}

