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
//#define __STDCPP_MATH_SPEC_FUNCS__

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
#include <grd/TfUtils.h>
#include <grd/CorrespondenceGraph.h>
#include <grd/Profiler.h>
#include <grd/fileutils.h>
#include <grd/DistributionSignature.h>

struct ScanItem {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // General data
    int id;
    Eigen::Vector3d poseVec;
    grd::Transformation2d poseMat;
    grd::VectorPoint2d features;

    // Signatures
    grd::GlarotSignature glarot;
    grd::DistributionSignature grdBR; // GRD with biased Rayleigh
    grd::DistributionSignature grdEr; // GRD with Erlang
};
typedef std::vector<ScanItem, Eigen::aligned_allocator<ScanItem> > VectorScanItem;

enum SignatureType {
    GLAROT, GRD_BIASED_RAYLEIGH, GRD_ERLANG
};

struct LocalizationResult {
    std::string name; // method name (glarot, distrBiasray, distrErlang)
    int id1; // associated scan item id1
    int id2; // associated scan item id2
    int matchedNum; // number of associated keypoint features between id1 and id2
    double residual; // error residual
    double errPos; // position error of registration w.r.t. ground truth 
    double errAng; // angular error of registration w.r.t. ground truth
    bool valid;
};
typedef std::vector<LocalizationResult> VectorLocalizationResult;


Gnuplot gp("gnuplot -persist");

/**
 * Finds the closest matching items (scans) to item i-th and returns the distance 
 * and the indices of the candidates. 
 * @param i
 * @param items
 * @param k
 * @param candidates
 * @param indexDiffMin
 * @param online
 */
void findCandidates(int i, const VectorScanItem& items, int k, SignatureType signType, std::vector<std::pair<double, int> >& candidates, int indexDiffMin, bool online);
//void findCandidatesGlarot(int i, const VectorScanItem& items, int k, std::vector<std::pair<double, int> >& candidates, int indexDiffMin, bool online);

void associatePointToPoint(const VectorScanItem& items, int i1, int i2, double distCorrespMin, double distCorrespTol, LocalizationResult& res);

void associatePointToPoint2(const VectorScanItem& items, int i1, int i2, LocalizationResult& res);

double computeTransform(const grd::VectorPoint2d& points1, const grd::VectorPoint2d& points2, const std::vector<std::pair<int, int> >& indices, grd::Transformation2d& transform);

void computeNNAssociation(const grd::VectorPoint2d& points1, const grd::VectorPoint2d& points2, double distMax, const grd::Transformation2d& transform, std::vector<std::pair<int, int> >& indices);

void distAffine2d(const grd::Transformation2d& pose1, const grd::Transformation2d& pose2, double& distPosition, double& distAngle);

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
    //grd::VectorPoint2d features;
    double featureNumAvg;
    // GLAROT
    grd::GlarotSignature glarot;
    int glarotThetaNum, glarotRhoNum;
    double glarotRhoStep;
    int glarotCandidateNum;
    std::vector<std::pair<double, int> > candidateIds;
    VectorLocalizationResult glarotResults;
    // CG
    double cgDistMin;
    double cgDistTol;
    // GRD
    int grdThetaNum; // maximum order of Fourier series
    int grdRangeNum; // maximum order of Laguerre polynomial series
    double grdVMKappa;
    double grdBRSigma;
    double grdErSigma;
    double grdThetaTol;
    double grdCorrTol;
    int grdCandidateNum;
    //  General params
    ScanItem item;
    VectorScanItem items;
    bool enableGlarot;
    bool enableGrdBR;
    bool enableGrdEr;
    //int candidateInitNum;
    int indexDiffMin;
    bool onlineLoopclose;
    int evalBeg;
    int evalEnd;
    LocalizationResult res, resBest;
    std::ofstream resultFileGlarot;
    std::ofstream resultFileGrdBR;
    std::ofstream resultFileGrdEr;


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

    params.getParam<std::string>("in", filenameIn, std::string(""));
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.getParam<std::string>("out", filenameOut, std::string(""));

    // SKIP parameters
    params.getParam<double>("skipThres", skipThres, double(0.15));
    params.getParam<std::string>("skipScoreType", skipScoreType, std::string("triangle_diff"));
    params.getParam<double>("skipGapQ", skipGapQ, double(0.25));
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
    params.getParam<int>("glarotRhoNum", glarotRhoNum, int(80));
    params.getParam<double>("glarotRhoStep", glarotRhoStep, double(0.1));
    params.getParam<int>("glarotCandidateNum", glarotCandidateNum, int(10));

    // GRD parameters
    params.getParam<int>("grdThetaNum", grdThetaNum, int(15));
    params.getParam<int>("grdRangeNum", grdRangeNum, int(20));
    params.getParam<double>("grdVMKappa", grdVMKappa, double(20.0));
    params.getParam<double>("grdBRSigma", grdBRSigma, double(0.2));
    params.getParam<double>("grdErSigma", grdErSigma, double(0.2));
    params.getParam<double>("grdThetaTol", grdThetaTol, double(0.5));
    grdThetaTol *= M_PI / 180.0;
    params.getParam<double>("grdCorrTol", grdCorrTol, double(0.001));
    params.getParam<int>("grdCandidateNum", grdCandidateNum, int(10));

    // CG
    params.getParam<double>("cgDistMin", cgDistMin, double(0.20));
    params.getParam<double>("cgDistTol", cgDistTol, double(0.20));

    // General options
    params.getParam<bool>("plot", plotOn, bool(false));
    params.getParam<double>("rate", rate, double(20.0));
    params.getParam<bool>("enableGlarot", enableGlarot, false);
    params.getParam<bool>("enableGrdBR", enableGrdBR, false);
    params.getParam<bool>("enableGrdEr", enableGrdEr, false);
    //params.getParam<int>("candidateInitNum", candidateInitNum, int(10));
    params.getParam<int>("indexDiffMin", indexDiffMin, int(2));
    params.getParam<bool>("onlineLoopclose", onlineLoopclose, false);
    params.getParam<int>("evalBeg", evalBeg, int(0));
    params.getParam<int>("evalEnd", evalEnd, int(-1));

    std::cout << "\nParams:\n";
    params.write(std::cout);
    std::cout << std::endl;

    std::cout << "Reading from file \"" << filenameIn << "\"" << std::endl;
    reader.readCarmen(filenameIn);

    grd::Rate looprate(rate);
    featureNumAvg = 0.0;
    for (int i = 0; i < reader.scans.size(); ++i) {
        grd::ScopedTimer timerSignatures("all signatures");
        // Sets general information about current scan item
        item.id = i;
        item.poseVec = reader.scans[i].wTl;
        grd::TfUtils::vector3dToAffine2d(reader.scans[i].wTl, item.poseMat);

        // Computes the features of current scan item
        {
            grd::ScopedTimer timer("skip");
            item.features.clear();
            skip.extract(reader.scans[i].scan, item.features);
        }
        featureNumAvg += item.features.size();

        grd::Profiler::getProfiler().updateStat("featureNum", item.features.size());

        // Computes GLAROT for current scan item
        if (enableGlarot) {
            grd::ScopedTimer timer("glarot");
            item.glarot.setSize(glarotThetaNum, glarotRhoNum, glarotRhoStep);
            item.glarot.clear();
            item.glarot.setPoints(item.features);
        }

        // Computes GRD with biased Rayleigh
        if (enableGrdBR) {
            item.grdBR.setRangeDistribution(grd::DistributionSignature::BIASED_RAYLEIGH);
            item.grdBR.setSize(grdThetaNum, grdRangeNum);
            item.grdBR.setVonMisesKappa(grdVMKappa);
            item.grdBR.setBiasedRayleighSigma(grdBRSigma);
            item.grdBR.setCorrelationThetaTolerance(grdThetaTol);
            item.grdBR.setCorrelationFourierTolerance(grdCorrTol);
            grd::ScopedTimer timer("grd_br");
            item.grdBR.setPoints(item.features);
            //            std::cout << "  " << i << " getAutoCorrel() " << item.grdBR.getAutoCorrel() 
            //                    << " correlate(item.grdBR) " << item.grdBR.correlate(item.grdBR) << "\n";
        }

        // Computes GRD with biased Rayleigh
        if (enableGrdEr) {
            item.grdEr.setRangeDistribution(grd::DistributionSignature::BIASED_RAYLEIGH);
            item.grdEr.setSize(grdThetaNum, grdRangeNum);
            item.grdEr.setVonMisesKappa(grdVMKappa);
            item.grdEr.setBiasedRayleighSigma(grdBRSigma);
            item.grdEr.setCorrelationThetaTolerance(grdThetaTol);
            item.grdEr.setCorrelationFourierTolerance(grdCorrTol);
            grd::ScopedTimer timer("grd_er");
            item.grdEr.setPoints(item.features);
        }

        // Adds item to item vector
        items.push_back(item);

        std::cout << "scan " << i << ": found " << item.features.size() << " features" << std::endl;
        if (i % 30 == 0) {
            std::cout << "----\nProfiler stats (label mean stddev min max count):" << std::endl;
            grd::Profiler::getProfiler().printStats(std::cout);
            std::cout << "----\n";
        }


        if (plotOn) {
            std::stringstream ss;
            ss << "scan-" << std::setfill('0') << std::setw(5) << (i);
            grd::g_gnuplot.setWinParam(ss.str(), std::make_pair(0.0, 20.0), std::make_pair(-10.0, 20.0));
            grd::g_gnuplot.addPointVector("scan", reader.scans[i].scan.points, 13, 1, 2);
            grd::g_gnuplot.addPointVector("features", items[i].features, 7, 1, 1);
            grd::g_gnuplot.plot();
        }

        looprate.sleep();
    }

    featureNumAvg = featureNumAvg / items.size();
    std::cout << "\nLoaded " << items.size() << " scans (avg feature num " << featureNumAvg << ")\n" << std::endl;

    if (enableGlarot) {
        std::stringstream ss;
        ss << grd::getWithoutExtension(filenameIn) << "_glarot_";
        std::string filename = grd::generateStampedString(ss.str(), ".txt");
        std::cout << "opening result file \"" << filename << "\"" << std::endl;
        resultFileGlarot.open(filename);
        resultFileGlarot << "# Parameters:\n";
        params.write(resultFileGlarot, "#  ");
        resultFileGlarot << "#\n";
        resultFileGlarot << "# scanInputId method scanBestMatchId matchingNum residual errPos[m] errAng[rad]\n";
    }

    if (enableGrdBR) {
        std::stringstream ss;
        ss << grd::getWithoutExtension(filenameIn) << "_grd_br_";
        std::string filename = grd::generateStampedString(ss.str(), ".txt");
        std::cout << "opening result file \"" << filename << "\"" << std::endl;
        resultFileGrdBR.open(filename);
        resultFileGrdBR << "# Parameters:\n";
        params.write(resultFileGrdBR, "#  ");
        resultFileGrdBR << "#\n";
        resultFileGrdBR << "# scanInputId method scanBestMatchId matchingNum residual errPos[m] errAng[rad]\n";
    }

    if (enableGrdEr) {
        std::stringstream ss;
        ss << grd::getWithoutExtension(filenameIn) << "_grd_er_";
        std::string filename = grd::generateStampedString(ss.str(), ".txt");
        std::cout << "opening result file \"" << filename << "\"" << std::endl;
        resultFileGrdEr.open(filename);
        resultFileGrdEr << "# Parameters:\n";
        params.write(resultFileGrdEr, "#  ");
        resultFileGrdEr << "#\n";
        resultFileGrdEr << "# scanInputId method scanBestMatchId matchingNum residual errPos[m] errAng[rad]\n";
    }

    //    for (int i = 0; i < items.size(); ++i) {
    //        std::cout << "i " << i << ", items[i].grdBR.getAutoCorrel() " << items[i].grdBR.getAutoCorrel() 
    //                << ", items[i].grdBR.correlate(items[i].grdBR) " << items[i].grdBR.correlate(items[i].grdBR) << "\n";
    //    }

    // Visits and tests all the items with indices [evalBeg, evalEnd]. 
    // If evalEnd < 0, it means that we have to reach the end of the items vector.
    for (int i = evalBeg; i < items.size() && (evalEnd < 0 || i < evalEnd); ++i) {
        // Best GLAROT matches are the signatures that minimize distance
        if (enableGlarot && items[i].features.size() > 4) {
            // Finds the 10 best matching item with current scan item i
            findCandidates(i, items, glarotCandidateNum, GLAROT, candidateIds, indexDiffMin, onlineLoopclose);
            std::cout << "---\nscan " << i << "(" << items[i].features.size() << " features): " << candidateIds.size() << " candidates: ";
            for (auto& cand : candidateIds) {
                double distPos, distAng;
                distAffine2d(items[i].poseMat, items[cand.first].poseMat, distPos, distAng);
                std::cout << "(" << cand.second << ", " << cand.first << ": pos " << distPos << " ang[deg] " << (180.0 / M_PI * distAng) << ") ";
            }
            std::cout << std::endl;
            // The features/keypoints of i-th scan item are matched with the features of each scan item
            resBest.id1 = -1;
            resBest.id2 = -1;
            resBest.name = "glarot";
            resBest.matchedNum = 0;
            resBest.residual = resBest.errPos = resBest.errAng = 1e+6;
            for (auto& cand : candidateIds) {
                res.name = "glarot";
                associatePointToPoint(items, i, cand.second, cgDistMin, cgDistTol, res);
                //                glarotResults.push_back(res);
                //                std::cout << res.id1 << " (" << items[res.id1].features.size() << " points) to "
                //                        << res.id2 << " (" << items[res.id2].features.size() << " points): "
                //                        << " matched " << res.matchedNum << ", residual " << res.residual
                //                        << ", errPos " << res.errPos << ", errAng[deg] " << (180.0 / M_PI * res.errAng)
                //                        << "\n";
                //                std::cout << " " << res.id1 << " " << res.name << " " << res.id2 << " " << res.matchedNum << " "
                //                        << res.residual << " " << res.errPos << " " << res.errAng << "\n";
                if (res.matchedNum > resBest.matchedNum || (res.matchedNum == resBest.matchedNum && res.errPos < resBest.errPos)) {
                    resBest.id1 = res.id1;
                    resBest.id2 = res.id2;
                    resBest.matchedNum = res.matchedNum;
                    resBest.residual = res.residual;
                    resBest.errPos = res.errPos;
                    resBest.errAng = res.errAng;
                }
            }
            std::cout << "best: "
                    << resBest.id1 << " (" << items[resBest.id1].features.size() << " points) to "
                    << resBest.id2 << " (" << items[resBest.id2].features.size() << " points): "
                    << " matched " << resBest.matchedNum << ", residual " << resBest.residual
                    << ", errPos " << resBest.errPos << ", errAng[deg] " << (180.0 / M_PI * resBest.errAng)
                    << "\n";
            resultFileGlarot << resBest.id1 << " " << resBest.name << " " << resBest.id2 << " " << resBest.matchedNum << " "
                    << resBest.residual << " " << resBest.errPos << " " << resBest.errAng << "\n";
        }

        // Best GSR matches are the signatures that maximize correlation
        if (enableGrdBR && items[i].features.size() > 4) {
            // Finds the 10 best matching item with current scan item i
            findCandidates(i, items, grdCandidateNum, GRD_BIASED_RAYLEIGH, candidateIds, indexDiffMin, onlineLoopclose);
            std::cout << "---\nscan " << i << "(" << items[i].features.size() << " features): " << candidateIds.size() << " candidates: ";
            for (auto& cand : candidateIds) {
                double distPos, distAng;
                distAffine2d(items[i].poseMat, items[cand.first].poseMat, distPos, distAng);
                std::cout << "(" << cand.second << ", " << cand.first << ": pos " << distPos << " ang[deg] " << (180.0 / M_PI * distAng) << ") ";
            }
            std::cout << std::endl;
            // The features/keypoints of i-th scan item are matched with the features of each scan item
            resBest.id1 = -1;
            resBest.id2 = -1;
            resBest.name = "grd_br";
            resBest.matchedNum = 0;
            resBest.residual = resBest.errPos = resBest.errAng = 1e+6;
            for (auto& cand : candidateIds) {
                res.name = "grd_br";

                associatePointToPoint(items, i, cand.second, cgDistMin, cgDistTol, res);
                //                glarotResults.push_back(res);
                //                std::cout << res.id1 << " (" << items[res.id1].features.size() << " points) to "
                //                        << res.id2 << " (" << items[res.id2].features.size() << " points): "
                //                        << " matched " << res.matchedNum << ", residual " << res.residual
                //                        << ", errPos " << res.errPos << ", errAng[deg] " << (180.0 / M_PI * res.errAng)
                //                        << "\n";
                //                std::cout << " " << res.id1 << " " << res.name << " " << res.id2 << " " << res.matchedNum << " "
                //                        << res.residual << " " << res.errPos << " " << res.errAng << "\n";
                //if (res.matchedNum > resBest.matchedNum || (res.matchedNum == resBest.matchedNum && res.errPos < resBest.errPos)) {
                if (res.errPos < resBest.errPos) {
                    resBest.id1 = res.id1;
                    resBest.id2 = res.id2;
                    resBest.matchedNum = res.matchedNum;
                    resBest.residual = res.residual;
                    resBest.errPos = res.errPos;
                    resBest.errAng = res.errAng;
                }
            }
            std::cout << "best: "
                    << resBest.id1 << " (" << items[resBest.id1].features.size() << " points) to "
                    << resBest.id2 << " (" << items[resBest.id2].features.size() << " points): "
                    << " matched " << resBest.matchedNum << ", residual " << resBest.residual
                    << ", errPos " << resBest.errPos << ", errAng[deg] " << (180.0 / M_PI * resBest.errAng)
                    << "\n";
            resultFileGrdBR << resBest.id1 << " " << resBest.name << " " << resBest.id2 << " " << resBest.matchedNum << " "
                    << resBest.residual << " " << resBest.errPos << " " << resBest.errAng << "\n";
        }

        // Best GSR matches are the signatures that maximize correlation
        if (enableGrdEr && items[i].features.size() > 4) {
            // Finds the 10 best matching item with current scan item i
            findCandidates(i, items, grdCandidateNum, GRD_ERLANG, candidateIds, indexDiffMin, onlineLoopclose);
            std::cout << "---\nscan " << i << "(" << items[i].features.size() << " features): " << candidateIds.size() << " candidates: ";
            for (auto& cand : candidateIds) {
                double distPos, distAng;
                distAffine2d(items[i].poseMat, items[cand.first].poseMat, distPos, distAng);
                std::cout << "(" << cand.second << ", " << cand.first << ": pos " << distPos << " ang[deg] " << (180.0 / M_PI * distAng) << ") ";
            }
            std::cout << std::endl;
            // The features/keypoints of i-th scan item are matched with the features of each scan item
            resBest.id1 = -1;
            resBest.id2 = -1;
            resBest.name = "grd_er";
            resBest.matchedNum = 0;
            resBest.residual = resBest.errPos = resBest.errAng = 1e+6;
            for (auto& cand : candidateIds) {
                res.name = "grd_er";
                associatePointToPoint(items, i, cand.second, cgDistMin, cgDistTol, res);
                //                glarotResults.push_back(res);
                //                std::cout << res.id1 << " (" << items[res.id1].features.size() << " points) to "
                //                        << res.id2 << " (" << items[res.id2].features.size() << " points): "
                //                        << " matched " << res.matchedNum << ", residual " << res.residual
                //                        << ", errPos " << res.errPos << ", errAng[deg] " << (180.0 / M_PI * res.errAng)
                //                        << "\n";
                //                std::cout << " " << res.id1 << " " << res.name << " " << res.id2 << " " << res.matchedNum << " "
                //                        << res.residual << " " << res.errPos << " " << res.errAng << "\n";
                if (res.matchedNum > resBest.matchedNum || (res.matchedNum == resBest.matchedNum && res.errPos < resBest.errPos)) {
                    resBest.id1 = res.id1;
                    resBest.id2 = res.id2;
                    resBest.matchedNum = res.matchedNum;
                    resBest.residual = res.residual;
                    resBest.errPos = res.errPos;
                    resBest.errAng = res.errAng;
                }
            }
            std::cout << "best: "
                    << resBest.id1 << " (" << items[resBest.id1].features.size() << " points) to "
                    << resBest.id2 << " (" << items[resBest.id2].features.size() << " points): "
                    << " matched " << resBest.matchedNum << ", residual " << resBest.residual
                    << ", errPos " << resBest.errPos << ", errAng[deg] " << (180.0 / M_PI * resBest.errAng)
                    << "\n";
            resultFileGrdEr << resBest.id1 << " " << resBest.name << " " << resBest.id2 << " " << resBest.matchedNum << " "
                    << resBest.residual << " " << resBest.errPos << " " << resBest.errAng << "\n";
        }
    }

    if (enableGlarot) {
        resultFileGlarot.close();
    }
    if (enableGrdBR) {
        resultFileGrdBR.close();
    }
    if (enableGrdEr) {
        resultFileGrdEr.close();
    }

    std::cout << "----\nProfiler stats (label mean stddev min max count):" << std::endl;
    grd::Profiler::getProfiler().printStats(std::cout);
    std::cout << "----\n";

    return 0;
}

void findCandidates(int i, const VectorScanItem& items, int k, SignatureType signType, std::vector<std::pair<double, int> >& candidates, int indexDiffMin, bool online) {
    // Pushes the distances between scans[id].glare and other glare
    typedef std::pair<double, int> Type;
    typedef std::vector<Type> Container;
    typedef std::less<Type> Comparator;
    typedef std::priority_queue<Type, Container, std::less<Type> > QueueMin;
    typedef std::priority_queue<Type, Container, std::greater<Type> > QueueMax;
    //std::priority_queue<Type, Container, Comparator> q;
    QueueMin queueGlarot;
    QueueMax queueGrd;
    double dist, corr;

    grd::ScopedTimer timerFind("findCandidates");

    //std::cout << "i " << i << ", items[i].grdBR.correlate(items[i].grdBR) " << items[i].grdBR.correlate(items[i].grdBR) << "\n";

    // Finds potentially matching points 
    for (int j = 0; j < items.size(); ++j) {
        // In online learning, it can only check the poses in the past, i.e. j < i. 
        // We also requires that the items to be compares are not too close in the sequence:
        // loop closure is for matching not temporally sequential 
        if (std::abs(j - i) > indexDiffMin && (!online || j < i)) {
            if (signType == GLAROT) {
                grd::ScopedTimer timerFind("glarot.distanceL1Min");
                dist = items[j].glarot.distanceL1Min(items[i].glarot);
                if (queueGlarot.size() < k) {
                    queueGlarot.push(std::make_pair(dist, j));
                } else if (dist < queueGlarot.top().first) {
                    queueGlarot.pop();
                    queueGlarot.push(std::make_pair(dist, j));
                }
            } else if (signType == GRD_BIASED_RAYLEIGH || signType == GRD_ERLANG) {
                if (signType == GRD_BIASED_RAYLEIGH) {
                    grd::ScopedTimer timerFind("grdBR.correlate");
                    corr = items[i].grdBR.correlate(items[j].grdBR);
                } else if (signType == GRD_ERLANG) {
                    grd::ScopedTimer timerFind("grdEr.correlate");
                    corr = items[i].grdEr.correlate(items[j].grdEr);
                }
                if (queueGrd.size() < k) {
                    queueGrd.push(std::make_pair(corr, j));
                } else if (corr > queueGrd.top().first) {
                    queueGrd.pop();
                    queueGrd.push(std::make_pair(corr, j));
                }
                //                std::cout << "  checking " << i << " to " << j << " std::abs(j - i) " << std::abs(j - i) << ": corr " << corr << " ";
                //                if (!queueGrd.empty())
                //                    std::cout << " (" << queueGrd.size() << "-best: " << queueGrd.top().second << " corr " << queueGrd.top().first << ")\n";
                //                else
                //                    std::cout << "\n";
            }
        }
    }

    // Copies the k-best minimum values into candidates
    candidates.clear();
    candidates.resize(k);
    for (int m = k - 1; m >= 0; --m) {
        if (signType == GLAROT && !queueGlarot.empty()) {
            candidates[m].first = queueGlarot.top().first;
            candidates[m].second = queueGlarot.top().second;
            queueGlarot.pop();
        } else if ((signType == GRD_BIASED_RAYLEIGH || signType == GRD_ERLANG) && !queueGrd.empty()) {
            candidates[m].first = queueGrd.top().first;
            candidates[m].second = queueGrd.top().second;
            queueGrd.pop();
        }
        //std::cout << " " << m << ") scan " << i << " matching with " << candidates[m].second << ": dist/corr " << candidates[m].first << "\n";
    }
}

void associatePointToPoint(const VectorScanItem& items, int i1, int i2, double distCorrespMin, double distCorrespTol, LocalizationResult& res) {
    grd::Transformation2d inputTtargetTrue, targetTinputTrue;
    grd::Transformation2d targetTinputCG, inputTtargetCG;
    //grd::Transformation2d targetTinputRefined, inputTtargetRefined;
    std::vector<std::pair<int, int> > assocCG;
    std::vector<std::pair<int, int> > assocNN;
    //    double distCorrespTol = 0.20;
    //    double distCorrespMin = 0.05;
    //double dist_nn_assoc_max = 1.0;

    // The groundtruth relative poses bewteen the reference frames of the two 
    // candicate matching items
    inputTtargetTrue = items[i1].poseMat.inverse() * items[i2].poseMat;
    targetTinputTrue = inputTtargetTrue.inverse();

    grd::ScopedTimer timerFind("associatePointToPoint");

    // Associates the features between the items and estimates the transformation
    grd::CorrespondenceGraph corrGraph;
    corrGraph.setDistanceToll(distCorrespTol);
    corrGraph.setDistanceMin(distCorrespMin);
    corrGraph.setInputs(items[i1].features);
    corrGraph.setTargets(items[i2].features);
    corrGraph.computeTransform(targetTinputCG, assocCG);
    inputTtargetCG = targetTinputCG.inverse();

    //    std::cout << "   CG associated " << i1 << " (" << items[i1].features.size() << " features)"
    //            << " vs " << i2 << " (" << items[i2].features.size() << " features): " << assocCorrGraph.size() << std::endl;

    // Refines the transformation by matching all the points... OMITTED!
    //    computeNNAssociation(items[i1].features, items[i2].features, dist_nn_assoc_max, inputTtargetCG, assocNN);
    //    std::cout << "   NN associated " << i1 << " vs " << i2 << ": " << assocNN.size() << std::endl;
    //
    //    res.residual = computeTransform(items[i1].features, items[i2].features, assocNN, inputTtargetRefined);
    //    targetTinputRefined = inputTtargetRefined.inverse();
    res.id1 = i1;
    res.id2 = i2;
    res.matchedNum = assocCG.size();

    distAffine2d(inputTtargetCG, inputTtargetTrue, res.errPos, res.errAng);
    res.valid = true;
}

void associatePointToPoint2(const VectorScanItem& items, int i1, int i2, LocalizationResult& res) {
    grd::Transformation2d inputTtargetTrue, targetTinputTrue;
    //grd::Transformation2d targetTinputCorr, inputTtargetCorr;
    grd::Transformation2d targetTinputRefined, inputTtargetRefined;
    //std::vector<std::pair<int, int> > assocCorrGraph;
    std::vector<std::pair<int, int> > assocNN;
    //double distCorrespTol = 0.08;
    //double distCorrespMin = 0.20;
    double dist_nn_assoc_max = 1.0;

    // The groundtruth relative poses bewteen the reference frames of the two 
    // candicate matching items
    inputTtargetTrue = items[i2].poseMat.inverse() * items[i1].poseMat;
    targetTinputTrue = inputTtargetTrue.inverse();

    // Refines the transformation by matching all the points... OMITTED!
    computeNNAssociation(items[i2].features, items[i1].features, dist_nn_assoc_max, inputTtargetTrue, assocNN);
    //    std::cout << "   NN " << i1 << " (" << items[i1].features.size() << " features)"
    //            << " vs " << i2 << " (" << items[i2].features.size() << " features): associated " << assocNN.size() << " features" << std::endl;

    res.residual = computeTransform(items[i1].features, items[i2].features, assocNN, inputTtargetRefined);
    res.id1 = i1;
    res.id2 = i2;
    res.matchedNum = assocNN.size();
    targetTinputRefined = inputTtargetRefined.inverse();
    distAffine2d(items[i1].poseMat, items[i2].poseMat * inputTtargetRefined, res.errPos, res.errAng);
    res.valid = true;

    //    gp << "set term wxt 2\n";
    //    gp << "plot "
    //            << "  '-' w p pt 7 ps 1 lc rgb \"red\" notitle"
    //            << ", '-' w p pt 7 ps 1 lc rgb \"green\" notitle"
    //            << ", '-' w l notitle"
    //            << std::endl;
    //    for (auto& f : items[i1].features) {
    //        Point2d ft = items[i1].poseMat * f;
    //        gp << ft.x() << " " << ft.y() << "\n";
    //    }
    //    gp << "e\n";
    //    for (auto& f : items[i2].features) {
    //        Point2d ft = items[i2].poseMat * f;
    //        gp << ft.x() << " " << ft.y() << "\n";
    //    }
    //    gp << "e\n";
    //    for (auto& a : assocNN) {
    //        Point2d f1 = items[i1].poseMat * items[i1].features[a.first];
    //        Point2d f2 = items[i2].poseMat * items[i2].features[a.second];
    //        gp << f1.x() << " " << f1.y() << "\n"
    //           << f2.x() << " " << f2.y() << "\n\n";
    //    }
    //    gp << "e\n";
}

double computeTransform(const grd::VectorPoint2d& points1, const grd::VectorPoint2d& points2,
        const std::vector<std::pair<int, int> >& indices, grd::Transformation2d& transform) {
    grd::Point2d t1 = Eigen::Vector2d::Zero();
    grd::Point2d t2 = Eigen::Vector2d::Zero();
    Eigen::Matrix2d S = Eigen::Matrix2d::Zero();
    int n = 0;
    for (int i = 0; i < (int) indices.size(); ++i) {
        if (0 <= indices[i].first && indices[i].first < (int) points1.size() &&
                0 <= indices[i].second && indices[i].second < (int) points2.size()) {
            t1 += points1[indices[i].first];
            t2 += points2[indices[i].second];
            n++;
        }
    }
    if (n == 0) {
        return 1e+6;
    }
    t1 = (1.0 / n) * t1;
    t2 = (1.0 / n) * t2;
    for (int i = 0; i < (int) indices.size(); ++i) {
        if (0 <= indices[i].first && indices[i].first < (int) points1.size() &&
                0 <= indices[i].second && indices[i].second < (int) points2.size()) {
            S += (points2[indices[i].second] - t2) * (points1[indices[i].first] - t1).transpose();
        }
    }
    double theta = atan2(S(0, 1) - S(1, 0), S(0, 0) + S(1, 1));

    Eigen::Rotation2Dd rot(theta);
    Eigen::Vector2d transl = t1 - (rot * t2);
    transform = grd::Transformation2d::Identity();
    transform.prerotate(rot);
    transform.pretranslate(transl);

    double residual = 0.0;
    for (int i = 0; i < (int) indices.size(); ++i) {
        if (0 <= indices[i].first && indices[i].first < (int) points1.size() &&
                0 <= indices[i].second && indices[i].second < (int) points2.size()) {
            residual += (points1[indices[i].first] - transform * points2[indices[i].second]).norm();
        }
    }
    residual = residual / n;
    return residual;
}

void computeNNAssociation(const grd::VectorPoint2d& points1, const grd::VectorPoint2d& points2, double distMax, const grd::Transformation2d& frame1Tframe2, std::vector<std::pair<int, int> >& indices) {
    indices.clear();
    int i2min;
    double dmin, d;
    double davg = 0.0;
    int counter = 0;
    for (int i1 = 0; i1 < points1.size(); ++i1) {
        i2min = -1;
        dmin = 1.05 * distMax;
        for (int i2 = 0; i2 < points2.size(); ++i2) {
            Eigen::Vector2d p2trans = frame1Tframe2 * points2[i2];
            d = (points1[i1] - p2trans).norm();
            //std::cout << "i1 " << i1 << ": " << points1[i1].transpose() << "; i2 " << i2 << ": " << p2trans.transpose() << "; d " << d << std::endl;
            if (d < dmin) {
                i2min = i2;
                dmin = d;
            }
        }
        if (i2min >= 0 && dmin < distMax) {
            indices.push_back(std::make_pair(i1, i2min));
        }
    }
}

void distAffine2d(const grd::Transformation2d& pose1, const grd::Transformation2d& pose2, double& distPosition, double& distAngle) {
    Eigen::Vector3d vec1 = grd::TfUtils::affine2dToVector3d(pose1);
    Eigen::Vector3d vec2 = grd::TfUtils::affine2dToVector3d(pose2);
    distPosition = (vec1.block<2, 1>(0, 0) - vec2.block<2, 1>(0, 0)).norm();
    distAngle = std::abs(vec1(2) - vec2(2));
    if (distAngle > M_PI) {
        distAngle = 2 * M_PI - distAngle;
    }
}

