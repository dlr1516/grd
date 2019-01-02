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
#include <grd/functions.h>

using namespace grd;

struct ScanItem {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // General data
    int id;
    Eigen::Vector3d poseVec;
    Transformation2d poseMat;
    VectorPoint2d features;

    // Signatures
    GlarotSignature glarot;
};
typedef std::vector<ScanItem, Eigen::aligned_allocator<ScanItem> > VectorScanItem;


Gnuplot gp("gnuplot -persist");
//std::ostream& gp = std::cout;

///**
// * Finds the closest matching items (scans) to item i-th and returns the distance 
// * and the indices of the candidates. 
// * @param i
// * @param items
// * @param k
// * @param candidates
// * @param indexDiffMin
// * @param online
// */
//void findCandidatesGlarot(int i, const VectorScanItem& items, int k, std::vector<std::pair<double, int> >& candidates, int indexDiffMin, bool online);
//
//void associatePointToPoint(const VectorScanItem& items, int i1, int i2, LocalizationResult& res);
//
//void associatePointToPoint2(const VectorScanItem& items, int i1, int i2, LocalizationResult& res);
//
//double computeTransform(const VectorPoint2d& points1, const VectorPoint2d& points2, const std::vector<std::pair<int, int> >& indices, Transformation2d& transform);
//
//void computeNNAssociation(const VectorPoint2d& points1, const VectorPoint2d& points2, double distMax, const Transformation2d& transform, std::vector<std::pair<int, int> >& indices);
//
//void distAffine2d(const Transformation2d& pose1, const Transformation2d& pose2, double& distPosition, double& distAngle);

void plotPoints(std::ostream& out, const VectorPoint2d& points);

void plotCompatible(std::ostream& out, const VectorPoint2d& points, const std::vector<std::pair<int, int> >& associations, bool useFirst);

void plotAssociations(std::ostream& out, const VectorPoint2d& points1, const VectorPoint2d& points2, const std::vector<std::pair<int, int> >& associations);

int main(int argc, char** argv) {
    // General application parameters
    ParamMap params;
    std::string filenameIn;
    std::string filenameCfg;
    std::string filenameOut;
    CarmenV2Reader reader;
    bool plotOn;
    bool savePlot;
    double rate;
    // Skip feature extraction
    SkipExtractor skip;
    double skipThres;
    std::string skipScoreType;
    double skipGapQ, skipGapM;

    // GLAROT
    GlarotSignature glarot;
    int glarotThetaNum, glarotRhoNum;
    double glarotRhoStep;
    int glarotCandidateNum;
    bool enableGlarot;

    // CG
    double cgDistMin;
    double cgDistTol;
    Transformation2d srcTdstCG;
    std::vector<std::pair<int, int> > assocCG;
    VectorPoint2d srcFeaturesTransf;

    // GRD
    DistributionSignature grdSrc;
    DistributionSignature grdDst;
    std::vector<double> corrCoeffs;
    bool enableGRD;
    int grdThetaNum;
    int grdRangeNum;
    double grdKappa; // von Mises Kappa
    double grdSigma; // bias Rayleigh width
    std::vector<double> plotCorrelGrd;

    // 
    ScanItem itemSrc;
    ScanItem itemDst;
    CorrespondenceGraph corrGraph;
    int srcId;
    int dstId;
    double plotXMin = 1e+6;
    double plotXMax =-1e+6;
    double plotYMin = 1e+6;
    double plotYMax =-1e+6;


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
    params.getParam<double>("skipGapQ", skipGapQ, double(0.25));
    params.getParam<double>("skipGapM", skipGapM, double(0.1));
    if (skipScoreType == "triangle_diff") {
        skip.setScoreType(SkipExtractor::SCORE_TRIANGLE_DIFF, skipThres);
    } else if (skipScoreType == "dce") {
        skip.setScoreType(SkipExtractor::SCORE_DCE, skipThres);
    } else {
        skip.setScoreType(SkipExtractor::SCORE_TRIANGLE_DIFF, skipThres);
    }
    skip.setGapThres(skipGapQ, skipGapM);

    // GLAROT parameters
    params.getParam<int>("glarotThetaNum", glarotThetaNum, int(8));
    params.getParam<int>("glarotRhoNum", glarotRhoNum, int(80));
    params.getParam<double>("glarotRhoStep", glarotRhoStep, double(0.1));
    params.getParam<int>("glarotCandidateNum", glarotCandidateNum, int(10));

    // CG
    params.getParam<double>("cgDistMin", cgDistMin, double(0.02));
    params.getParam<double>("cgDistTol", cgDistTol, double(0.20));
    corrGraph.setDistanceMin(cgDistMin);
    corrGraph.setDistanceToll(cgDistTol);

    // GRD
    params.getParam<bool>("enableGRD", enableGRD, bool(false));
    params.getParam<int>("grdThetaNum", grdThetaNum, int(5));
    params.getParam<int>("grdRangeNum", grdRangeNum, int(8));
    params.getParam<double>("grdKappa", grdKappa, double(1.00));
    params.getParam<double>("grdSigma", grdSigma, double(0.30));


    // General options
    params.getParam<bool>("plotOn", plotOn, bool(false));
    params.getParam<bool>("savePlot", savePlot, bool(false));
    params.getParam<double>("rate", rate, double(20.0));
    params.getParam<bool>("enableGlarot", enableGlarot, true);
    params.getParam<int>("srcId", srcId, int(-1));
    params.getParam<int>("dstId", dstId, int(-1));

    std::cout << "\nParams:\n";
    params.write(std::cout);
    std::cout << std::endl;

    std::cout << "Reading from file \"" << filenameIn << "\"" << std::endl;
    reader.readCarmen(filenameIn);

    std::cout << "compare scans " << srcId << " and " << dstId << std::endl;
    if (srcId < 0 || srcId >= reader.scans.size()) {
        std::cerr << "invalid srcId " << srcId << " not in [0, " << reader.scans.size() << std::endl;
        return -1;
    }
    if (dstId < 0 || dstId >= reader.scans.size()) {
        std::cerr << "invalid srcId " << dstId << " not in [0, " << reader.scans.size() << std::endl;
        return -1;
    }

    itemSrc.id = srcId;
    itemSrc.poseVec = reader.scans[srcId].wTl;
    TfUtils::vector3dToAffine2d(reader.scans[srcId].wTl, itemSrc.poseMat);
    itemSrc.features.clear();
    skip.extract(reader.scans[srcId].scan, itemSrc.features);
    itemSrc.glarot.setSize(glarotThetaNum, glarotRhoNum, glarotRhoStep);
    itemSrc.glarot.clear();
    itemSrc.glarot.setPoints(itemSrc.features);
    std::cout << "computing " << itemSrc.features.size() << " features and glarot for " << srcId << std::endl;


    itemDst.id = dstId;
    itemDst.poseVec = reader.scans[dstId].wTl;
    TfUtils::vector3dToAffine2d(reader.scans[dstId].wTl, itemDst.poseMat);
    itemDst.features.clear();
    skip.extract(reader.scans[dstId].scan, itemDst.features);
    itemDst.glarot.setSize(glarotThetaNum, glarotRhoNum, glarotRhoStep);
    itemDst.glarot.clear();
    itemDst.glarot.setPoints(itemSrc.features);
    std::cout << "computing " << itemDst.features.size() << " features and glarot for " << dstId << std::endl;

    if (plotOn) {
        std::stringstream ss;
        ss << "scan-" << std::setfill('0') << std::setw(5) << srcId;
        g_gnuplot.setTerminalWxt(0);
        g_gnuplot.setWinParam(ss.str(), std::make_pair(0.0, 30.0), std::make_pair(-10.0, 10.0));
        g_gnuplot.addPointVector("scan src", reader.scans[srcId].scan.points, 13, 1, 2);
        g_gnuplot.addPointVector("features src", itemSrc.features, 7, 1, 1);
        g_gnuplot.plot();

        ss.str("");
        ss << "scan-" << std::setfill('0') << std::setw(5) << dstId;
        g_gnuplot.setTerminalWxt(1);
        g_gnuplot.setWinParam(ss.str(), std::make_pair(0.0, 30.0), std::make_pair(-10.0, 10.0));
        g_gnuplot.addPointVector("scan dst", reader.scans[dstId].scan.points, 13, 1, 2);
        g_gnuplot.addPointVector("features dst", itemDst.features, 7, 1, 1);
        g_gnuplot.plot();
    }
    
    for (auto& p : itemSrc.features) {
        if (p.x() < plotXMin) plotXMin = p.x();
        if (p.x() > plotXMax) plotXMax = p.x();
        if (p.y() < plotYMin) plotYMin = p.y();
        if (p.y() > plotYMax) plotYMax = p.y();
    }

    corrGraph.setInputs(itemSrc.features);
    corrGraph.setTargets(itemDst.features);
    //corrGraph.associate(associations);
    corrGraph.computeTransform(srcTdstCG, assocCG);

    std::cout << "\nAssociations\n";
    for (auto& a : assocCG) {
        std::cout << "  " << a.first << " " << a.second << std::endl;
    }

    for (auto& f : itemSrc.features) {
        srcFeaturesTransf.push_back(srcTdstCG * f);
    }

    if (plotOn) {
        std::cout << " xmin " << plotXMin << ", xmax " << plotXMax << ", ymin " << plotYMin << ", ymax " << plotYMax << std::endl;
        // Plots
        //Gnuplot gp("gnuplot -persist");
        //  gp << "set terminal pngcairo size 550,500 enhanced font 'Verdana,10'\n"
        //     << "set output \"correspondence_graph.png\"\n";
        //gp << "unset border; unset xtics; unset ytics; unset x2tics; unset y2tics;\n";
        if (savePlot) {
            gp << "set terminal postscript eps enhanced color font 'Arial,24'\n"
               << "set output \"cg_" << srcId << "_" << dstId << ".eps\"\n";
        }
        //gp << "set key font \",12\"\n"
        gp << "set key right bottom\n";
        //gp << "set key vertical maxrows 2\n";
        gp << "set title \"correspondence graph\"\n";
//        gp << "set xrange[" << plotXMin << ":" << (plotXMax + 0.2 * (plotXMax-plotXMin)) << "]\n";
//        gp << "set yrange[" << (plotYMin - 0.2 * (plotYMax-plotYMin)) << ":" << plotYMax << "]\n";
        gp << "set size ratio -1\n";
        gp << "plot ";
        gp << "'-' notitle w l lt 0 lc 3 lw 0.6";
        gp << ", '-' notitle w l lt 1 lc 1 lw 0.4";
        gp << ", '-' notitle w l lt 2 lc 2 lw 0.4";
        //gp << ", '-' title \"cg scan " << srcId << "\" w p pt 7 lc 1 ps 0.9 ";
        //gp << ", '-' title \"cg scan " << dstId << "\" w p pt 7 lc 2 ps 0.9 ";
        //gp << ", '-' title \"scan " << srcId << " aligned\" w p pt 7 lc 3 ps 0.6 ";
        gp << ", '-' title \"cg features src\" w p pt 7 lc 1 ps 0.9 ";
        gp << ", '-' title \"cg features dst\" w p pt 7 lc 2 ps 0.9 ";
        gp << ", '-' title \"features src aligned\" w p pt 7 lc 3 ps 0.6 ";
        gp << std::endl;
        plotAssociations(gp, itemSrc.features, itemDst.features, assocCG);
        gp << "e\n";
        plotCompatible(gp, itemSrc.features, assocCG, true);
        gp << "e\n";
        plotCompatible(gp, itemDst.features, assocCG, false);
        gp << "e\n";
        plotPoints(gp, itemSrc.features);
        gp << "e\n";
        plotPoints(gp, itemDst.features);
        gp << "e\n";
        plotPoints(gp, srcFeaturesTransf);
        gp << "e\n";
    }

    if (enableGRD) {
        std::cout << "compute signature GRD for src\n";
        grdSrc.setVonMisesKappa(grdKappa);
        grdSrc.setBiasedRayleighSigma(grdSigma);
        grdSrc.setSize(grdThetaNum, grdRangeNum);
        grdSrc.setPoints(itemSrc.features);
        std::cout << "src coeffs:\n" << grdSrc.getSignatureCoeffs() << "\n"
                << "  auto-correlation: " << grdSrc.getAutoCorrel() << "\n";

        std::cout << "compute signature GRD for dst\n";
        grdDst.setVonMisesKappa(grdKappa);
        grdDst.setBiasedRayleighSigma(grdSigma);
        grdDst.setSize(grdThetaNum, grdRangeNum);
        grdDst.setPoints(itemDst.features);
        std::cout << "dst coeffs:\n" << grdDst.getSignatureCoeffs() << "\n"
                << "  auto-correlation: " << grdDst.getAutoCorrel() << "\n";

        std::cout << "correlate signatures GRD src/dst\n";
        double thetaShift;
        double corrMax = grdSrc.correlate(grdDst, thetaShift, corrCoeffs);
        std::cout << "  correlation max " << corrMax << " at angle[deg]: " << (180.0 / M_PI * thetaShift) << std::endl;

        std::cout << "plotOn " << plotOn << std::endl;
        if (plotOn) {
            std::cout << "plotting corrlation..." << std::endl;
            gp << "set term wxt 3\n";
            gp << "set title \"GRD correlation scans " << srcId << " " << dstId << "\"\n";
            gp << "set xlabel \"phi [deg]\"\n";
            gp << "set ylabel \"correlation\"\n";
            gp << "plot ";
            gp << "'-' title \"correl\" w l lt 1 lc 3 lw 1.0";
            gp << ", '-' title \"correl max\" w p pt 7 lc 1 ps 1.3";
            gp << std::endl;
            for (int t = 0; t < 360; ++t) {
                double corr = evaluateFourier(corrCoeffs, M_PI / 180.0 * t);
                gp << t << " " << corr << "\n";
            }
            gp << "e\n";
            gp << (180.0 / M_PI * thetaShift) << " " << corrMax << "\n";
            gp << "e\n";
        }
    }

    //    Rate looprate(rate);
    //    featureNumAvg = 0.0;
    //    for (int i = 0; i < reader.scans.size(); ++i) {
    //        ScopedTimer timerSignatures("all signatures");
    //        // Sets general information about current scan item
    //        item.id = i;
    //        item.poseVec = reader.scans[i].wTl;
    //        TfUtils::vector3dToAffine2d(reader.scans[i].wTl, item.poseMat);
    //
    //        // Computes the features of current scan item
    //        {
    //            ScopedTimer timer("skip");
    //            item.features.clear();
    //            skip.extract(reader.scans[i].scan, item.features);
    //            std::cout << "scan " << i << ": found " << item.features.size() << " features" << std::endl;
    //        }
    //        featureNumAvg += item.features.size();
    //
    //        // Computes GLAROT for current scan item
    //        {
    //            ScopedTimer timer("glarot");
    //            item.glarot.setSize(glarotThetaNum, glarotRhoNum, glarotRhoStep);
    //            item.glarot.clear();
    //            item.glarot.setPoints(reader.scans[i].scan.points);
    //        }
    //        if (i % 10 == 0) {
    //            std::cout << "----\nProfiler stats (label mean stddev min max count):" << std::endl;
    //            Profiler::getProfiler().printStats(std::cout);
    //            std::cout << "----\n";
    //        }
    //
    //        // Adds item to item vector
    //        items.push_back(item);
    //
    //        if (plotOn) {
    //            std::stringstream ss;
    //            ss << "scan-" << std::setfill('0') << std::setw(5) << (i);
    //            g_gnuplot.setWinParam(ss.str(), std::make_pair(0.0, 20.0), std::make_pair(-10.0, 20.0));
    //            g_gnuplot.addPointVector("scan", reader.scans[i].scan.points, 13, 1, 2);
    //            g_gnuplot.addPointVector("features", items[i].features, 7, 1, 1);
    //            g_gnuplot.plot();
    //        }
    //
    //        looprate.sleep();
    //    }

    return 0;
}

void plotPoints(std::ostream& out, const VectorPoint2d& points) {
    for (int i = 0; i < points.size(); ++i) {
        out << points[i].x() << " " << points[i].y() << std::endl;
    }
}

void plotCompatible(std::ostream& out, const VectorPoint2d& points, const std::vector<std::pair<int, int> >& associations, bool useFirst) {
    int n = associations.size();
    int i1, i2;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (useFirst) {
                i1 = associations[i].first;
                i2 = associations[j].first;
            } else {
                i1 = associations[i].second;
                i2 = associations[j].second;
            }
            if (0 <= i1 && i1 < points.size() && 0 <= i2 && i2 < points.size()) {
                out << points[i1].x() << " " << points[i1].y() << "\n"
                        << points[i2].x() << " " << points[i2].y() << "\n"
                        << "\n";
            }
        }
    }
}

void plotAssociations(std::ostream& out, const VectorPoint2d& points1, const VectorPoint2d& points2, const std::vector<std::pair<int, int> >& associations) {
    for (auto& a : associations) {
        out << points1[a.first].x() << " " << points1[a.first].y() << "\n"
                << points2[a.second].x() << " " << points2[a.second].y() << "\n"
                << "\n";
    }
}
