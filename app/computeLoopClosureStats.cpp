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
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <queue>
#include <vector>
#include <boost/filesystem/path.hpp>

#include <grd/GnuplotVisualizer.h>
#include <grd/TfUtils.h>
#include <grd/ParamMap.h>

#define PRINT_OUT_VARIABLE(X) std::cout << "  " <<  #X << ": " << X << std::endl;

struct LocalizationResult {
    std::string name;
    int id1;
    int id2;
    int matchedNum;
    double residual;
    double errPos;
    double errAng;
    bool valid;
};

struct LocalizationStatGroup {
    std::string name;
    std::vector<LocalizationResult> results;
};

class LocalizationStatCollector {
public:
    std::vector<LocalizationStatGroup> statGroups;
    int match_num_thr_min;
    int match_num_thr_max;
    double err_pos_valid;
    double err_ang_valid;

    LocalizationStatCollector(grd::ParamMap& params) {
        std::string stat_filename;

        params.getParam<std::string>("stat_filename", stat_filename, "");
        params.getParam<int>("match_num_thr_min", match_num_thr_min, 0);
        params.getParam<int>("match_num_thr_max", match_num_thr_max, 25);
        params.getParam<double>("err_pos_valid", err_pos_valid, 0.5);
        params.getParam<double>("err_ang_valid_deg", err_ang_valid, 10.0);
        err_ang_valid = (M_PI / 180.0) * err_ang_valid;

        std::cout << "\nParams:\n";
        params.write(std::cout);
        std::cout << std::endl;

        readFile(stat_filename);

        boost::filesystem::path stat_path(stat_filename);
        std::cout << "\nstat_path: \"" << stat_path.string() << "\"\n"
                << "stat_path.parent_path(): \"" << stat_path.parent_path().string() << "\"\n"
                << "stat_path.stem(): \"" << stat_path.stem().string() << "\"\n"
                << std::endl;
        for (auto& group : statGroups) {
            // Name of stat file
            std::stringstream ss;
            //ss << stat_path.parent_path().string() << "/" << stat_path.stem().string() << "_" << group.name << ".dat";
            ss << "pr_" << stat_path.stem().string() << "_" << group.name << ".dat";
            std::cout << "Save stat for " << group.name << " in file \"" << ss.str() << "\"" << std::endl;
            // Open new file
            std::ofstream statFile(ss.str());
            computePrecisionRecall(group, statFile);
            statFile.close();
        }
    }

    void readFile(std::string filename) {
        std::string line;
        LocalizationResult item;

        std::ifstream statFile(filename.c_str());
        if (!statFile) {
            std::cerr << __PRETTY_FUNCTION__ << ": cannot open stat file \"" << filename << "\"" << std::endl;
            return;
        }
        std::cout << "Reading file \"" << filename << "\"" << std::endl;

        while (!statFile.eof()) {
            std::getline(statFile, line);
            if (line.length() > 0) {
                std::stringstream ss(line);
                if (ss >> item.id1) {
                    for (int i = 0; readLocalizationResult(ss, item); ++i) {
                        // If the stat group has not been created, it adds it to statGroups
                        if (i < statGroups.size()) {
                            if (statGroups[i].name != item.name) {
                                std::cerr << "Warning: Potential inconsistency stat group name: "
                                        << "\"" << item.name << "\" instead of \"" << statGroups[i].name << "\"" << std::endl;
                            }
                            statGroups[i].results.push_back(item);
                        } else {
                            LocalizationStatGroup group;
                            group.name = item.name;
                            group.results.push_back(item);
                            statGroups.push_back(group);
                            std::cout << "added group \"" << group.name << "\"" << std::endl;
                        }
                    }
                }
            }
        }
        statFile.close();
        for (auto& group : statGroups) {
            std::cout << "group " << group.name << " " << group.results.size() << " results" << std::endl;
        }
    }

    bool readLocalizationResult(std::istream& in, LocalizationResult& locres) {
        if (in >> locres.name >> locres.id2 >> locres.matchedNum >> locres.residual >> locres.errPos >> locres.errAng) {
            locres.valid = true;
        } else {
            locres.valid = false;
        }
        return locres.valid;
    }

    void computePrecisionRecall(const LocalizationStatGroup& group, std::ostream& out) {
        double tp, tn, fp, fn;
        double errPosAvg, errAngAvg, precision, recall;
        double precisionPrev = 0.0, recallPrev = 0.0;

        std::cout << "Compute stats for \"" << group.name << "\" num threshold " << match_num_thr_min << " -> " << match_num_thr_max << std::endl;
        out << "# Thres \tTP \tTN \tFP \tFN \trecall \tprecision \tErrPosAvg[m] \tErrAngAvg[deg]" << std::endl;
        for (int thr = match_num_thr_min; thr < match_num_thr_max; ++thr) {
            tp = tn = fp = fn = 0;
            errPosAvg = errAngAvg = 0.0;
            for (auto& item : group.results) {
                // Positive: item.matchedNum >= thr
                // Negative: item.matchedNum < thr
                if (item.matchedNum >= thr) {
                    if (item.errPos < err_pos_valid && item.errAng < err_ang_valid) {
                        errPosAvg += item.errPos;
                        errAngAvg += item.errAng;
                        tp++;
                    } else {
                        fp++;
                    }
                } else {
                    if (item.errPos < err_pos_valid && item.errAng < err_ang_valid) {
                        fn++;
                    } else {
                        tn++;
                    }
                }
            }
            if (tp + tn > 0) {
                recall = 1.0 * tp / (tp + fn);
            } else {
                recall = 0.0;
            }
            if (tp + fp > 0) {
                precision = 1.0 * tp / (tp + fp);
            } else {
                precision = 0.0;
            }
            if (tp > 0) {
                errPosAvg = errPosAvg / tp;
                errAngAvg = errAngAvg / tp;
            }
            if (precision >= precisionPrev || recall >= recallPrev) {
                out << thr << " \t" << std::setprecision(6) << tp << " \t" << tn << " \t" << fp << " \t" << fn << " \t"
                        << std::setprecision(3) << recall << " \t"
                        << std::setprecision(3) << precision << " \t"
                        << std::setprecision(4) << errPosAvg << " \t"
                        << std::setprecision(2) << (180.0 / M_PI * errAngAvg)
                        << std::endl;
                precisionPrev = precision;
                recallPrev = recall;
            }
        }
    }
};

int main(int argc, char** argv) {
    grd::ParamMap params;
    std::string filenameIn;
    std::string filenameCfg;

    // Reads file params from command line
    params.read(argc, argv);
    params.getParam<std::string>("in", filenameIn, std::string(""));
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));

    // Reads potential parameter value from file and, then again, from command line
    if (!params.read(filenameCfg)) {
        std::cout << "Cannot open configuration file \"" << filenameCfg << "\": using default values" << std::endl;
    }
    params.read(argc, argv);

    LocalizationStatCollector lsc(params);

    std::cout << "Closing " << __FILE__ << std::endl;
    return (0);
}
