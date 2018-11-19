#include <grd/GnuplotVisualizer.h>

namespace grd {

    GnuplotVisualizer::GnuplotVisualizer() :
    gp("gnuplot -rv -persist"), first(true), disabled(false), tableEnabled(false) {
        //		gp.set_auto_close(true);
        std::cout << __FILE__ << "," << __LINE__ << ": first " << first << ", labels.size() " << labels.size() << std::endl;
    };

    GnuplotVisualizer::GnuplotVisualizer(const GnuplotVisualizer& orig) {
    }

    GnuplotVisualizer::~GnuplotVisualizer() {
        gp.clearTmpfiles();
        gp.close();
    };

    void GnuplotVisualizer::createFilePng(std::string filename, std::string title) {
        //    gp << "set terminal pngcairo size 1024,768 enhanced font 'Times-Roman,14'\n"
        //       << "set grid back ls 12 \n"
        //            << "set output \"" << filename << "\"\n";
        gp << "set title \"" << title << "\"" << std::endl;
    }

    void GnuplotVisualizer::createFileEps(std::string filename, std::string title) {
        gp << "set terminal postscript eps size 3.5,3.5 enhanced color font 'Arial,14'\n"
                << "set output \"" << filename << "\"\n";
        gp << "set title \"" << title << "\"" << std::endl;
    }

    void GnuplotVisualizer::setTable(std::string filename) {
        gp << "set table '" << filename << "'\n";
        tableEnabled = true;
    }

    void GnuplotVisualizer::setTerminalWxt(int wxtNum) {
        gp << "set term wxt " << wxtNum << "\n";
    }

    bool GnuplotVisualizer::checkInit() {
        if (disabled) return false;
        //std::cout << __FILE__ << "," << __LINE__ << ": first " << first << std::endl;
        if (first) {
            first = false;
            //			gp << "unset label" << std::endl;
            gp << "plot ";
        } else {
            gp << ", ";
        }
        return true;
    }

    void GnuplotVisualizer::addDownLeftText(const std::string& text, double xOffset, double yOffset) {
        if (disabled) return;
        gp << std::endl;
        gp << "set label 100000 \"" << text << "\" at graph " << xOffset << "," << yOffset << " left" << std::endl;
    }

    void GnuplotVisualizer::plot() {
        if (disabled) return;
        //std::cout << "first " << first << ", labels.size() " << labels.size() << std::endl;
        gp << std::endl;
        for (int i = 0; i < labels.size(); ++i) {
            gp << "set label \"" << labels[i].get<0>() << "\" at " << labels[i].get<1>() << "," << labels[i].get<2>() << " tc lt " << labels[i].get<3>() << std::endl;
        }
        labels.clear();
        first = true;

        if (tableEnabled) {
            gp << "unset table\n";
            tableEnabled = false;
        }
    }

    void GnuplotVisualizer::setWinParam(std::string t, const std::pair<float, float>& xrange,
            const std::pair<float, float>& yrange) {
        if (disabled) return;
        gp << "set term wxt size 1600,900 title '" << t << "' noraise\n";
        //gp << "set term postscript eps enhanced color font 'Helvetica,10'\n";
        //gp << "set output \"" << t << ".eps\"\n";
        gp << "set size ratio -1\n";
        gp << "unset key\n";
        if (std::abs(xrange.first - xrange.second) > 1e-5) {
            gp << "set xrange [" << xrange.first << ":" << xrange.second << "]\n";
        }
        if (std::abs(yrange.first - yrange.second) > 1e-5) {
            gp << "set yrange [" << yrange.first << ":" << yrange.second << "]\n";
        }
    }

    GnuplotVisualizer g_gnuplot;

} // end of namespace 