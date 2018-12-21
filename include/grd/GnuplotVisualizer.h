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
#ifndef GRD_GNUPLOT_VISUALIZER_H
#define GRD_GNUPLOT_VISUALIZER_H

#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <boost/lexical_cast.hpp>

#include <grd/thirdparty/gnuplot-iostream.h>

namespace grd {

    enum class PlotType {
        Point,
        Line,
        Vector
    };

    class GnuplotVisualizer {
    public:
        Gnuplot gp;

        GnuplotVisualizer();

        GnuplotVisualizer(const GnuplotVisualizer& orig);

        virtual ~GnuplotVisualizer();

        void createFilePng(std::string filename, std::string title);

        void createFileEps(std::string filename, std::string title);

        void setTable(std::string filename);

        void setTerminalWxt(int wxtNum);

        bool checkInit();

        template <typename VectorT>
        void addPointVector(const std::string& name, const VectorT& v, int pt = 1, double ps = 1, int lt = 1) {
            //std::cout << __FILE__ << "," << __LINE__ << "," << __FUNCTION__ << ": adding data \"" << name << "\", first " << first << std::endl;
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double> > points;
            for (auto& el : v) {
                points.push_back(boost::make_tuple(el(0), el(1)));
            }
            gp << gp.file1d(points, "/tmp/" + name) << " w p ";
            gp << "pt " << pt << " ";
            gp << "ps " << ps << " ";
            gp << "lt " << lt << " ";
            //		gp << "title '" << name << "' ";
            gp << "notitle ";
        }

        template <typename T>
        void addPointVectorLines(const std::string& name, const std::vector<T>& v, int pt = 1, int ps = 1, int lt = 1) {
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double> > points;
            for (auto& el : v) {
                points.push_back(boost::make_tuple(el(0), el(1)));
            }
            gp << gp.file1d(points, "/tmp/" + name) << " w l ";
            //		gp << "lt " << pt << " ";
            //gp << "ps " << ps << " ";
            gp << "lc " << lt << " ";
            gp << "title '" << name << "' ";
        }

        template <typename T>
        void addPointVectorLinesBezier(const std::string& name, const std::vector<T>& v, int pt = 1, int ps = 1, int lt = 1) {
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double> > points;
            for (auto& el : v) {
                points.push_back(boost::make_tuple(el(0), el(1)));
            }
            gp << gp.file1d(points, "/tmp/" + name) << " smooth cspline ";
            //		gp << "lt " << pt << " ";
            //gp << "ps " << ps << " ";
            gp << "lc " << lt << " ";
            gp << "title '" << name << "' ";
        }

        template <typename T>
        void addPointVectorDottedLines(const std::string& name, const std::vector<T>& v, int pt = 1, int pi = 1, int lc = 1) {
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double> > points;
            for (auto& el : v) {
                points.push_back(boost::make_tuple(el(0), el(1)));
            }
            gp << gp.file1d(points, "/tmp/" + name) << " with linespoints ";
            gp << "lt " << -1 << " ";
            gp << "pt " << pt << " ";
            gp << "pi " << pi << " ";
            gp << "lc " << lc << " ";
            gp << "title '" << name << "' ";
        }

        template <typename T, typename R>
        void addPointVector(const std::string& name, const std::vector<T>& v, const R& trans, double pt = 1, double ps = 1, double lt = 1) {
            //std::cout << __FILE__ << "," << __LINE__ << "," << __FUNCTION__ << ": adding data \"" << name << "\", first " << first << std::endl;
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double> > points;
            for (auto& el : v) {
                const T& p = trans * el;
                points.push_back(boost::make_tuple(p(0), p(1)));
            }
            gp << gp.file1d(points, "/tmp/" + name) << " w p ";
            gp << "pt " << pt << " ";
            gp << "ps " << ps << " ";
            gp << "lt " << lt << " ";
            //		gp << "title '" << name << "' ";
            gp << "notitle ";
        }

        template <typename T>
        void addPoint(const std::string& name, const T& p, int pt = 1, int ps = 1, int lt = 1) {
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double> > zero;
            zero.push_back(boost::make_tuple(p(0), p(1)));
            //		for(auto& el : v){
            //				points.push_back(boost::make_tuple(-el(1), el(0)));
            //		}
            gp << gp.file1d(zero, "/tmp/" + name) << " w p ";
            gp << "pt " << pt << " ";
            gp << "ps " << ps << " ";
            gp << "lt " << lt << " ";
            //		gp << "title '" << name << "' ";
            gp << "notitle ";
        }

        template <typename T>
        void addTrajectory(const std::string& name, const std::vector<T>& v, const T& trans, double lc = 1) {
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double> > points;
            for (auto& el : v) {
                const Eigen::Vector2d& p = trans * el.translation();
                points.push_back(boost::make_tuple(p(0), p(1)));
            }
            gp << gp.file1d(points, "/tmp/" + name) << " w l ";
            gp << "lc " << lc << " ";
            //		gp << "title '" << name << "' ";
            gp << "notitle ";
        }

        template <typename T>
        void addTrajectory(const std::string& name, const std::vector<T>& v, const T& trans, const std::vector<bool>& boolcolor, double lc1 = 2, double lc2 = 1) {
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double, double> > points;
            for (int i = 0; i < v.size(); ++i) {
                const T& el = v[i];
                const Eigen::Vector2d& p = trans * el.translation();
                points.push_back(boost::make_tuple(p(0), p(1), (boolcolor[i] ? lc1 : lc2)));
            }
            gp << gp.file1d(points, "/tmp/" + name) << " w l ";
            gp << "lc var ";
            //		gp << "title '" << name << "' ";
            gp << "notitle ";
        }

        template <typename T>
        void addRect(const std::string& name, const T& min, const T& max, double lc = 1) {
            if (!checkInit()) return;
            std::vector<boost::tuple<double, double> > points;
            //		for (auto& el : v) {
            //			const Eigen::Vector2d& p = trans * el.translation();
            //			points.push_back(boost::make_tuple(p(0), p(1)));
            //		}
            points.push_back(boost::make_tuple(min(0), min(1)));
            points.push_back(boost::make_tuple(max(0), min(1)));
            points.push_back(boost::make_tuple(max(0), max(1)));
            points.push_back(boost::make_tuple(min(0), max(1)));
            points.push_back(boost::make_tuple(min(0), min(1)));
            gp << gp.file1d(points, "/tmp/" + name) << " w l ";
            gp << "lc " << lc << " ";
            //		gp << "title '" << name << "' ";
            gp << "notitle ";
        }

        template <typename VectorT>
        void addOrientedPointVector(const std::string& name, const VectorT& v, int lt = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double> > circles;
            for (auto& el : v) {
                circles.push_back(boost::make_tuple(el.point(0), el.point(1), el.radius));
            }

            gp << gp.file1d(circles, "/tmp/" + name) << " w circle ";
            gp << "lt " << lt << " ";
            gp << "title '" << name << "', ";

            std::vector<boost::tuple<double, double, double, double> > vecs;
            //		for (auto& el : v) {
            //			Eigen::Affine2d tr;
            //			tr.prerotate(Eigen::Rotation2Df())
            //			Eigen::Rotation2Df rot(theta);
            //		Eigen::Vector2f transl
            //			vecs.push_back(boost::make_tuple(el(0), el(1), radius));
            //		}

            decltype(v[0].point) versor, vec;
            for (int i = 0; i < v.size(); ++i) {
                Eigen::Affine2d tr = Eigen::Affine2d::Identity();
                tr.prerotate(Eigen::Rotation2Dd(v[i].orientation));
                //			tr.pretranslate(v[i]);
                versor << v[i].radius, 0;
                vec = tr*versor;
                vecs.push_back(boost::make_tuple(v[i].point(0), v[i].point(1), vec(0), vec(1)));
            }

            gp << gp.file1d(vecs, "/tmp/vec_" + name) << " w vector ";
            gp << "lt " << lt << " ";
            gp << "title 'vec_" << name << "' ";
        }

        template <typename T>
        void addOrientedPointVector(const std::string& name, const std::vector<T>& v, const std::vector<double>& o, double radius = 1, int lt = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double> > circles;
            for (auto& el : v) {
                circles.push_back(boost::make_tuple(el(0), el(1), radius));
            }

            gp << gp.file1d(circles, "/tmp/" + name) << " w circle ";
            gp << "lt " << lt << " ";
            gp << "title '" << name << "', ";

            T versor;
            versor << radius, 0;
            std::vector<boost::tuple<double, double, double, double> > vecs;
            //		for (auto& el : v) {
            //			Eigen::Affine2d tr;
            //			tr.prerotate(Eigen::Rotation2Df())
            //			Eigen::Rotation2Df rot(theta);
            //		Eigen::Vector2f transl
            //			vecs.push_back(boost::make_tuple(el(0), el(1), radius));
            //		}

            for (int i = 0; i < v.size(); ++i) {
                Eigen::Affine2d tr = Eigen::Affine2d::Identity();
                tr.prerotate(Eigen::Rotation2Dd(o[i]));
                //			tr.pretranslate(v[i]);
                T vec = tr*versor;
                vecs.push_back(boost::make_tuple(v[i](0), v[i](1), vec(0), vec(1)));
            }

            gp << gp.file1d(vecs, "/tmp/vec_" + name) << " w vector ";
            gp << "lt " << lt << " ";
            gp << "title 'vec_" << name << "' ";
        }

        template <typename T>
        void addTraiectoryPairLines(const std::string& name, const std::vector<T>& v1, const std::vector<T>& v2, const T& trans, int lc = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double, double> > vecs;
            for (int i = 0; i < v1.size(); ++i) {
                const Eigen::Vector2d& p1 = trans * v1[i].translation();
                const Eigen::Vector2d& p2 = trans * v2[i].translation();
                vecs.push_back(boost::make_tuple(p1(0), p1(1), p2(0) - p1(0), p2(1) - p1(1)));
            }

            gp << gp.file1d(vecs, "/tmp/" + name) << " w vectors nohead ";
            gp << "lc " << lc << " ";
            gp << "title '" << name << "' ";
        }

        template <typename T>
        void addOrientedPointVector(const std::string& name, const std::vector<T>& v, const std::vector<double>& o, const std::vector<double>& r, int lt = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double> > circles;
            for (int i = 0; i < v.size(); ++i) {
                circles.push_back(boost::make_tuple(v[i](0), v[i](1), r[i]));
            }

            gp << gp.file1d(circles, "/tmp/" + name) << " w circle ";
            gp << "lt " << lt << " ";
            gp << "title '" << name << "', ";


            std::vector<boost::tuple<double, double, double, double> > vecs;
            //		for (auto& el : v) {
            //			Eigen::Affine2d tr;
            //			tr.prerotate(Eigen::Rotation2Df())
            //			Eigen::Rotation2Df rot(theta);
            //		Eigen::Vector2f transl
            //			vecs.push_back(boost::make_tuple(el(0), el(1), radius));
            //		}

            for (int i = 0; i < v.size(); ++i) {
                Eigen::Affine2d tr = Eigen::Affine2d::Identity();
                tr.prerotate(Eigen::Rotation2Dd(o[i]));
                //			tr.pretranslate(v[i]);
                T versor;
                versor << r[i], 0;
                T vec = tr*versor;
                vecs.push_back(boost::make_tuple(v[i](0), v[i](1), vec(0), vec(1)));
            }

            gp << gp.file1d(vecs, "/tmp/vec_" + name) << " w vector ";
            gp << "lt " << lt << " ";
            gp << "title 'vec_" << name << "' ";
        }

        template <typename T>
        void addOrientedPoint(const std::string& name, const T& v, const double o, double radius = 1, int lt = 1) {
            if (!checkInit()) return;
            //
            std::vector<boost::tuple<double, double, double> > circles;
            //				for (auto& el : v) {
            circles.push_back(boost::make_tuple(v(0), v(1), radius));
            //				}

            gp << gp.file1d(circles, "/tmp/circ_" + name) << " w circle ";
            gp << "lt " << lt << " ";
            gp << "title '" << name << "', ";

            T versor;
            versor << radius, 0;
            std::vector<boost::tuple<double, double, double, double> > vecs;

            Eigen::Affine2d tr = Eigen::Affine2d::Identity();
            tr.prerotate(Eigen::Rotation2Dd(o));
            //			tr.pretranslate(v[i]);
            T vec = tr*versor;
            vecs.push_back(boost::make_tuple(v(0), v(1), vec(0), vec(1)));

            gp << gp.file1d(vecs, "/tmp/vec_" + name) << " w vector ";
            gp << "lt " << lt << " ";
            gp << "title 'vec_" << name << "' ";
        }

        template <typename T>
        void addTextPoint(const T& v, const std::string& text, int lt = 1) {
            if (disabled) return;
            labels.push_back(boost::make_tuple(text, v(0), v(1), lt));
        }

        template <typename T>
        void addCircle(const std::string& name, const T& p, double radius = 1, int lt = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double> > circles;
            circles.push_back(boost::make_tuple(p(0), p(1), radius));

            gp << gp.file1d(circles, "/tmp/" + name) << " w circle ";
            gp << "lt " << lt << " ";
            gp << "title '" << name << "' ";
        }

        template <typename T>
        void addHistogramCorner(const std::string& name, const T& p, const std::vector<T>& L, const std::vector<T>& R, double radius = 1, int slices = 16, int lt = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double> > circles;
            circles.push_back(boost::make_tuple(p(0), p(1), radius));

            gp << gp.file1d(circles, "/tmp/" + name) << " w circle ";
            gp << "lt " << 9 << " ";
            gp << "title '" << name << "', ";


            T versor;
            versor << radius, 0;
            std::vector<boost::tuple<double, double, double, double> > slice;

            for (int i = 0; i < slices; ++i) {
                Eigen::Affine2d tr = Eigen::Affine2d::Identity();
                tr.prerotate(Eigen::Rotation2Dd(i * 2 * M_PI / slices));
                //			tr.pretranslate(v[i]);
                T vec = tr*versor;
                slice.push_back(boost::make_tuple(p(0), p(1), vec(0), vec(1)));
            }

            gp << gp.file1d(slice, "/tmp/vec_" + name) << " w vector ";
            gp << "nohead lt " << 9 << " ";
            gp << "title 'slice_" << name << "', ";

            std::vector<boost::tuple<double, double, double, double> > vecs;
            vecs.push_back(boost::make_tuple(p(0), p(1), L[L.size() - 1](0) - p(0), L[L.size() - 1](1) - p(1)));
            for (int i = L.size() - 2; i >= 0; --i) {
                //if ((p - L[i]).norm() < radius)
                vecs.push_back(boost::make_tuple(p(0), p(1), L[i](0) - p(0), L[i](1) - p(1)));
            }
            gp << gp.file1d(vecs, "/tmp/vecL_" + name) << " w vector ";
            gp << "nohead lt " << lt << " ";
            gp << "title 'vec_" << name << "', ";

            vecs.clear();
            vecs.push_back(boost::make_tuple(p(0), p(1), R[0](0) - p(0), R[0](1) - p(1)));
            for (int i = 1; i < R.size(); ++i) {
                //if ((p - R[i]).norm() < radius)
                vecs.push_back(boost::make_tuple(p(0), p(1), R[i](0) - p(0), R[i](1) - p(1)));
            }
            gp << gp.file1d(vecs, "/tmp/vecR_" + name) << " w vector ";
            gp << "nohead lt " << lt + 3 << " ";
            gp << "title 'vec_" << name << "' ";
        }

        template <typename T>
        void addVectorMatrix(const std::string& name, const T& p) {
            if (!checkInit()) return;
            gp << gp.file1d(p, "/tmp/" + name) << " matrix ";
            gp << "title 'mat_" << name << "' ";
        }

        template <typename T>
        void addHistogramDescriptor(const std::string& name, const T& p, const std::vector<double>& h, double radius = 1, int slices = 16, int lt = 1, int lw = 2) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double> > circles;
            circles.push_back(boost::make_tuple(p(0), p(1), radius));

            gp << gp.file1d(circles, "/tmp/" + name) << " w circle ";
            gp << "lt " << 9 << " ";
            gp << "lw " << lw << " ";
            gp << "title '" << name << "', ";


            T versor;
            versor << radius, 0;
            std::vector<boost::tuple<double, double, double, double> > slice;
            assert(h.size() == slices && "error with histogram lenght");
            for (int i = 0; i < slices; ++i) {
                Eigen::Affine2d tr = Eigen::Affine2d::Identity();
                tr.prerotate(Eigen::Rotation2Dd(i * 2 * M_PI / slices));
                //			tr.pretranslate(v[i]);
                T vec = tr*versor;
                slice.push_back(boost::make_tuple(p(0), p(1), vec(0), vec(1)));

                //			labels.push_back(boost::make_tuple(boost::lexical_cast<std::string>(h[i]), p(0) + vec(0), p(1) + vec(1), lt));
            }

            gp << gp.file1d(slice, "/tmp/vec_" + name) << " w vectors ";
            gp << "nohead lt " << 9 << " ";
            gp << "lw " << lw << " ";
            gp << "title 'slice_" << name << "' ";


        }

        template <typename T>
        void addBinaryDescriptor(const std::string& name, const T& p, double radius = 1, int slices = 24, int circles = 8, int lt = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double> > cir;
            double rstep = radius / circles;
            for (int i = 0; i < circles; ++i) {
                cir.push_back(boost::make_tuple(p(0), p(1), rstep * (i + 1)));
            }

            gp << gp.file1d(cir, "/tmp/" + name) << " w circle ";
            gp << "lt " << 9 << " ";
            gp << "title '" << name << "', ";


            T versor;
            versor << radius, 0;
            std::vector<boost::tuple<double, double, double, double> > slice;
            //		assert(h.size() == slices && "error with histogram lenght");
            for (int i = 0; i < slices; ++i) {
                Eigen::Affine2d tr = Eigen::Affine2d::Identity();
                tr.prerotate(Eigen::Rotation2Dd(i * 2 * M_PI / slices));
                //			tr.pretranslate(v[i]);
                T vec = tr*versor;
                slice.push_back(boost::make_tuple(p(0), p(1), vec(0), vec(1)));

                //			labels.push_back(boost::make_tuple(boost::lexical_cast<std::string>(h[i]), p(0) + vec(0), p(1) + vec(1), lt));
            }

            gp << gp.file1d(slice, "/tmp/vec_" + name) << " w vector ";
            gp << "nohead lt " << 9 << " ";
            gp << "title 'slice_" << name << "' ";


        }

        template <typename T>
        void addTriangularCorner(const std::string& name, const T& p, const T& L, const T& R, double radius = 1, int slices = 16, int lt = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double, double> > vecs;
            vecs.push_back(boost::make_tuple(p(0), p(1), L(0) - p(0), L(1) - p(1)));
            //		for (int i = L.size()-2; i >= 0; --i) {
            //			//if ((p - L[i]).norm() < radius)
            //				vecs.push_back(boost::make_tuple(L[i+1](0), L[i+1](1), L[i](0) - L[i+1](0), L[i](1) - L[i+1](1)));
            //		}
            gp << gp.file1d(vecs, "/tmp/vecL_" + name) << " w vector ";
            gp << "nohead lt " << lt << " ";
            gp << "title 'vecL_" << name << "', ";

            vecs.clear();
            vecs.push_back(boost::make_tuple(p(0), p(1), R(0) - p(0), R(1) - p(1)));
            //		for (int i = 1; i < R.size(); ++i) {
            //			//if ((p - R[i]).norm() < radius)
            //				vecs.push_back(boost::make_tuple(R[i-1](0), R[i-1](1), R[i](0) - R[i-1](0), R[i](1) - R[i-1](1)));
            //		}
            gp << gp.file1d(vecs, "/tmp/vecR_" + name) << " w vector ";
            gp << "nohead lt " << lt + 3 << " ";
            gp << "title 'vecR_" << name << "', ";

            vecs.clear();
            vecs.push_back(boost::make_tuple(L(0), L(1), R(0) - L(0), R(1) - L(1)));
            //		for (int i = 1; i < R.size(); ++i) {
            //			//if ((p - R[i]).norm() < radius)
            //				vecs.push_back(boost::make_tuple(R[i-1](0), R[i-1](1), R[i](0) - R[i-1](0), R[i](1) - R[i-1](1)));
            //		}
            gp << gp.file1d(vecs, "/tmp/vecD_" + name) << " w vector ";
            gp << "nohead lt " << 3 << " ";
            gp << "title 'vecD_" << name << "' ";
        }

        template <typename T>
        void addConstellationVector(const std::string& name, T p, const std::vector<T>& f2, int lt = 1, int lc = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double, double> > vecs;
            for (auto& i : f2) {
                vecs.push_back(boost::make_tuple(p(0), p(1), i(0) - p(0), i(1) - p(1)));
            }

            gp << gp.file1d(vecs, "/tmp/" + name) << " w vector ";
            gp << "nohead lt " << lt << " ";
            gp << "lc " << lc << " ";
            gp << "title 'vec_" << name << "' ";
        }

        template <typename VectorT>
        void addAssociationVector(const std::string& name, const VectorT& f1, const VectorT& f2, std::vector<std::pair<int, int> >& match, int lt = 1, int lc = 1, int lw = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double, double, double> > vecs;
            for (auto& i : match) {
                if (i.second >= 0) {
                    vecs.push_back(boost::make_tuple(f1[i.first](0), f1[i.first](1), f2[i.second](0) - f1[i.first](0), f2[i.second](1) - f1[i.first](1)));
                }
            }

            gp << gp.file1d(vecs, "/tmp/" + name) << " w vectors ";
            gp << "lt " << lt << " ";
            gp << "lw " << lw << " ";
            gp << "lc " << lc << " ";
            gp << "title 'vec_" << name << "' ";
        }

        template <typename T, typename R>
        void addPointCliqueLines(const std::string& name, const std::vector<T>& v, const R& trans, int lt = 1, int ps = 1, int lc = 1) {
            if (!checkInit()) return;
            //std::vector<boost::tuple<double, double, double, double> > lines;
            std::vector<boost::tuple<double, double> > lines;
            for (int i = 0; i < v.size(); ++i) {
                T vt1 = trans * v[i];
                for (int j = i + 1; j < v.size(); ++j) {
                    //lines.push_back(boost::make_tuple(v[i](0), v[i](1), v[j](0)-v[i](0), v[j](1)-v[i](1)));
                    T vt2 = trans * v[j];
                    lines.push_back(boost::make_tuple(vt1(0), vt1(1)));
                    lines.push_back(boost::make_tuple(vt2(0), vt2(1)));
                }
            }
            gp << gp.file1d(lines, "/tmp/" + name) << " w l ";
            gp << "lt " << lt << " ";
            //		gp << "ps " << ps << " ";
            gp << "lc " << lc << " ";
            gp << "title '" << name << "' ";
        }

        template <typename T>
        void addLine(const std::string& name, const T& p1, const T& p2, int lt = 1, int lc = 1) {
            if (!checkInit()) return;

            std::vector<boost::tuple<double, double> > points;
            points.push_back(boost::make_tuple(p1(0), p1(1)));
            points.push_back(boost::make_tuple(p2(0), p2(1)));

            gp << gp.file1d(points, "/tmp/" + name) << " w l ";
            gp << "lt " << lt << " ";
            gp << "lc " << lc << " ";
            gp << "title '" << name << "' ";
        }

        void addDownLeftText(const std::string& text, double xOffset, double yOffset);

        void plot();

        void setWinParam(std::string t, const std::pair<float, float>& xrange = std::make_pair(0.f, 0.f),
                const std::pair<float, float>& yrange = std::make_pair(0.f, 0.f));

        inline void setEnable(bool enable) {
            disabled = !enable;
        }


    private:
        bool first;
        std::vector<boost::tuple<std::string, double, double, double> > labels;
        //	std::string winTitle;
        bool disabled;
        bool tableEnabled;
    };

    extern GnuplotVisualizer g_gnuplot;
}

#endif