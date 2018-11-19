#include <grd/TfUtils.h>
#include <grd/CorrespondenceGraph.h>

namespace grd {

    namespace TfUtils {

        //----------------------------------------------------------------------
        // CONVERSION UTILS
        //----------------------------------------------------------------------

        void vector3dToAffine2d(const Eigen::Vector3d& x, Transformation2d& t) {
            t = Transformation2d::Identity();
            t.prerotate(Eigen::Rotation2Dd(x[2]));
            t.pretranslate(Eigen::Vector2d(x[0], x[1]));
        }

        void vector3dToAffine2d(double x, double y, double theta, Transformation2d& t) {
            t = Transformation2d::Identity();
            t.prerotate(Eigen::Rotation2Dd(theta));
            t.pretranslate(Eigen::Vector2d(x, y));
        }

        void vector3dToAffine2dPtr(const Eigen::Vector3d& x, Transformation2d* t) {
            *t = Transformation2d::Identity();
            t->prerotate(Eigen::Rotation2Dd(x[2]));
            t->pretranslate(Eigen::Vector2d(x[0], x[1]));
        }

        Transformation2d vector3dToAffine2d(const Eigen::Vector3d& x) {
            Transformation2d t;
            vector3dToAffine2d(x, t);
            return t;
        }

        void affine2dToVector3d(const Transformation2d& t, Eigen::Vector3d& x) {
            x << t.translation()[0], t.translation()[1], std::atan2(t.matrix()(1, 0), t.matrix()(0, 0));
        }

        Eigen::Vector3d affine2dToVector3d(const Transformation2d& t) {
            Eigen::Vector3d x;
            affine2dToVector3d(t, x);
            return x;
        }

        void vector3dToIsometry2d(const Eigen::Vector3d& x, Eigen::Isometry2d& t) {
            t = Eigen::Isometry2d::Identity();
            t.prerotate(Eigen::Rotation2Dd(x[2]));
            t.pretranslate(Eigen::Vector2d(x[0], x[1]));
        }

        void vector3dToIsometry2d(double x, double y, double theta, Eigen::Isometry2d& t) {
            t = Eigen::Isometry2d::Identity();
            t.prerotate(Eigen::Rotation2Dd(theta));
            t.pretranslate(Eigen::Vector2d(x, y));
        }

        void vector3dToIsometry2dPtr(const Eigen::Vector3d& x, Eigen::Isometry2d* t) {
            *t = Eigen::Isometry2d::Identity();
            t->prerotate(Eigen::Rotation2Dd(x[2]));
            t->pretranslate(Eigen::Vector2d(x[0], x[1]));
        }

        Eigen::Isometry2d vector3dToIsometry2d(const Eigen::Vector3d& x) {
            Eigen::Isometry2d t;
            vector3dToIsometry2d(x, t);
            return t;
        }

        void Isometry2dToVector3d(const Eigen::Isometry2d& t, Eigen::Vector3d& x) {
            x << t.translation()[0], t.translation()[1], std::atan2(t.matrix()(1, 0), t.matrix()(0, 0));
        }

        Eigen::Vector3d Isometry2dToVector3d(const Eigen::Isometry2d& t) {
            Eigen::Vector3d x;
            Isometry2dToVector3d(t, x);
            return x;
        }

        //----------------------------------------------------------------------
        // TRANSFORMATION  UTILS
        //----------------------------------------------------------------------

        double computeTransform(const VectorPoint2d& points1, const VectorPoint2d& points2,
                const std::vector<std::pair<int, int> >& indices, Transformation2d& frame1Tframe2) {
            Point2d t1 = Eigen::Vector2d::Zero();
            Point2d t2 = Eigen::Vector2d::Zero();
            Matrix2d S = Eigen::Matrix2d::Zero();
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
            frame1Tframe2 = Transformation2d::Identity();
            frame1Tframe2.prerotate(rot);
            frame1Tframe2.pretranslate(transl);

            double residual = 0.0;
            for (int i = 0; i < (int) indices.size(); ++i) {
                if (0 <= indices[i].first && indices[i].first < (int) points1.size() &&
                        0 <= indices[i].second && indices[i].second < (int) points2.size()) {
                    residual += (points1[indices[i].first] - frame1Tframe2 * points2[indices[i].second]).norm();
                }
            }
            residual = residual / n;
            return residual;
        }

        void computeNNAssociation(const VectorPoint2d& points1, const VectorPoint2d& points2, double distMax, const Transformation2d& frame1Tframe2, std::vector<std::pair<int, int> >& indices) {
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

        void distAffine2d(const Transformation2d& pose1, const Transformation2d& pose2, double& distPosition, double& distAngle) {
            Eigen::Vector3d vec1 = TfUtils::affine2dToVector3d(pose1);
            Eigen::Vector3d vec2 = TfUtils::affine2dToVector3d(pose2);
            distPosition = (vec1.block<2, 1>(0, 0) - vec2.block<2, 1>(0, 0)).norm();
            distAngle = std::abs(vec1(2) - vec2(2));
            if (distAngle > M_PI) {
                distAngle = 2 * M_PI - distAngle;
            }
        }

        double associateTwoSteps(const VectorPoint2d& inputs, const VectorPoint2d& targets,
                double distCorrespTol, double distCorrespMin, double distNNMax,
                std::vector<std::pair<int, int> >& assocPP, Transformation2d& inputsTtargets) {
            //Transformation2d inputTtargetTrue, targetTinputTrue;
            Transformation2d targetTinputCG, inputTtargetCG;
            Transformation2d targetTinputNN, inputTtargetNN;
            std::vector<std::pair<int, int> > assocCG;
            std::vector<std::pair<int, int> > assocNN;
            double residual;

            // Associates the features between the items and estimates the transformation
            CorrespondenceGraph corrGraph;
            corrGraph.setDistanceToll(distCorrespTol);
            corrGraph.setDistanceMin(distCorrespMin);
            corrGraph.setInputs(inputs);
            corrGraph.setTargets(targets);
            corrGraph.computeTransform(targetTinputCG, assocCG);
            inputTtargetCG = targetTinputCG.inverse();
            //        std::cout << "   CG associated " << i1 << " (" << items[i1].features.size() << " features)"
            //                << " vs " << i2 << " (" << items[i2].features.size() << " features): " << assocCorrGraph.size() << std::endl;

            // Recomputes association using NN not to 
            computeNNAssociation(inputs, targets, distNNMax, inputTtargetCG, assocNN);
            residual = computeTransform(inputs, targets, assocNN, inputTtargetNN);
            targetTinputNN = inputTtargetNN.inverse();

            assocPP = assocNN;
            inputsTtargets = inputTtargetNN;
            return residual;
        }

    } // end of namespace 

} // end of namespace 
