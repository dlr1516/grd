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
#ifndef GRD_TF_UTILS_H
#define GRD_TF_UTILS_H

#include <Eigen/Dense>
#include <grd/Point.h>

namespace grd {

    namespace TfUtils {

        /**
         * @brief converts from 3D vector to tf transorm
         * @param x vector in the form [x, y, theta]^T
         * @param tran transform in tf namespace
         */
        //	void vector3dToTf(const Eigen::Vector3d& x, tf::Transform& tran) {
        //		tran.setOrigin(tf::Vector3(x[0], x[1], 0.0));
        //		tf::Quaternion q;
        //		q.setRPY(0, 0, x[2]);
        //		tran.setRotation(q);
        //	}

        /**
         * @brief converts between 3D vector and tf tranform
         * @param x vector in the form [x, y, theta]^T
         * @return transform in tf namespace
         */
        //	tf::Transform vector3dToTf(const Eigen::Vector3d& x) {
        //		tf::Transform tran;
        //		vector3dToTf(x, tran);
        //		return tran;
        //	}

        /**
         * @brief converts from 3D vector to 2d affine matrix
         * @param x vector in the form [x, y, theta]^T
         * @param t 2d affine matrix
         */
        void vector3dToAffine2d(const Eigen::Vector3d& x, Transformation2d& t);

        /**
         * @brief converts from 3D vector to 2d affine matrix
         * @param x coordinate x
         * @param y coordinate y
         * @param theta angle theta
         * @param t 2d affine matrix
         */
        void vector3dToAffine2d(double x, double y, double theta, Transformation2d& t);

        /**
         * @brief converts from 3D vector to 2d affine matrix
         * @param x vector in the form [x, y, theta]^T
         * @param t 2d affine matrix
         */
        void vector3dToAffine2dPtr(const Eigen::Vector3d& x, Transformation2d* t);

        /**
         * @brief converts from 3D vector to 2d affine matrix
         * @param x vector in the form [x, y, theta]^T
         * @return 2d affine matrix
         */
        Transformation2d vector3dToAffine2d(const Eigen::Vector3d& x);

        /**
         * @brief converts from 2d affine matrix to 3D vector.
         * @param t 2d affine matrix
         * @return vector in the form [x, y, theta]^T
         */
        void affine2dToVector3d(const Transformation2d& t, Eigen::Vector3d& x);

        /**
         * @brief converts from 2d affine matrix to 3D vector.
         * @param t 2d affine matrix
         * @return vector in the form [x, y, theta]^T
         */
        Eigen::Vector3d affine2dToVector3d(const Transformation2d& t);

        /**
         * @brief converts from 3D vector to 2d isometry matrix
         * @param x vector in the form [x, y, theta]^T
         * @param t 2d isometry matrix
         */
        void vector3dToIsometry2d(const Eigen::Vector3d& x, Eigen::Isometry2d& t);

        /**
         * @brief converts from 3D vector to 2d isometry matrix
         * @param x coordinate x
         * @param y coordinate y
         * @param theta angle theta
         * @param t 2d isometry matrix
         */
        void vector3dToIsometry2d(double x, double y, double theta, Eigen::Isometry2d& t);

        /**
         * @brief converts from 3D vector to 2d isometry matrix
         * @param x vector in the form [x, y, theta]^T
         * @param t 2d isometry matrix
         */
        void vector3dToIsometry2dPtr(const Eigen::Vector3d& x, Eigen::Isometry2d* t);

        /**
         * @brief converts from 3D vector to 2d isometry matrix
         * @param x vector in the form [x, y, theta]^T
         * @return 2d isometry matrix
         */
        Eigen::Isometry2d vector3dToIsometry2d(const Eigen::Vector3d& x);

        /**
         * @brief converts from 2d isometry matrix to 3D vector.
         * @param t 2d isometry matrix
         * @return vector in the form [x, y, theta]^T
         */
        void Isometry2dToVector3d(const Eigen::Isometry2d& t, Eigen::Vector3d& x);

        /**
         * @brief converts from 2d isometry matrix to 3D vector.
         * @param t 2d isometry matrix
         * @return vector in the form [x, y, theta]^T
         */
        Eigen::Vector3d Isometry2dToVector3d(const Eigen::Isometry2d& t);

        /**
         * Computes the transformation from reference frame points1 to reference frame points2
         * @param points1 
         * @param points2
         * @param indices
         * @param transform
         * @return 
         */
        double computeTransform(const VectorPoint2d& points1, const VectorPoint2d& points2, const std::vector<std::pair<int, int> >& indices, Transformation2d& transform);

        /**
         * Computes the nearest neighbor association between points1 and points2. 
         * @param points1
         * @param points2
         * @param distMax
         * @param transform raw estimate of transformation between 
         * @param indices association pairs of indices (i1,i2) where i1 index of points1 and i2 of points2
         */
        void computeNNAssociation(const VectorPoint2d& points1, const VectorPoint2d& points2, double distMax, const Transformation2d& transform, std::vector<std::pair<int, int> >& indices);

        /**
         * Computes the distance between the origins of reference frames pose1 and pose2 and their relative rotation angle
         * @param pose1
         * @param pose2
         * @param distPosition
         * @param distAngle
         */
        void distAffine2d(const Transformation2d& pose1, const Transformation2d& pose2, double& distPosition, double& distAngle);

        /**
         * Computes point-to-point association between point sets inputs and targets 
         * using first CorrespondenceGraph method and then nearest neighbor policy 
         * to include other potential associations. 
         * @param inputs
         * @param targets
         * @param distCorrespTol
         * @param distCorrespMin
         * @param distNNMax
         * @param assocPP
         * @param inputsTtargets
         */
        double associateTwoSteps(const VectorPoint2d& inputs, const VectorPoint2d& targets,
                double distCorrespTol, double distCorrespMin, double distNNMax,
                std::vector<std::pair<int, int> >& assocPP, Transformation2d& inputsTtargets);

        /**
         * @brief converts from 2d affine matrix to tf transorm 
         * @param t 2d affine matrix
         * @param tran transform in tf namespace
         */
        //	void affine2dToTf(const Transformation2d& t, tf::Transform& tran) {
        //		tran.setOrigin(tf::Vector3(t.translation()[0], t.translation()[1], 0.0));
        //		tf::Quaternion q;
        //		q.setRPY(0, 0, atan2(t.rotation()(1, 0), t.rotation()(0, 0)));
        //		tran.setRotation(q);
        //	}

        /**
         * @brief converts from 2d affine matrix to tf transorm 
         * @param t 2d affine matrix
         * @return transform in tf namespace
         */
        //	tf::Transform affine2dToTf(const Transformation2d& t) {
        //		tf::Transform tran;
        //		affine2dToTf(t, tran);
        //		return tran;
        //	}

        /**
         * @brief converts from tf transform to 2d affine matrix
         * @param tran transform in tf namespace
         * @param t 2d affine matrix
         */
        //	void TfToAffine2d(const tf::Transform& tran, Transformation2d& t) {
        //		t = Transformation2d::Identity();
        //		t.prerotate(Eigen::Rotation2Dd(tf::getYaw(tran.getRotation())));
        //		t.pretranslate(Eigen::Vector2d(tran.getOrigin().getX(), tran.getOrigin().getY()));
        //	}

        /**
         * @brief converts from tf transform to 2d affine matrix
         * @param tran transform in tf namespace
         * @return 2d affine matrix
         */
        //	Transformation2d TfToAffine2d(const tf::Transform& tran) {
        //		Transformation2d t;
        //		TfToAffine2d(tran, t);
        //		return t;
        //	}

        /**
         * @bried compute difference between two 2d transforms wT1 and wT0 as wT0^-1 * wT1 
         * @param odom1 tranform 1 wrt the world
         * @param odom0 transfom 0 wrt the world
        //     * @param dx difference between transforms in the form [x, y, theta]^T
        //     */
        //	void computeOdometryDiff2d(const nav_msgs::Odometry& odom1, const nav_msgs::Odometry& odom0, Eigen::Vector3d & dx) {
        //		tf::Pose tf1;
        //		tf::poseMsgToTF(odom1.pose.pose, tf1);
        //		Transformation2d wT1;
        //		
        //		tf::Pose tf0;
        //		tf::poseMsgToTF(odom0.pose.pose, tf0);
        //		Transformation2d wT0;
        //		
        //		TfToAffine2d(tf1, wT1);
        //		TfToAffine2d(tf0, wT0);
        //
        //		Transformation2d odomDisp = wT0.inverse() * wT1;
        //		dx(0) = odomDisp.translation()(0);
        //		dx(1) = odomDisp.translation()(1);
        //		dx(2) = atan2(odomDisp.rotation()(1, 0), odomDisp.rotation()(0, 0));
        //	}

        /**
         * @bried compute difference between two 2d transforms wT1 and wT0 as wT0^-1 * wT1 
         * @param odom1 tranform 1 wrt the world
         * @param odom0 transfom 0 wrt the world
         * @return difference between transforms in the form [x, y, theta]^T
         */
        //	Eigen::Vector3d computeOdometryDiff2d(const nav_msgs::Odometry& odom1, const nav_msgs::Odometry& odom0) {
        //		Eigen::Vector3d dx;
        //		computeOdometryDiff2d(odom1, odom0, dx);
        //		return dx;
        //	}
    } // end of TfUtils

} // end of namespace grd

#endif
