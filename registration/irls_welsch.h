#ifndef TCF_REGISTRATION_IRLS_WELSCH_H_
#define TCF_REGISTRATION_IRLS_WELSCH_H_
#include "utils/for_cloud.h"
#include "registration/ransac_1pt2pt3pt.h"

Eigen::Matrix4f saCauchyIRLSRigidModel(const Matf3D& src, const Matf3D& dst, const float& tau);
Eigen::Matrix4f rigidTrans(const Matf3D& A, const Matf3D& B, Eigen::MatrixXf& weights);
#endif