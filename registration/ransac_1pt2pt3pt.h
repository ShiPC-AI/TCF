#ifndef TCF_REGISTRATION_RANSAC1PT2PT3PT_H_
#define TCF_REGISTRATION_RANSAC1PT2PT3PT_H_
#include "utils/for_cloud.h"

// One-point RANSAC
Matf6D ransac1Pt(Matf6D& x, float t);
Mati1D getNonZeroColumnIndicesFromRowVector(const Mati1D& flags);
void computeDistanceMatrix(const Matf3D& data, Eigen::MatrixXf& dist_matrix);
void sortRowVectorDescending(const Mati1D& data, std::vector<int>& sorted_column_indices);

// Two-point RANSAC
Matf6D ransac2Pt(Matf6D& x, float t);

// Three-point RANSAC
Eigen::Matrix4f ransac3Pt(const Matf6D& x, int s, float t);
bool isDegenerate(const Eigen::Matrix<float, 6, 3>& x);
bool iscolinear(const Eigen::Matrix<float, 3, 1>& p1, const Eigen::Matrix<float, 3, 1>& p2,
    const Eigen::Matrix<float, 3, 1>& p3);
Eigen::Matrix4f rigidMotion(const Matf3D& A, const Matf3D& B);
Mati1D dist3d(const Eigen::Matrix4f& trans, const Matf6D& x, const float t);
#endif