#include "registration/tcf.h"

Eigen::Matrix4f twoStageConsensusFilter(MatfD3 match_1, MatfD3 match_2, float t) {
    Eigen::Matrix4f trans = Eigen::Matrix4f::Identity();
    int source_num = match_1.rows();
    int target_num = match_2.rows();

    if (source_num != target_num) {
        PCL_ERROR("Correspondence must have the same dimension.\n"); 
        return trans;
    }
    if (source_num < 3) {
        PCL_WARN("Must have at least 3 points to fitTING.\n"); 
        return trans;
    }

    //  two point cloud are 'stacked' to create a 6xN array for ransac
    Matf6D stacked_mat(6, source_num);
    stacked_mat.topRows(3) = match_1.transpose();
    stacked_mat.bottomRows(3) = match_2.transpose();

    // ONe-point RANSAC
    Matf6D xinliers_1 = ransac1Pt(stacked_mat, t);
    if (xinliers_1.cols() < 3) {
        PCL_WARN("ransac 1 matches less than 3\n");
        return trans;
    }
   
    // Two-point RANSAC
    Matf6D xinliers_2 = ransac2Pt(xinliers_1, t);
    if (xinliers_2.cols() < 3) {
        PCL_WARN("ransac 2 matches less than 3\n");
        return trans;
    }

    // Three-point RANSAC
    trans = ransac3Pt(xinliers_2, 3, t); 
    Mati1D inlier_column = dist3d(trans, xinliers_2, t);
    Matf6D xinliers_3 = xinliers_2(Eigen::all, inlier_column);
    if (xinliers_3.cols() < 3) {
        PCL_WARN("ransac 3 matches less than 3\n");
        return trans;
    }

    // SA-Cauchy IRLS;
    trans = saCauchyIRLSRigidModel(xinliers_3.topRows(3), xinliers_3.bottomRows(3), 1.3);
    return trans;
}