#include "registration/irls_welsch.h"

Eigen::Matrix4f saCauchyIRLSRigidModel(const Matf3D& src, const Matf3D& dst, const float& tau) {
    float prev_cost = std::pow(10, 15); // initial energy cost
    int maxIter = 100; // maximum iteration times
    int n = src.cols(); // measurement num
    Eigen::MatrixXf weights = Eigen::MatrixXf::Ones(1, n); // weights vector
    float alpha = 0.0; 
    
    Mati1D inlier_column = Eigen::VectorXi::LinSpaced(n, 0, n - 1).transpose(); // inlier col indices
    Eigen::Matrix4f trans = Eigen::Matrix4f::Identity(); // output pose

    Matf3D src_current = src;
    Matf3D dst_current = dst;
    for (int i = 1; i <= maxIter; ++i) {
        Matf3D src_temp = src_current(Eigen::all, inlier_column);
        Matf3D des_temp = dst_current(Eigen::all, inlier_column);
        src_current = src_temp;
        dst_current = des_temp;
          
        trans = rigidTrans(src_current, dst_current, weights); // compute transform
        Matf3D fit = (trans.block<3, 3>(0, 0) * src_current).colwise() + trans.block<3, 1>(0, 3);
        Matf1D residuals = (fit - dst_current).colwise().norm(); // compute error

        if (i == 1) {
            alpha = residuals.cwiseAbs().maxCoeff(); // intial alpha value
        }
        Matf1D residuals2 = residuals.array().square();
        float cost = weights.cwiseProduct(residuals2).sum(); // current cost sum

        Mati1D flags = (residuals.array() < (3 * alpha)).cast<int>();
        inlier_column = getNonZeroColumnIndicesFromRowVector(flags); // update inlier indcies
        
        Matf1D E = residuals(Eigen::all, inlier_column); // update error of matching points
        weights = (E.array().square().array() + std::pow(alpha, 2)).cwiseInverse() * std::pow(alpha, 2); // update SA-Cauch Weights
        float cost_diff = std::abs(cost - prev_cost); // cost diff

        alpha = alpha / tau; // update alpha
        prev_cost = cost; // update cost
        if (cost_diff < 0.01 || alpha < 1.0) {
            break;
        } 
    }

    return trans;
}

Eigen::Matrix4f rigidTrans(const Matf3D& A, const Matf3D& B, Eigen::MatrixXf& weights) {
    float sw = weights.sum();
    if (sw < std::numeric_limits<float>::epsilon()) {
        weights = Eigen::MatrixXf::Ones(1, A.cols());
        sw = weights.sum();
    }

    Eigen::MatrixXf w = weights / sw;
    Eigen::Matrix<float, 3, 1> lc = A * w.transpose();
    Eigen::Matrix<float, 3, 1> rc = B * w.transpose();
    Eigen::MatrixXf w2 = w.cwiseSqrt();
    Eigen::MatrixXf w2_repmat(A.rows(), A.cols());
    for (int i = 0; i < w2_repmat.rows(); ++i) {
        w2_repmat.row(i) = w2;
    }

    Eigen::MatrixXf left = (A.colwise() - lc).cwiseProduct(w2_repmat);
    Eigen::MatrixXf right = (B.colwise() - rc).cwiseProduct(w2_repmat);
    Eigen::MatrixXf M = left * right.transpose();

    float Sxx = M(0,0); float Syx = M(1,0); float Szx = M(2,0);
    float Sxy = M(0,1); float Syy = M(1,1); float Szy = M(2,1);
    float Sxz = M(0,2); float Syz = M(1,2); float Szz = M(2,2);
    Eigen::Matrix4f N = Eigen::Matrix4f::Identity();
    N << Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx,
        Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz,
        Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy,
        Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz;

    // eigen values and vectors
    Eigen::EigenSolver<Eigen::Matrix4f> es(N);
	Eigen::Vector4f evalue = es.eigenvalues().real();
	Eigen::Matrix4f evector = es.eigenvectors().real();
    
    // max evalue and vector
    int maxRow = 0, maxCol = 0; 
    evalue.maxCoeff(&maxRow, &maxCol);
    Eigen::Vector4f q = evector.col(maxRow);

    q.cwiseAbs().maxCoeff(&maxRow, &maxCol); 
    float max_vec = q(maxRow, 0);
    if (max_vec < 0) {
        q = q * (-1);
    }
    q.normalize();
    float q0 = q(0);
    float qx = q(1);
    float qy = q(2);
    float qz = q(3);
    Eigen::Vector3f v = q.tail(3);

    Eigen::Matrix3f Z = Eigen::Matrix3f::Identity();
    Z << q0, -qz, qy,
        qz, q0, -qx,
        -qy, qx, q0;

    Eigen::Matrix3f R = v * v.transpose() + Z * Z;
    Eigen::Vector3f t = rc - R * lc;

    Eigen::Matrix4f T = Eigen::Matrix4f::Identity();
    T.block<3, 3>(0, 0) = R;
    T.block<3, 1>(0, 3) = t;
    return T;
}