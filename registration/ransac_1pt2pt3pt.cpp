#include "registration/ransac_1pt2pt3pt.h"

Matf6D ransac1Pt(Matf6D& x, float t) {
    int s = 1;
    int max_trials = 10000;
    int npts = x.cols();

    float p = 0.99; // Desired probability of choosing at least one samplefree from outliers
    int trialcount = 0;
    int bestscore = 0;
    Matf6D bestinliers; //

    float N = 1; // Dummy initialisation for number of trials.            
    float t2 = 2.0 * t; // 
    float eps = std::numeric_limits<float>::epsilon();
    
    while (N > trialcount) {
        int ind = std::rand() % npts;
        Eigen::Matrix<float, 6, 1> seedpoint = x.col(ind);
        Matf6D lineset = x.colwise() - seedpoint;

        Matf1D D1 = lineset.topRows(3).colwise().norm();
        Matf1D D2 = lineset.bottomRows(3).colwise().norm();
        Matf1D len = (D1 - D2).array().abs();
       
        Mati1D flag = (len.array() < t2).cast<int>();
        Mati1D inlier_column = getNonZeroColumnIndicesFromRowVector(flag);
        Matf6D inliers = x(Eigen::all, inlier_column);

        int s1 = inliers.cols();
        int inlier_size = 0; // 
        for (int i = 1; i <= 50; ++i) {
            Matf3D src = inliers.topRows(3);
            Matf3D dst = inliers.bottomRows(3);

            Eigen::MatrixXf src_dist_matrix(src.cols(), src.cols());
            Eigen::MatrixXf dst_dist_matrix(dst.cols(), dst.cols());
            computeDistanceMatrix(src, src_dist_matrix);
            computeDistanceMatrix(dst, dst_dist_matrix);
            Eigen::MatrixXf Z = (src_dist_matrix - dst_dist_matrix).array().abs();
            Eigen::MatrixXi F = (Z.array() < t2).cast<int>();
            inlier_size = std::ceil(std::sqrt(F.sum()));
            
            Mati1D F_colwise_sum = F.colwise().sum();
            std::vector<int> sorted_column_indices_total;
            sortRowVectorDescending(F_colwise_sum, sorted_column_indices_total);

            std::vector<int> sorted_column_indices_inlier(sorted_column_indices_total.begin(), 
                sorted_column_indices_total.begin() + inlier_size);
            Matf6D selected_inliers = inliers(Eigen::all, sorted_column_indices_inlier);
            inliers = selected_inliers;

            if ((s1 - inlier_size) < 5) {
                break;
            }
            s1 = inlier_size;
        }

        if (inlier_size > bestscore) {
            bestscore = inlier_size;
            bestinliers = inliers;
            float fracinliers = static_cast<float>(inlier_size) / npts;
            float pNoOutliers = 1 - std::pow(fracinliers, s);
            pNoOutliers = std::max(eps, pNoOutliers); // Avoid division by -Inf
            pNoOutliers = std::min(1 - eps, pNoOutliers); // Avoid division by 0.
            N = log(1-p)/log(pNoOutliers);
            N = std::max(N, static_cast<float>(RansacN1)); // at least try 
        }
        ++trialcount; 

        if (trialcount > max_trials) {
            break;
        }
    }
  
    return bestinliers;
}

Mati1D getNonZeroColumnIndicesFromRowVector(const Mati1D& flags) {
    int count = 0;
    Mati1D nonzero_column(1, flags.count());
    for (int i = 0; i < flags.cols(); ++i) {
        if (flags(0, i) > 0) {
            nonzero_column(0, count++) = i;
        } 
    }
    return nonzero_column;
}

void computeDistanceMatrix(const Matf3D& data, Eigen::MatrixXf& dist_matrix) {
    for (int i = 0; i < data.cols(); ++i)
        dist_matrix.row(i) = (data.colwise() - data.col(i)).colwise().norm();
}

void sortRowVectorDescending(const Mati1D& data, std::vector<int>& sorted_column_indices) {
    std::vector<int> sorted_data(data.cols());
    sorted_column_indices.resize(data.cols());
    for (int i = 0; i < data.cols(); ++i) {
        sorted_data[i] = data(0, i);
        sorted_column_indices[i] = i;
    }

    std::sort(sorted_column_indices.begin(), sorted_column_indices.end(), [&sorted_data](int i, int j) {
        return sorted_data[i] > sorted_data[j]; });  
}


Matf6D ransac2Pt(Matf6D& x, float t) {
    int s = 2;
    int max_trials = 10000;
    int npts = x.cols();

    float p = 0.99; // Desired probability of choosing at least one samplefree from outliers
    int trialcount = 0;
    int bestscore = 0;
    Matf6D bestinliers;

    float N = 1; // Dummy initialisation for number of trials.         
    float t2 = 2.0 * t; //
    float eps = std::numeric_limits<float>::epsilon();

    while (N > trialcount) {
        int id_1 = std::rand() % npts;
        int id_2 = std::rand() % npts;
        if (id_1 == id_2) {
            continue;
        }
        Eigen::Matrix<float, 6, 1> diff = x.col(id_1) - x.col(id_2);
        float diff_dist = std::abs(diff.head(3).norm() - diff.tail(3).norm());
        
        if (diff_dist > t2) {
            continue;
        }

        Matf6D xinliers = x;
        Matf6D lineset1 = x.colwise() - x.col(id_1);
        Matf6D lineset2 = x.colwise() - x.col(id_2);
        Matf1D len1 = (lineset1.topRows(3).colwise().norm() - lineset1.bottomRows(3).colwise().norm()).cwiseAbs();
        Matf1D len2 = (lineset2.topRows(3).colwise().norm() - lineset2.bottomRows(3).colwise().norm()).cwiseAbs();;

        Mati1D flag1 = (len1.array() < t2).cast<int>();
        Mati1D flag2 = (len2.array() < t2).cast<int>();
        Mati1D flag = flag1.cwiseMin(flag2);
        Mati1D inlier_column_rough = getNonZeroColumnIndicesFromRowVector(flag);
        
        if (flag.sum() < 3) {
            continue;
        }

        Matf6D lineset1_f = lineset1(Eigen::all, inlier_column_rough);
        Matf6D lineset2_f = lineset2(Eigen::all, inlier_column_rough);
        Matf6D xinliers_temp = xinliers(Eigen::all, inlier_column_rough);//
        xinliers = xinliers_temp;
        Mati1D inlier_flag = Mati1D::Zero(1, lineset1_f.cols());

        for (int i = 0; i < lineset1_f.cols(); ++i) {
            Eigen::Matrix<float, 3, 1> xk = lineset1_f.col(i).head(3);
            Eigen::Matrix<float, 3, 1> yk = lineset1_f.col(i).tail(3);
            Eigen::Matrix<float, 3, 1> xl = lineset2_f.col(i).head(3);
            Eigen::Matrix<float, 3, 1> yl = lineset2_f.col(i).tail(3);
            
            float ty = yk.dot(yl) / (yk.norm() * yl.norm() + 0.00001);
            if (ty < -1) {
                ty = -1;
            } else if (ty > 1) {
                ty = 1;
            } 

            float tx = xk.dot(xl) / (xk.norm() * xl.norm() + 0.00001);
            if (tx < -1) {
                tx = -1;
            } else if (tx > 1) {
                tx = 1;
            } 

            float thetay = std::acos(ty);
            float thetax = std::acos(tx);
            float alpha = std::asin(std::min(0.9999, 2 * t / (xk.norm() + 0.00001))); // 4
            float beta = std::asin(std::min(0.9999, 2 * t / (xl.norm() + 0.00001))); // 4
            float dtheta = std::abs(thetay - thetax);

            if (dtheta < alpha + beta) {
                inlier_flag(0, i) = 1;
            }
        }

        // Find the number of inliers to this model.
        Mati1D inlier_column_refined = getNonZeroColumnIndicesFromRowVector(inlier_flag);
        Matf6D updated_inliers = xinliers(Eigen::all, inlier_column_refined);
        Eigen::MatrixXf src_dist_matrix(updated_inliers.cols(), updated_inliers.cols());
        Eigen::MatrixXf dst_dist_matrix(updated_inliers.cols(), updated_inliers.cols());
        Matf3D src = updated_inliers.topRows(3);
        Matf3D dst = updated_inliers.bottomRows(3);
        computeDistanceMatrix(src, src_dist_matrix);
        computeDistanceMatrix(dst, dst_dist_matrix);

        Eigen::MatrixXi F = ((src_dist_matrix - dst_dist_matrix).array().abs() <= t2).cast<int>();
        int inlier_size = std::ceil(std::sqrt(F.sum()));

        if (inlier_size >= bestscore) {
            bestscore = inlier_size;
            bestinliers = xinliers(Eigen::all, inlier_column_refined);
            float fracinliers = (float)inlier_size / npts;
            float pNoOutliers = 1 - std::pow(fracinliers, s);
            pNoOutliers = std::max(eps, pNoOutliers); //Avoid division by -Inf
            pNoOutliers = std::min(1 - eps, pNoOutliers); //Avoid division by 0.
            N = log(1-p)/log(pNoOutliers);
            N = std::max(N, (float)RansacN2); // at least try RansacN2 times
        }
        ++trialcount; 

        if (trialcount > max_trials) {
            break;
        } 
    }
    return bestinliers;
}

/////////////////////////////////////////////////////////////
Eigen::Matrix4f ransac3Pt(const Matf6D& x, int s, float t) {
    int max_data_trials = 1000;
    int max_trials = 10000;
    int npts = x.cols();
    
    float p = 0.99; //Desired probability of choosing at least one sample
    int trialcount = 0;
    int bestscore = 0;
    Mati1D bestinliers; //

    float N = 1; //Dummy initialisation for number of trials.
    float eps = std::numeric_limits<float>::epsilon();
    
    Eigen::Matrix4f trans = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f trans_best = Eigen::Matrix4f::Identity();
    while (N > trialcount) {
        int count = 1;
        bool degenerate = true;

        while (degenerate) {
            int id_1 = std::rand() % npts;
            int id_2 = std::rand() % npts;
            int id_3 = std::rand() % npts;
            if (id_1 == id_2 || id_1 == id_3 || id_2 == id_3) {
                continue;
            }

            Eigen::Vector3i ind(id_1, id_2, id_3);
            degenerate = isDegenerate(x(Eigen::all, ind));

            if (!degenerate) {
                trans = rigidMotion(x(Eigen::seq(0, 2), ind), x(Eigen::seq(3, 5), ind));
            } 

            ++count;
            if (count > max_data_trials) {
                break;
            }
        }

        Mati1D inlier_column = dist3d(trans, x, t);
        int inlier_size = inlier_column.cols();

        if (inlier_size >= bestscore) {
            bestscore = inlier_size;
            bestinliers = inlier_column;
            trans_best = trans;
            
            float fracinliers = 0, pNoOutliers = 0;  

            Mati1D loInliers = bestinliers;
            int loIter = 0;
            int loRansacMaxIter = 50;//
            int NOSample = std::min(s * 7, (int)loInliers.cols());

            if (inlier_size < s) {
                fracinliers = static_cast<float>(inlier_size) / npts;
                pNoOutliers = 1 - std::pow(fracinliers, s);
                pNoOutliers = std::max(eps, pNoOutliers); //Avoid division by -Inf
                pNoOutliers = std::min(1 - eps, pNoOutliers); //Avoid division by 0.
                N = log(1-p)/log(pNoOutliers);
                continue;                    
            }

            if (trialcount < 50) {
                loIter = loRansacMaxIter;
            }

            std::vector<int> loInliers_vec(&loInliers(0, 0), loInliers.data() + loInliers.size());
            while (loIter < loRansacMaxIter) {
                loIter = loIter + 1;
                std::random_shuffle(loInliers_vec.begin(), loInliers_vec.end());
                std::vector<int> loind(loInliers_vec.begin(), loInliers_vec.begin() + NOSample);
                trans = rigidMotion(x(Eigen::seq(0, 2), loind), x(Eigen::seq(3, 5), loind));
                Mati1D loUpdatedInliers = dist3d(trans, x, t);
                
                if (loUpdatedInliers.cols() > bestscore) {
                    bestscore = loUpdatedInliers.cols();
                    inlier_size = bestscore;
                    bestinliers = loUpdatedInliers;
                    trans_best = trans;
                }
            }
            
            fracinliers =  inlier_size / static_cast<float>(npts);
            pNoOutliers = 1.0 - std::pow(fracinliers, s);
            pNoOutliers = std::max(eps, pNoOutliers);
            pNoOutliers = std::min(1 - eps, pNoOutliers);
            N = log(1-p)/log(pNoOutliers);
            N = std::max(N, (float)RansacN3); // at least try
        }
        ++trialcount;

        if (trialcount > max_trials) {
            break; // Safeguard against being stuck in this loop forever
        }
    }

    if (trans_best.array().isNaN().any()) {
        std::cout << "Three-POint RANSAC's aatrix exist NaN.\n";
        trans_best.setIdentity();
    }
    return trans_best;
}

bool isDegenerate(const Eigen::Matrix<float, 6, 3>& x) {
    Eigen::Matrix<float, 3, 3> x1 = x.topRows(3);
    Eigen::Matrix<float, 3, 3> x2 = x.bottomRows(3);
    bool flag_1 = iscolinear(x1.col(0), x1.col(1), x1.col(2));
    bool flag_2 = iscolinear(x2.col(0), x2.col(1), x2.col(2));
    if (flag_1 || flag_2) {
        return true;
    } else {
        return false;
    }
}
bool iscolinear(const Eigen::Vector3f& p1, const Eigen::Vector3f& p2, 
    const Eigen::Vector3f& p3) {
    Eigen::Vector3f p1p2 = p2 - p1;
    Eigen::Vector3f p1p3 = p3 - p1;
    return (p1p2.cross(p1p3)).norm() < std::numeric_limits<float>::epsilon();
}

Eigen::Matrix4f rigidMotion(const Matf3D& A, const Matf3D& B) {
    Eigen::Matrix<float, 3, 1> lc = A.rowwise().mean();
    Eigen::Matrix<float, 3, 1> rc = B.rowwise().mean();
    Eigen::Matrix3f M = (A.colwise() - lc) * (B.colwise() - rc).transpose();
    
    float Sxx = M(0,0); float Syx = M(1,0); float Szx = M(2,0);
    float Sxy = M(0,1); float Syy = M(1,1); float Szy = M(2,1);
    float Sxz = M(0,2); float Syz = M(1,2); float Szz = M(2,2);
    Eigen::Matrix4f N = Eigen::Matrix4f::Identity();
    N << Sxx+Syy+Szz, Syz-Szy, Szx-Sxz, Sxy-Syx,
        Syz-Szy, Sxx-Syy-Szz, Sxy+Syx, Szx+Sxz,
        Szx-Sxz, Sxy+Syx, -Sxx+Syy-Szz, Syz+Szy,
        Sxy-Syx, Szx+Sxz, Syz+Szy, -Sxx-Syy+Szz;

    // eigen values and vectors
    Eigen::EigenSolver<Eigen::Matrix4f> es(N);
	Eigen::Vector4f evalue = es.eigenvalues().real();
	Eigen::Matrix4f evector = es.eigenvectors().real();
    
    // max evalue and vector
    int max_row = 0, max_col = 0; 
    evalue.maxCoeff(&max_row, &max_col);
    Eigen::Vector4f q = evector.col(max_row);

    q.cwiseAbs().maxCoeff(&max_row, &max_col); 
    float max_vec = q(max_row, 0);
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


Mati1D dist3d(const Eigen::Matrix4f& trans, const Matf6D& x, const float t) {
    Eigen::Matrix3f R = trans.block<3, 3>(0, 0);
    Eigen::Vector3f T = trans.block<3, 1>(0, 3);
    Matf3D x2_ = (R * x.topRows(3)).colwise() + T;
    Matf1D d2 = (x2_ - x.bottomRows(3)).colwise().norm();
    Mati1D flag = (d2.array() < t).cast<int>();
    return getNonZeroColumnIndicesFromRowVector(flag);
}