#include "utils/for_cloud.h"

void randomSampleCloud(const CloudPtr& cloud_in, 
    CloudPtr& cloud_out, int N) {
    pcl::RandomSample<pcl::PointXYZ> random_sampler;
    random_sampler.setInputCloud(cloud_in);
    random_sampler.setSample(N);
    random_sampler.filter(*cloud_out);
}

void voxelSampleCloud(const CloudPtr& cloud_in, CloudPtr& cloud_out, const float leaf) {
    pcl::VoxelGrid<pcl::PointXYZ> grid;
    grid.setInputCloud(cloud_in);
    grid.setLeafSize(leaf, leaf, leaf);
    grid.filter(*cloud_out);
}

float pcResolution(CloudPtr& cloud) {
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
    tree->setInputCloud(cloud);

    // sample size
    int search_num = std::floor((int)cloud->size() / 20) * 10;
    int pt_num = cloud->size();

    std::vector<float> distances(search_num);
    for (int i = 0; i < search_num; ++i) {
        int idx = std::rand() % pt_num;
        std::vector<int> idx_nknsearch;
        std::vector<float> sqdist_nknsearch;
        tree->nearestKSearch(cloud->points[idx], 2, idx_nknsearch, sqdist_nknsearch);
        distances[i] = sqdist_nknsearch[1];
    }
    std::sort(distances.begin(), distances.end());
    int id_mid = (int)(search_num - 1) / 2;
    return std::sqrt(distances[id_mid]);
}

void matrixN3toCloud(const MatfD3& mat, CloudPtr& cloud) {
    cloud->resize(mat.rows());
    for (int i = 0; i < mat.rows(); ++i) {
        cloud->points[i].x = mat(i, 0);
        cloud->points[i].y = mat(i, 1);
        cloud->points[i].z = mat(i, 2);
        std::cout << cloud->points[i].getVector3fMap().transpose() << "\n";
    }
}

std::pair<double, double> computeTransError(const Eigen::Matrix4f& trans, const Eigen::Matrix4f& gt) {
    Eigen::Matrix3f eR = gt.block<3, 3>(0, 0) * trans.block<3, 3>(0, 0).transpose();
    double er_cos = std::max(-1.0, std::min(1.0, (eR.trace() - 1.0)/2.0));
    double er_degree = std::acos(er_cos) * 180.0 / M_PI;
    double et = (trans - gt).block<3, 1>(0, 3).norm();
    return std::make_pair(er_degree, et);
}