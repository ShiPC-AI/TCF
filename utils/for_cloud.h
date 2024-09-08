#ifndef TCF_UTILS_CLOUD_H_
#define TCF_UTILS_CLOUD_H_

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/random_sample.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <pcl/search/kdtree.h>

// ransac
#define RansacN1 0
#define RansacN2 0
#define RansacN3 0

// point cloud
typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;
typedef pcl::PointCloud<pcl::PointXYZ>::Ptr CloudPtr;

// eigen matrix
// Mat<type><rows><cols>, D:dynamic
typedef Eigen::Matrix<float, 6, Eigen::Dynamic> Matf6D;
typedef Eigen::Matrix<float, 3, Eigen::Dynamic> Matf3D;
typedef Eigen::Matrix<float, Eigen::Dynamic, 3> MatfD3;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic> Mati1D;
typedef Eigen::Matrix<float, 1, Eigen::Dynamic> Matf1D;

// Random downsample the point cloud
void randomSampleCloud(const CloudPtr& cloud_in, CloudPtr& cloud_out,  int N);

// Voxel downsample the point cloud
void voxelSampleCloud(const CloudPtr& cloud_in, CloudPtr& cloud_out, const float leaf);

// Estimate the resolution of the point cloud
float pcResolution(CloudPtr& cloud);

// Dynamic float N*3-type matrix into point cloud 
void matrixN3toCloud(const MatfD3& mat, CloudPtr& cloud);

// Estimate registraiton error
std::pair<double, double> computeTransError(const Eigen::Matrix4f& trans, const Eigen::Matrix4f& gt);

#endif
