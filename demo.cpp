#include "registration/tcf.h"
#include "utils/for_cloud.h"
#include "utils/for_io.h"
#include "utils/for_time.h"
#include <pcl/io/pcd_io.h>
#include <nlohmann/json.hpp>  // for reading json file

using json = nlohmann::json;
int main(int argc, char** argv) {
    std::cout << "==========================================\n";
    std::cout << "======== Demo of TCF Registration ========\n";
    std::cout << "==========================================\n";

    // Open the JSON file
    std::ifstream ifs("../config/config_eth.json"); // eth 
    // std::ifstream ifs("../config/config_kitti.json"); // kitti 
    if (!ifs.is_open()) {
        std::cout << "Cannot open config.json file.\n";
        return 1;
    }

    // Parse the JSON file
    json config;
    ifs >> config;
    std::string path_source_cloud = config["path_source_cloud"];
    std::string path_target_cloud = config["path_target_cloud"];
    std::string path_matches = config["path_matches"];
    std::string path_gt = config["path_gt"];

    // Load ground-truth pose
    Eigen::Matrix4f gt = Eigen::Matrix4f::Identity(); 
    loadMatrix44(path_gt, gt);
    std::cout << "GT pose: \n" << gt << "\n";

    // Load point cloud and resolution
    // this can be replaced by a user-defined value
    CloudPtr source_cloud(new PointCloud), target_cloud(new PointCloud);
    pcl::io::loadPCDFile(path_source_cloud, *source_cloud);
    pcl::io::loadPCDFile(path_target_cloud, *target_cloud);
    float rs = pcResolution(source_cloud);
    float rt = pcResolution(target_cloud);
    float th = std::max(rs, rt);
    std::cout << "Resolution: " << th << " m\n";

    // Load correspondences
    Eigen::MatrixXf matches;
    loadMatrixDynamic(path_matches, matches);
    MatfD3 source_match = matches.leftCols(3);
    MatfD3 target_match = matches.rightCols(3);

    // Start registration
    TicToc tic_tcf;
    std::srand(unsigned(std::time(nullptr)));
    Eigen::Matrix4f trans = twoStageConsensusFilter(source_match, target_match, 3*th);
    double time_registration = tic_tcf.toc();
    std::cout << "Runtime: " << time_registration << " ms.\n";

    // compute error
    std::pair<double, double> error = computeTransError(trans, gt);
    std::cout << "RE: " << error.first << " deg, TE: " << error.second << " m.\n";
    
    return 0;
}