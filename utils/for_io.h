#ifndef TCF_UTILS_IO_H_
#define TCF_UTILS_IO_H_

#include <iomanip> // for setprecision
#include <cstdlib> // for EXIT_FAILURE
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

template <typename T>
void saveMatrixDynamic(const std::string& filename, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, const int& N) {
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::exit((std::cout << "Failed to save:" << filename << "\n", EXIT_FAILURE));
    ofs << std::setprecision(N); 
    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) 
            ofs << mat(i, j) << " ";
        ofs << "\n";
    }
    ofs.close();
}

template <typename T>
void loadMatrixDynamic(const std::string& filename, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) 
        std::exit((std::cout << "Failed to open:" << filename << "\n", EXIT_FAILURE));

    std::vector<std::vector<T>> values;
    std::string line;
    while (std::getline(ifs, line)) {
        std::istringstream ss(line);
        values.emplace_back(std::istream_iterator<T>{ss}, std::istream_iterator<T>{}); // C++ >=17
    }

    if (values.empty() || values[0].empty()) 
        std::exit((std::cout << "Matrix data is empty or malformed.", EXIT_FAILURE));
    mat.resize(values.size(), values[0].size());
    for (int i = 0; i < values.size(); ++i) 
        for (int j = 0; j < values[0].size(); ++j) 
            mat(i, j) = values[i][j];
    ifs.close();
}

template <typename T>
void loadMatrix44(const std::string& filename, Eigen::Matrix<T, 4, 4>& matrix) {
    matrix.setIdentity();
    std::ifstream ifs(filename);
    if (!ifs.is_open()) 
        std::exit((std::cout << "Failed to open 4*4 matrix file:" << filename << "\n", EXIT_FAILURE));

    ifs >> matrix(0, 0) >> matrix(0, 1) >> matrix(0, 2) >> matrix(0, 3)
        >> matrix(1, 0) >> matrix(1, 1) >> matrix(1, 2) >> matrix(1, 3) 
        >> matrix(2, 0) >> matrix(2, 1) >> matrix(2, 2) >> matrix(2, 3)
        >> matrix(3, 0) >> matrix(3, 1) >> matrix(3, 2) >> matrix(3, 3);
    ifs.close();
};

template <typename T>
void saveMatrix44(const std::string& filename, const Eigen::Matrix<T, 4, 4>& matrix, const int& N) {
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::exit((std::cout << "Failed to save 4*4 matrix file:" << filename << "\n", EXIT_FAILURE));
    ofs << std::setprecision(N);
    ofs << double(matrix(0, 0)) << " " << double(matrix(0, 1)) << " " << double(matrix(0, 2)) << " " << double(matrix(0, 3)) << std::endl
        << double(matrix(1, 0)) << " " << double(matrix(1, 1)) << " " << double(matrix(1, 2)) << " " << double(matrix(1, 3)) << std::endl 
        << double(matrix(2, 0)) << " " << double(matrix(2, 1)) << " " << double(matrix(2, 2)) << " " << double(matrix(2, 3)) << std::endl
        << double(matrix(3, 0)) << " " << double(matrix(3, 1)) << " " << double(matrix(3, 2)) << " " << double(matrix(3, 3));
    ofs.close();
};
#endif