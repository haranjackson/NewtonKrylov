#ifndef TYPES_H
#define TYPES_H

#include "include/eigen3/Eigen"


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;

typedef Eigen::Ref<Mat> Matr;
typedef Eigen::Ref<Vec> Vecr;

typedef Eigen::HouseholderQR<Mat> DecQR;

typedef std::function<Vec(Vec)> VecFunc;


#endif
