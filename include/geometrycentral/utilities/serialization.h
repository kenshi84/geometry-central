#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>
#include <type_traits>

namespace cereal {

template<class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void save(Archive& archive, const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& m);

template<class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void load(Archive& archive, Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& m);

template<class Archive, typename Scalar, typename StorageIndex>
void save(Archive& archive, const Eigen::Triplet<Scalar, StorageIndex>& t);

template<class Archive, typename Scalar, typename StorageIndex>
void load(Archive& archive, Eigen::Triplet<Scalar, StorageIndex>& t);

template<class Archive, typename Scalar, int Options, typename StorageIndex>
void save(Archive& archive, const Eigen::SparseMatrix<Scalar,Options,StorageIndex>& m);

template<class Archive, typename Scalar, int Options, typename StorageIndex>
void load(Archive& archive, Eigen::SparseMatrix<Scalar,Options,StorageIndex>& m);

} // namespace cereal

namespace geometrycentral {

template <typename T>
std::string toSerializedBlob(const T& obj);

template <typename T>
void fromSerializedBlob(const std::string& serializedBlob, T& obj);

} // namespace geometrycentral

#include "geometrycentral/utilities/serialization.ipp"
