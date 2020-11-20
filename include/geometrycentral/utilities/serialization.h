#pragma once

#include <Eigen/Core>
#include <string>
#include <type_traits>

namespace cereal {

template<class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void save(Archive& archive, const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& m);

template<class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void load(Archive& archive, Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& m);

} // namespace cereal

namespace geometrycentral {

template <typename T>
std::string toSerializedBlob(const T& obj);

template <typename T>
void fromSerializedBlob(const std::string& serializedBlob, T& obj);

} // namespace geometrycentral

#include "geometrycentral/utilities/serialization.ipp"
