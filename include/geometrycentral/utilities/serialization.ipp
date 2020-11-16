#pragma once

#include <sstream>
#include <cereal/archives/portable_binary.hpp>

namespace cereal {

template<class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void save(Archive& archive, const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& m) {
    archive(cereal::make_nvp("rows", m.rows()));
    archive(cereal::make_nvp("cols", m.cols()));
    for (int i=0; i<m.rows(); ++i)
    for (int j=0; j<m.cols(); ++j)
        archive(cereal::make_nvp("row_" + std::to_string(i) + "_col_" + std::to_string(j), m(i,j)));
}

template<class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void load(Archive& archive, Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& m) {
    Eigen::Index rows, cols;
    archive(cereal::make_nvp("rows", rows));
    archive(cereal::make_nvp("cols", cols));
    m.resize(rows,cols);
    for (int i=0; i<m.rows(); ++i)
    for (int j=0; j<m.cols(); ++j)
        archive(m(i,j));
}

} // namespace cereal

namespace geometrycentral {

template <typename T>
inline std::string toSerializedBlob(const T& obj) {
  std::ostringstream oss;
  cereal::PortableBinaryOutputArchive oarchive(oss);
  oarchive(obj);
  return oss.str();
}

template <typename T>
inline void fromSerializedBlob(const std::string& serializedBlob, T& obj) {
  std::istringstream iss(serializedBlob);
  cereal::PortableBinaryInputArchive iarchive(iss);
  iarchive(obj);
}

} // namespace geometrycentral
