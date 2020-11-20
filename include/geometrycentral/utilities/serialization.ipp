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

namespace detail {

// https://stackoverflow.com/a/9154394
template <typename T>
inline auto toSerializedBlob_impl(const T& obj, std::ostringstream& oss, cereal::PortableBinaryOutputArchive& oarchive, std::string& out, int dummy) -> decltype(obj.toSerializedBlob(), void()) {
  out = obj.toSerializedBlob();
}

template <typename T>
inline auto toSerializedBlob_impl(const T& obj, std::ostringstream& oss, cereal::PortableBinaryOutputArchive& oarchive, std::string& out, int long) -> decltype(oarchive(obj), void()) {
  oarchive(obj);
  out = oss.str();
}

}

template <typename T>
inline std::string toSerializedBlob(const T& obj) {
  std::ostringstream oss;
  cereal::PortableBinaryOutputArchive oarchive(oss);
  std::string out;
  detail::toSerializedBlob_impl(obj, oss, oarchive, out, 0);
  return out;
}

template <typename T>
inline void fromSerializedBlob(const std::string& serializedBlob, T& obj) {
  std::istringstream iss(serializedBlob);
  cereal::PortableBinaryInputArchive iarchive(iss);
  iarchive(obj);
}

} // namespace geometrycentral
