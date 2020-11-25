#pragma once

#include <sstream>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>

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

template<class Archive, typename Scalar, typename StorageIndex>
void save(Archive& archive, const Eigen::Triplet<Scalar, StorageIndex>& t) {
  archive(cereal::make_nvp("row", t.row()));
  archive(cereal::make_nvp("col", t.col()));
  archive(cereal::make_nvp("value", t.value()));
}

template<class Archive, typename Scalar, typename StorageIndex>
void load(Archive& archive, Eigen::Triplet<Scalar, StorageIndex>& t) {
  StorageIndex row, col;
  Scalar value;
  archive(row, col, value);
  t = {row, col, value};
}

template<class Archive, typename Scalar, int Options, typename StorageIndex>
void save(Archive& archive, const Eigen::SparseMatrix<Scalar,Options,StorageIndex>& m) {
  Eigen::Index rows = m.rows();
  Eigen::Index cols = m.cols();
  archive(CEREAL_NVP(rows));
  archive(CEREAL_NVP(cols));

  std::vector<Eigen::Triplet<Scalar, StorageIndex>> triplets;
  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<Scalar,Options,StorageIndex>::InnerIterator it(m, k); it; ++it) {
      triplets.emplace_back(it.row(), it.col(), it.value());
    }
  }
  archive(CEREAL_NVP(triplets));
}

template<class Archive, typename Scalar, int Options, typename StorageIndex>
void load(Archive& archive, Eigen::SparseMatrix<Scalar,Options,StorageIndex>& m) {
  Eigen::Index rows, cols;
  archive(cereal::make_nvp("rows", rows));
  archive(cereal::make_nvp("cols", cols));

  std::vector<Eigen::Triplet<Scalar, StorageIndex>> triplets;
  archive(CEREAL_NVP(triplets));

  m = Eigen::SparseMatrix<Scalar,Options,StorageIndex>(rows, cols);
  m.setFromTriplets(triplets.begin(), triplets.end());
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
