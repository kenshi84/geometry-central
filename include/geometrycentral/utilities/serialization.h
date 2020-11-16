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

namespace detail {

// https://stackoverflow.com/a/257382
template <typename T>
class has_toSerializedBlob {
  typedef char one;
  struct two { char x[2]; };

  template <typename C> static one test( decltype(&C::toSerializedBlob) ) ;
  template <typename C> static two test(...);    

public:
  enum { value = sizeof(test<T>(0)) == sizeof(char) };
};

}

template <typename T>
typename std::enable_if<detail::has_toSerializedBlob<T>::value, std::string>::type toSerializedBlob(const T& obj);

template <typename T>
typename std::enable_if<!detail::has_toSerializedBlob<T>::value, std::string>::type toSerializedBlob(const T& obj);

template <typename T>
void fromSerializedBlob(const std::string& serializedBlob, T& obj);

} // namespace geometrycentral

#include "geometrycentral/utilities/serialization.ipp"
