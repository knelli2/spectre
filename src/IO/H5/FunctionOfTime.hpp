// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <hdf5.h>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/Object.hpp"
#include "IO/H5/OpenGroup.hpp"

namespace h5::detail {
class FunctionOfTime : public h5::Object {
 public:
  /// \cond HIDDEN_SYMBOLS
  static std::string extension() { return ".fot"; }

  FunctionOfTime(bool exists, detail::OpenGroup&& group, hid_t location,
                 const std::string& name,
                 const std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>&
                     initial_function_of_time = nullptr,
                 uint32_t version = 0);

  FunctionOfTime(const FunctionOfTime& /*rhs*/) = delete;
  FunctionOfTime& operator=(const FunctionOfTime& /*rhs*/) = delete;
  FunctionOfTime(FunctionOfTime&& /*rhs*/) = delete;             // NOLINT
  FunctionOfTime& operator=(FunctionOfTime&& /*rhs*/) = delete;  // NOLINT

  ~FunctionOfTime() override;
  /// \endcond HIDDEN_SYMBOLS

  /*!
   * Appends the \p time and \p update to the `Updates.dat` subfile
   */
  void append(double time, const DataVector& update);

  /*!
   * \returns the legend of the FunctionOfTime file
   */
  const std::vector<std::string>& get_legend() const { return legend_; }

  /*!
   * \returns all data stored in the FunctionOfTime file
   *
   * \example
   * \snippet Test_FunctionOfTime.cpp h5dat_get_data
   */
  template <typename T = Matrix>
  T get_data() const;

  /*!
   * \brief Get only some columns over a range of rows
   * \requires all members of `these_columns` have a value less than the number
   * of columns, `first_row < last_row` and `last_row` is less than or equal to
   * the number of rows
   * \returns a subset of the data from the FunctionOfTime file
   *
   * \example
   * \snippet Test_FunctionOfTime.cpp h5dat_get_subset
   */
  template <typename T = Matrix>
  T get_data_subset(const std::vector<size_t>& these_columns,
                    size_t first_row = 0, size_t num_rows = 1) const;

  /*!
   * \returns the number of rows (first index) and columns (second index)
   */
  const std::array<hsize_t, 2>& get_dimensions() const { return size_; }

  /*!
   * \returns the header of the FunctionOfTime file
   */
  const std::string& get_header() const { return header_; }

  /*!
   * \returns the user-specified version number of the FunctionOfTime file
   *
   * \note h5::Version returns a uint32_t, so we return one here too for the
   * version
   */
  uint32_t get_version() const { return version_; }

  const std::string& subfile_path() const override { return path_; }

 private:
  /// \cond HIDDEN_SYMBOLS
  detail::OpenGroup group_;
  std::string name_;
  std::string path_;
  uint32_t version_;
  std::optional<h5::Dat> updates_;
  std::vector<std::string> legend_;
  std::string header_;
  hid_t dataset_id_{-1};
  /// \endcond HIDDEN_SYMBOLS
};
}  // namespace h5::detail
