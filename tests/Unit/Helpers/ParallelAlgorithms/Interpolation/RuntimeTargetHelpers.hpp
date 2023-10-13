// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/AlignedLattice.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/Target.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace Frame {
struct Grid;
struct Distorted;
struct Inertial;
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace intrp::Targets::TestHelpers {
template <typename DerivedClass, size_t Dim, typename... ExtraCacheTags>
struct Metavars {
  using const_global_cache_tags =
      tmpl::list<domain::Tags::Domain<Dim>, ExtraCacheTags...>;

  struct factory_creation {
    using factory_classes =
        tmpl::map<tmpl::pair<Target<Dim>, tmpl::list<DerivedClass>>>;
  };

  using component_list = tmpl::list<>;
};

/*!
 * \brief Tests the factory creation of a runtime `intrp::Targets::Target` for a
 * single frame.
 *
 * \details Checks that the target points of the derived class are the expected
 * ones, along with the expected name and frame.
 *
 * \tparam DerivedClass Derived class of `intrp::Targets::Target`
 * \tparam Dim Spatial dimension. Can be deduced from \p expected_points
 * \tparam TargetFrame Frame of points. Can be deduced from \p expected_points
 * \tparam DbTags Tags list for \p box
 * \tparam Metavariables metavariables for \p cache
 * \param option_string Option string with a `Frame:` option passed to
 * `::TestHelpers::test_factory_creation()`
 * \param expected_points Expected points in the proper frame
 * \param expected_name Expected name of the \p DerivedClass
 * \param time The time to pass to
 * `intrp::Targets::Target::block_logical_coordinates()`
 * \param box `db::DataBox`
 * \param cache `Parallel::GlobalCache`
 */
template <typename DerivedClass, size_t Dim, typename TargetFrame,
          typename DbTags, typename Metavariables>
void test_target_single_frame(
    const std::string& option_string,
    const tnsr::I<DataVector, Dim, TargetFrame>& expected_points,
    const std::string& expected_name, const double time,
    const db::DataBox<DbTags>& box,
    const Parallel::GlobalCache<Metavariables>& cache) {
  const std::unique_ptr<Target<Dim>> target =
      ::TestHelpers::test_factory_creation<Target<Dim>, DerivedClass>(
          option_string);

  CHECK(target->name() == expected_name);

  const auto points = target->block_logical_coordinates(
      box, cache, time, get_output(TargetFrame{}));

  for (size_t i = 0; i < points.size(); i++) {
    REQUIRE(points[i].has_value());
    const auto& block_logical_point = points[i].value().data;
    for (size_t d = 0; d < Dim; d++) {
      // We multiply by 0.01 because our domain goes from -100 to 100 so the
      // points are scaled by a factor of 100
      CHECK(block_logical_point.get(d) ==
            approx(0.01 * expected_points.get(d)[i]));
    }
  }
}

/*!
 * \brief Tests the factory creation of a runtime `intrp::Targets::Target` for
 * all frames using `test_target_single_frame()`.
 *
 * \details Loops over `Frame::Grid`, `Frame::Distorted`, `Frame:Inertial` and
 * calls `test_target_single_frame()` for each one. A DataBox is constructed
 * using the passed in \p BoxTags. A Global cache is constructed with a
 * `domain::Tags::Domain` and the \p ExtraCacheTags. The domain that is
 * constructed is an `domain::creators::AlignedLattice` (because it's available
 * for all dimensions) with bounds -100 to 100 in each dimension and no time
 * dependence. This means that all test target points must be within this
 * domain.
 *
 * \tparam DerivedClass Derived class of `intrp::Targets::Target`
 * \tparam Dim Spatial dimension. Can be deduced from \p
 * expected_points_no_frame
 * \tparam BoxTags DataBox tags
 * \tparam ExtraCacheTags Extra DataBox tags for GlobalCache
 * \param option_string Option string to factory create a \p DerivedClass
 * \param expected_points_no_frame Target points but in a tensor of
 * `Frame::NoFrame`. A frame will be added in this function before it is passed
 * to `test_target_single_frame()`.
 * \param expected_name Expected name of the \p DerivedClass
 * \param time The time to pass to
 * `intrp::Targets::Target::block_logical_coordinates()`. Default 1.0
 * \param box_items Items to construct DataBox with. Default empty
 * `tuples::TaggedTuple`
 * \param extra_cache_items Extra items to construct GlobalCache with. The
 * GlobalCache is always constructed with a Domain. Default empty
 * `tuples::TaggedTuple`
 * \return std::unique_ptr<Target<Dim>> A unique pointer
 * to a newly factory created Target<Dim> using the \p option_string
 */
template <typename DerivedClass, size_t Dim, typename... BoxTags,
          typename... ExtraCacheTags>
std::unique_ptr<Target<Dim>> test_target_all_frames(
    const std::string& option_string,
    const tnsr::I<DataVector, Dim, Frame::NoFrame>& expected_points_no_frame,
    const std::string& expected_name, const double time = 1.0,
    tuples::TaggedTuple<BoxTags...> box_items = tuples::TaggedTuple<>{},
    tuples::TaggedTuple<ExtraCacheTags...> extra_cache_items =
        tuples::TaggedTuple<>{}) {
  CAPTURE(expected_name);
  CAPTURE(Dim);
  CAPTURE(option_string);

  const auto box = db::create<tmpl::list<BoxTags...>>(
      std::move(tuples::get<BoxTags>(box_items))...);

  domain::creators::AlignedLattice<Dim> lattice{
      make_array<Dim, std::vector<double>>(std::vector<double>{-100.0, 100.0}),
      make_array<Dim, size_t>(0_st),
      {},
      {},
      {},
      {},
      make_array<Dim, bool>(false)};
  const Parallel::GlobalCache<Metavars<DerivedClass, Dim, ExtraCacheTags...>>
      cache{{lattice.create_domain(),
             std::move(tuples::get<ExtraCacheTags>(extra_cache_items))...}};

  tmpl::for_each<tmpl::list<Frame::Grid, Frame::Distorted, Frame::Inertial>>(
      [&](const auto frame_v) {
        using frame = tmpl::type_from<decltype(frame_v)>;
        const tnsr::I<DataVector, Dim, frame> expected_points_in_frame{};
        for (size_t i = 0; i < Dim; i++) {
          make_const_view(make_not_null(&expected_points_in_frame.get(i)),
                          expected_points_no_frame.get(i), 0,
                          expected_points_no_frame.get(i).size());
        }

        test_target_single_frame<DerivedClass>(option_string,
                                               expected_points_in_frame,
                                               expected_name, time, box, cache);
      });

  return ::TestHelpers::test_factory_creation<Target<Dim>, DerivedClass>(
      option_string);
}
}  // namespace intrp::Targets::TestHelpers
