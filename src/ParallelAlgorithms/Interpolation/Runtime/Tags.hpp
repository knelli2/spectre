// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <deque>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
template <size_t VolumeDim>
class ElementId;
/// \endcond

namespace intrp2 {

namespace OptionTags {
/*!
 * \ingroup OptionGroupsGroup
 * \brief Groups option tags for InterpolationTargets.
 */
struct InterpolationTargets {
  static constexpr Options::String help{"Options for interpolation targets"};
};

/// Option tag that determines if volume data will be dumped from the
/// Interpolator upon a failure.
struct DumpVolumeDataOnFailure {
  using type = bool;
  static constexpr Options::String help{
      "Whether or not to dump all volume data currently stored by the "
      "interpolator. Volume data is written to the file corresponding to the "
      "node it was collected on."};
  using group = InterpolationTargets;
};
}  // namespace OptionTags

/// Tags for items held in the `DataBox` of `InterpolationTarget` or
/// `Interpolator`.
namespace Tags {
/// Tag that determines if volume data will be dumped form the Interpolator upon
/// failure
struct DumpVolumeDataOnFailure : db::SimpleTag {
  using type = bool;
  using option_tags = tmpl::list<OptionTags::DumpVolumeDataOnFailure>;
  static constexpr bool pass_metavariables = false;

  static bool create_from_options(const bool input) { return input; }
};

/// `temporal_id`s that we have already interpolated onto.
///  This is used to prevent problems with multiple late calls to
///  AddTemporalIdsToInterpolationTarget.
template <typename TemporalId>
struct CompletedTemporalIds : db::SimpleTag {
  using type = std::deque<TemporalId>;
};

///////////////////////////////////////////////////////////////////////////////
// Tags needed for runtime targets
/*!
 * \brief Map between the string of the InterpolationTarget component element
 * and the db::Access corresponding to that element.
 *
 * \note This is the expected way to retrieve and mutate tags on the target
 * array element.
 */
struct DbAccess : db::SimpleTag {
  using type = std::unique_ptr<db::Access>;
};

/*!
 * \brief `std::unordered_set` of `db::tag_name`s of quantities to observe
 */
struct VarsToObserve : db::SimpleTag {
  using type = std::unordered_set<std::string>;
};

/*!
 * \brief Points of the target surface in the frame specified from the input
 * file.
 */
template <size_t Dim>
struct Points : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame::NoFrame>;
};

/*!
 * \brief How many sets of points are stored in each quantity.
 *
 * \details It is possible for a target to specify multiple "sets" of its points
 * in a single contiguous buffer (think multiple radii of a sphere all one after
 * each other). This tag stores how many of these sets there are (i.e. how many
 * radii for the sphere example).
 *
 * This is needed because compute tags typically assume there is just one "set"
 * of points and cannot handle if there are multiple in the same buffer.
 * Therefore, this tag is intended to be used to know how many views to create
 * into the original buffer so we don't have to write special compute tags for
 * the interpolation framework.
 */
struct NumberOfSetsOfPoints : db::SimpleTag {
  using type = size_t;
};

/*!
 * \brief Frame of target points in `intrp2::Tags::Points`.
 */
struct Frame : db::SimpleTag {
  using type = std::string;
};

/*!
 * \brief For each time of an interpolation, this is a map between tensors and
 * the interpolated volume data
 */
template <typename TemporalId>
struct InterpolatedVars : db::SimpleTag {
  using type = std::unordered_map<
      TemporalId, std::unordered_map<std::string, std::vector<DataVector>>>;
};

struct InvalidPointsFillValue : db::SimpleTag {
  using type = double;
};

template <typename TemporalId>
struct NumberOfFilledPoints : db::SimpleTag {
  using type = std::unordered_map<TemporalId, size_t>;
};

template <typename TemporalId>
struct CurrentTemporalIds : db::SimpleTag {
  using type = std::unordered_set<TemporalId>;
};

template <typename Target>
struct Callbacks : db::SimpleTag {
  using type =
      std::vector<std::unique_ptr<intrp2::callbacks::Callback<Target>>>;
};

template <typename Target, size_t Dim>
using common_target_tags =
    tmpl::list<intrp2::Tags::Frame, intrp2::Tags::Points<Dim>,
               intrp2::Tags::VarsToObserve,
               intrp2::Tags::InvalidPointsFillValue>;
}  // namespace Tags
}  // namespace intrp2
