// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <pup.h>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/IdPair.hpp"
#include "DataStructures/LinkedMessageId.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/Callback.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Destination.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"

/// \cond
namespace Frame {
struct BlockLogical;
struct NoFrame;
}  // namespace Frame
/// \endcond

/*!
 * \brief Namespace for data structures that hold all the data of the horizon
 * finder.
 *
 * \details This includes volume data, interpolated variables, the surface of
 * the horizon, and more.
 */
namespace ah::Storage {
// NOTE: Could be a std::pair, but the variable names are more descriptive here
// and we need to define a hash anyways
struct NumberAndId {
  size_t number{0_st};
  LinkedMessageId<double> id;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

bool operator==(const NumberAndId& lhs, const NumberAndId& rhs);
bool operator<(const NumberAndId& lhs, const NumberAndId& rhs);

struct VolumeVars {
  /// The mesh that corresponds to the volume vars
  Mesh<3> mesh;

  /// A `Variables` of the volume tensors that are required to find a horizon.
  /// See `ah::volume_vars_for_horizon_finding` for the list of tensors.
  Variables<ah::volume_vars_for_horizon_finding<3, Frame::NoFrame>>
      volume_vars_for_horizon_finding;

  /// Any extra volume tensors that will need to be interpolated onto the final
  /// horizon surface. These are stored in a `DataVector`, so the Event that
  /// sends them and the post horizon find callback that receives them must know
  /// what order the tensors are in.
  DataVector extra_volume_vars;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

struct VolumeDataAndCallbacks {
  std::string frame;
  Destination destination;
  bool error_on_failure;
  std::unordered_map<ElementId<3>, VolumeVars> volume_vars_per_element;
  std::vector<std::unique_ptr<Parallel::Callback>> callbacks;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

// NOTE: The idea here is to actually separate out observations and control
// system stuff. The control systems are already ordered with the
// LinkedMessageId so I don't think they need to have this extra int we keep
// track of in the Events. The observations, on the other hand, do need this
// extra int, unless we change it so that they are also ordered with a
// LinkedMessageId.
//
// As for the variables. It'd be really easy if we could interpolate the same
// set of variables to both Obs and CS surfaces and then use compute tags to get
// the rest of what we want. CS we know exactly what tags/quantities we need so
// we may even be able to hard code that. However, for Obs we want to choose
// what quantities to observe at so runtime so we'll probably want compute tags.
// But then what to do about frames? For observation do we just interpolate in
// the inertial frame and then do frame transformations later? Does it matter
// since it's observations and doesn't slow down the algorithm?

// NOTE: 9/21/23 morning meeting: For control systems, always store data in
// Frame::NoFrame, for both the vars required for horizon finding (ing_g
// spatial, christoffel spatial, ext curve) and the size control vars that need
// to be interpolated onto the surface. Then, upon sending data to control
// system (or upon receiving it) put it into a Variables in the proper frame.
//
// For observations this is a bit harder because the frame we find the horizon
// in is not the frame of the quantities we are interpolating.
//
// NOTE: For observation, have event send along function for doing observing
// stuff. Then just send sequence of doubles, and then function knows what
// tensors the doubles corresponds to.

struct InterpolatedVars {
  /// `vars` holds the interpolated `Variables` on some subset of the
  /// points in `block_coord_holders`.  The grid points inside vars
  /// are indexed according to `global_offsets` below.  The size of
  /// `vars` changes as more `Element`s send data to this `Interpolator`.
  Variables<ah::volume_vars_for_horizon_finding<3, Frame::NoFrame>> vars{};
  /// Holds the `ElementId`s of `Element`s for which interpolation has
  /// already been done for this `Info`.
  std::unordered_set<ElementId<3>> interpolation_is_done_for_these_elements{};
  std::set<size_t> indicies_interpolated_to_thus_far{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

struct SurfaceAndPoints {
  ylm::Strahlkorper<Frame::NoFrame> strahlkorper{};
  std::string frame{};
  tnsr::I<DataVector, 3, Frame::NoFrame> cartesian_coords{};
  std::vector<std::optional<
      IdPair<domain::BlockId, tnsr::I<double, 3, ::Frame::BlockLogical>>>>
      block_logical_coords{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};
}  // namespace ah::Storage

namespace std {
template <>
struct hash<ah::Storage::NumberAndId> {
  size_t operator()(const ah::Storage::NumberAndId& number_and_id) const {
    return std::hash<size_t>{}(number_and_id.id);
  }
};
}  // namespace std
