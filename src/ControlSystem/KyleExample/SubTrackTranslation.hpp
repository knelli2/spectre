// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "ControlSystem/Actions/ReduceToRunCallbacks.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/Protocols/Submeasurement.hpp"
#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/DefiniteIntegral.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Reduction.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

#include <ostream>
#include <vector>
#include "Parallel/Printf.hpp"

namespace control_system::Actions {
template <typename ControlSystems>
struct ReduceToRunCallbacks;
}  // namespace control_system::Actions

namespace Frame {
struct Inertial;
struct Logical;
}  // namespace Frame

namespace control_system {
namespace QueueTags {
// tag goes in LinkedMessageQueue of all control systems that need this
// submeasurement
struct MeasureTranslation {
  using type = std::vector<double>;
};
}  // namespace QueueTags

struct SubTrackTranslation : tt::ConformsTo<protocols::Submeasurement> {
  template <typename ControlSystem>
  using interpolation_target_tag = void;

  struct VectorAppend {
    std::vector<DataVector> operator()(const std::vector<DataVector>& a,
                                       const std::vector<DataVector>& b) {
      std::vector<DataVector> result{};
      for (auto& dv : a) {
        result.push_back(dv);
      }
      for (auto& dv : b) {
        result.push_back(dv);
      }
      return result;
    }
  };

  using VolumeIntegralDatum =
      Parallel::ReductionDatum<std::vector<double>, funcl::VectorPlus>;

  using ReductionData =
      tmpl::wrap<tmpl::list<Parallel::ReductionDatum<LinkedMessageId<double>,
                                                     funcl::AssertEqual<>>,
                            VolumeIntegralDatum,
                            Parallel::ReductionDatum<std::vector<DataVector>,
                                                     VectorAppend>>,
                 Parallel::ReductionData>;

  // using argument_tags = tmpl::list<::domain::Tags::MeshVelocity<2>>;
  static constexpr size_t Dim = 1;

  using argument_tags = tmpl::list<
      ::domain::Tags::Mesh<Dim>,
      ::domain::Tags::DetInvJacobian<Frame::Logical, Frame::Inertial>,
      ::ScalarWave::Psi, ::domain::Tags::Coordinates<Dim, Frame::Inertial>>;

  template <typename Metavariables, typename ArrayIndex,
            typename ParallelComponent, typename ControlSystems>
  static void apply(
      const typename ::domain::Tags::Mesh<Dim>::type& mesh,
      const Scalar<DataVector>& det_inv_jacobian, const Scalar<DataVector>& psi,
      const typename ::domain::Tags::Coordinates<Dim, Frame::Inertial>::type&
          coords,
      const LinkedMessageId<double>& measurement_id,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ParallelComponent* /*meta*/,
      ControlSystems /*meta*/) {
    // Do reduction here. RunCallbacks is put in the actual reduction
    // Just need a dummy component to do this reduction which doesn't have to be
    // a reduction.
    const auto& control_component_proxy = Parallel::get_parallel_component<
        ControlComponent<Metavariables, tmpl::front<ControlSystems>>>(cache);
    const auto& component_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);
    const auto& self_proxy = component_proxy[array_index];

    const DataVector det_jacobian = 1.0 / get(det_inv_jacobian);
    // const double local_volume = definite_integral(det_jacobian, mesh);
    std::vector<double> integrals{};
    // just integral of psi squared
    integrals.push_back(
        definite_integral(det_jacobian * square(get(psi)), mesh));
    // integrals of psi squared weighted by coordinates to get a "center of
    // mass" kind of thing
    for (size_t i = 0; i < Dim; i++) {
      integrals.push_back(definite_integral(
          det_jacobian * square(get(psi)) * coords.get(i), mesh));
    }

    std::vector<DataVector> vec_coords{coords.get(0)};

    Parallel::contribute_to_reduction<
        Actions::ReduceToRunCallbacks<ControlSystems>>(
        ReductionData(measurement_id, integrals, vec_coords), self_proxy,
        control_component_proxy);
  }
};
}  // namespace control_system
