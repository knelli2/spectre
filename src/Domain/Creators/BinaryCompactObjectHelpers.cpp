// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/BinaryCompactObjectHelpers.hpp"

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>

#include "DataStructures/DataVector.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Wedge.hpp"
#include "Domain/DomainHelpers.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "Options/ParseError.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

#include "Parallel/Printf.hpp"

namespace domain::creators::bco {
std::unordered_map<std::string, tnsr::I<double, 3, Frame::Grid>>
create_grid_anchors(const std::array<double, 3>& center_a,
                    const std::array<double, 3>& center_b) {
  std::unordered_map<std::string, tnsr::I<double, 3, Frame::Grid>> result{};
  result["Center" + get_output(ObjectLabel::A)] =
      tnsr::I<double, 3, Frame::Grid>{center_a};
  result["Center" + get_output(ObjectLabel::B)] =
      tnsr::I<double, 3, Frame::Grid>{center_b};
  result["Center"] = tnsr::I<double, 3, Frame::Grid>{std::array{0.0, 0.0, 0.0}};

  return result;
}

template <bool IsCylindrical>
TimeDependentMapOptions<IsCylindrical>::TimeDependentMapOptions(
    double initial_time,
    std::optional<ExpansionMapOptions> expansion_map_options,
    std::optional<RotationMapOptions> rotation_options,
    std::optional<ShapeMapOptions<domain::ObjectLabel::A>> shape_options_A,
    std::optional<ShapeMapOptions<domain::ObjectLabel::B>> shape_options_B,
    const Options::Context& context)
    : initial_time_(initial_time),
      expansion_map_options_(expansion_map_options),
      rotation_options_(rotation_options),
      shape_options_A_(shape_options_A),
      shape_options_B_(shape_options_B) {
  if (not(expansion_map_options_.has_value() or rotation_options_.has_value() or
          shape_options_A_.has_value() or shape_options_B_.has_value())) {
    PARSE_ERROR(context,
                "Time dependent map options were specified, but all options "
                "were 'None'. If you don't want time dependent maps, specify "
                "'None' for the TimeDependentMapOptions. If you want time "
                "dependent maps, specify options for at least one map.");
  }

  const auto check_l_max = [&context](const auto& shape_option,
                                      const domain::ObjectLabel label) {
    if (shape_option.has_value() and shape_option.value().l_max <= 1) {
      PARSE_ERROR(context, "Initial LMax for object "
                               << label << " must be 2 or greater but is "
                               << shape_option.value().l_max << " instead.");
    }
  };

  check_l_max(shape_options_A_, domain::ObjectLabel::A);
  check_l_max(shape_options_B_, domain::ObjectLabel::B);
}

template <bool IsCylindrical>
std::unordered_map<std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
TimeDependentMapOptions<IsCylindrical>::create_functions_of_time(
    const std::unordered_map<std::string, double>& initial_expiration_times)
    const {
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      result{};

  // Get existing function of time names that are used for the maps and assign
  // their initial expiration time to infinity (i.e. not expiring)
  std::unordered_map<std::string, double> expiration_times{
      {expansion_name, std::numeric_limits<double>::infinity()},
      {rotation_name, std::numeric_limits<double>::infinity()},
      {gsl::at(size_names, 0), std::numeric_limits<double>::infinity()},
      {gsl::at(size_names, 1), std::numeric_limits<double>::infinity()},
      {gsl::at(shape_names, 0), std::numeric_limits<double>::infinity()},
      {gsl::at(shape_names, 1), std::numeric_limits<double>::infinity()}};

  // If we have control systems, overwrite these expiration times with the ones
  // supplied by the control system
  for (const auto& [name, expr_time] : initial_expiration_times) {
    expiration_times[name] = expr_time;
  }

  // ExpansionMap FunctionOfTime for the function \f$a(t)\f$ in the
  // domain::CoordinateMaps::TimeDependent::CubicScale map
  if (expansion_map_options_.has_value()) {
    result[expansion_name] =
        std::make_unique<FunctionsOfTime::PiecewisePolynomial<2>>(
            initial_time_,
            std::array<DataVector, 3>{
                {{gsl::at(expansion_map_options_.value().initial_values, 0)},
                 {gsl::at(expansion_map_options_.value().initial_values, 1)},
                 {0.0}}},
            expiration_times.at(expansion_name));

    // ExpansionMap FunctionOfTime for the function \f$b(t)\f$ in the
    // domain::CoordinateMaps::TimeDependent::CubicScale map
    result[expansion_outer_boundary_name] =
        std::make_unique<FunctionsOfTime::FixedSpeedCubic>(
            1.0, initial_time_,
            expansion_map_options_.value().outer_boundary_velocity,
            expansion_map_options_.value().outer_boundary_decay_time);
  }

  // RotationMap FunctionOfTime for the rotation angles about each
  // axis.  The initial rotation angles don't matter as we never
  // actually use the angles themselves. We only use their derivatives
  // (omega) to determine map parameters. In theory we could determine
  // each initial angle from the input axis-angle representation, but
  // we don't need to.
  if (rotation_options_.has_value()) {
    result[rotation_name] = std::make_unique<
        FunctionsOfTime::QuaternionFunctionOfTime<3>>(
        initial_time_,
        std::array<DataVector, 1>{DataVector{1.0, 0.0, 0.0, 0.0}},
        std::array<DataVector, 4>{
            {{3, 0.0},
             {gsl::at(rotation_options_.value().initial_angular_velocity, 0),
              gsl::at(rotation_options_.value().initial_angular_velocity, 1),
              gsl::at(rotation_options_.value().initial_angular_velocity, 2)},
             {3, 0.0},
             {3, 0.0}}},
        expiration_times.at(rotation_name));
  }

  // Size and Shape FunctionOfTime for objects A and B
  for (size_t i = 0; i < shape_names.size(); i++) {
    if (i == 0 ? shape_options_A_.has_value() : shape_options_B_.has_value()) {
      const auto make_initial_size_values = [](const auto& lambda_options) {
        return std::array<DataVector, 4>{
            {{gsl::at(lambda_options.value().initial_size_values, 0)},
             {gsl::at(lambda_options.value().initial_size_values, 1)},
             {gsl::at(lambda_options.value().initial_size_values, 2)},
             {0.0}}};
      };

      const size_t initial_l_max = i == 0 ? shape_options_A_.value().l_max
                                          : shape_options_B_.value().l_max;
      const std::array<DataVector, 4> initial_size_values =
          i == 0 ? make_initial_size_values(shape_options_A_)
                 : make_initial_size_values(shape_options_B_);
      const DataVector shape_zeros{
          ylm::Spherepack::spectral_size(initial_l_max, initial_l_max), 0.0};

      result[gsl::at(shape_names, i)] =
          std::make_unique<FunctionsOfTime::PiecewisePolynomial<2>>(
              initial_time_,
              std::array<DataVector, 3>{shape_zeros, shape_zeros, shape_zeros},
              expiration_times.at(gsl::at(shape_names, i)));
      result[gsl::at(size_names, i)] =
          std::make_unique<FunctionsOfTime::PiecewisePolynomial<3>>(
              initial_time_, initial_size_values,
              expiration_times.at(gsl::at(size_names, i)));
    }
  }

  return result;
}

template <bool IsCylindrical>
void TimeDependentMapOptions<IsCylindrical>::build_maps(
    const std::array<std::array<double, 3>, 2>& centers,
    const std::optional<std::array<double, IsCylindrical ? 2 : 3>>&
        object_A_radii,
    const std::optional<std::array<double, IsCylindrical ? 2 : 3>>&
        object_B_radii,
    const double domain_outer_radius) {
  if (expansion_map_options_.has_value()) {
    expansion_map_ = Expansion{domain_outer_radius, expansion_name,
                               expansion_outer_boundary_name};
  }
  if (rotation_options_.has_value()) {
    rotation_map_ = Rotation{rotation_name};
  }

  for (size_t i = 0; i < 2; i++) {
    const auto& radii_opt = i == 0 ? object_A_radii : object_B_radii;
    if (radii_opt.has_value()) {
      const auto& radii = radii_opt.value();
      if (not(i == 0 ? shape_options_A_.has_value()
                     : shape_options_B_.has_value())) {
        ERROR_NO_TRACE(
            "Trying to build the shape map for object "
            << (i == 0 ? domain::ObjectLabel::A : domain::ObjectLabel::B)
            << ", but no time dependent map options were specified "
               "for that object.");
      }

      const size_t initial_l_max = i == 0 ? shape_options_A_.value().l_max
                                          : shape_options_B_.value().l_max;

      std::unique_ptr<domain::CoordinateMaps::ShapeMapTransitionFunctions::
                          ShapeMapTransitionFunction>
          transition_func{};

      // Currently, we don't support different transition functions for the
      // cylindrical domain
      if constexpr (IsCylindrical) {
        transition_func =
            std::make_unique<domain::CoordinateMaps::
                                 ShapeMapTransitionFunctions::SphereTransition>(
                radii[0], radii[1]);

        gsl::at(shape_maps_, i) =
            Shape{gsl::at(centers, i),     initial_l_max,
                  initial_l_max,           std::move(transition_func),
                  gsl::at(shape_names, i), gsl::at(size_names, i)};
      } else {
        // These must match the order of orientations_for_sphere_wrappings() in
        // DomainHelpers.hpp
        const std::array<size_t, 6> axes{2, 2, 1, 1, 0, 0};

        const bool transition_ends_at_cube =
            i == 0 ? shape_options_A_->transition_ends_at_cube
                   : shape_options_B_->transition_ends_at_cube;

        using FalloffEnum =
            domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff;
        const FalloffEnum falloff =
            i == 0 ? shape_options_A_->falloff : shape_options_B_->falloff;

        const double inner_sphericity = 1.0;
        const double outer_sphericity = transition_ends_at_cube ? 0.0 : 1.0;

        const double inner_radius = radii[0];
        const double outer_radius =
            transition_ends_at_cube ? radii[2] : radii[1];

        for (size_t j = 0; j < axes.size(); j++) {
          Parallel::printf(
              "Object %s. Axis %d. Inner radius %.16f. Outer radius %.16f\n",
              i == 0 ? "A" : "B", gsl::at(axes, j), inner_radius, outer_radius);
          transition_func = std::make_unique<
              domain::CoordinateMaps::ShapeMapTransitionFunctions::Wedge>(
              inner_radius, outer_radius, inner_sphericity, outer_sphericity,
              falloff, gsl::at(axes, j));

          gsl::at(gsl::at(shape_maps_, i), j) =
              Shape{gsl::at(centers, i),     initial_l_max,
                    initial_l_max,           std::move(transition_func),
                    gsl::at(shape_names, i), gsl::at(size_names, i)};
        }
      }

    } else if (i == 0 ? shape_options_A_.has_value()
                      : shape_options_B_.has_value()) {
      ERROR_NO_TRACE(
          "No excision was specified for object "
          << (i == 0 ? domain::ObjectLabel::A : domain::ObjectLabel::B)
          << ", but ShapeMap options were specified for that object.");
    }
  }
}

template <bool IsCylindrical>
bool TimeDependentMapOptions<IsCylindrical>::has_distorted_frame_options(
    domain::ObjectLabel object) const {
  ASSERT(object == domain::ObjectLabel::A or object == domain::ObjectLabel::B,
         "object label for TimeDependentMapOptions must be either A or B, not"
             << object);
  return object == domain::ObjectLabel::A ? shape_options_A_.has_value()
                                          : shape_options_B_.has_value();
}

template <bool IsCylindrical>
template <domain::ObjectLabel Object>
typename TimeDependentMapOptions<IsCylindrical>::template MapType<
    Frame::Distorted, Frame::Inertial>
TimeDependentMapOptions<IsCylindrical>::distorted_to_inertial_map(
    const IncludeDistortedMapType& include_distorted_map) const {
  bool block_has_shape_map = false;

  if constexpr (IsCylindrical) {
    block_has_shape_map = include_distorted_map;
  } else {
    const bool transition_ends_at_cube =
        Object == domain::ObjectLabel::A
            ? shape_options_A_->transition_ends_at_cube
            : shape_options_B_->transition_ends_at_cube;
    block_has_shape_map =
        include_distorted_map.has_value() and
        (transition_ends_at_cube or include_distorted_map.value() < 6);
  }

  if (block_has_shape_map) {
    if (expansion_map_.has_value() and rotation_map_.has_value()) {
      return std::make_unique<detail::di_map<Expansion, Rotation>>(
          expansion_map_.value(), rotation_map_.value());
    } else if (expansion_map_.has_value()) {
      return std::make_unique<detail::di_map<Expansion>>(
          expansion_map_.value());
    } else if (rotation_map_.has_value()) {
      return std::make_unique<detail::di_map<Rotation>>(rotation_map_.value());
    } else {
      return std::make_unique<detail::di_map<Identity>>(Identity{});
    }
  } else {
    return nullptr;
  }
}

template <bool IsCylindrical>
template <domain::ObjectLabel Object>
typename TimeDependentMapOptions<IsCylindrical>::template MapType<
    Frame::Grid, Frame::Distorted>
TimeDependentMapOptions<IsCylindrical>::grid_to_distorted_map(
    const IncludeDistortedMapType& include_distorted_map) const {
  bool block_has_shape_map = false;

  if constexpr (IsCylindrical) {
    block_has_shape_map = include_distorted_map;
  } else {
    const bool transition_ends_at_cube =
        Object == domain::ObjectLabel::A
            ? shape_options_A_->transition_ends_at_cube
            : shape_options_B_->transition_ends_at_cube;
    block_has_shape_map =
        include_distorted_map.has_value() and
        (transition_ends_at_cube or include_distorted_map.value() < 6);
  }

  if (block_has_shape_map) {
    const size_t index = get_index(Object);
    const std::optional<Shape>* shape{};
    if constexpr (IsCylindrical) {
      shape = &gsl::at(shape_maps_, index);
    } else {
      if (include_distorted_map.value() >= 12) {
        ERROR(
            "Invalid 'include_distorted_map' argument. Max value allowed is "
            "11, but it is "
            << include_distorted_map.value());
      }
      shape = &gsl::at(gsl::at(shape_maps_, index),
                       include_distorted_map.value() % 6);
    }
    if (not shape->has_value()) {
      ERROR(
          "Requesting grid to distorted map with distorted frame but shape map "
          "options were not specified.");
    }
    return std::make_unique<detail::gd_map<Shape>>(shape->value());
  } else {
    return nullptr;
  }
}

template <bool IsCylindrical>
template <domain::ObjectLabel Object>
typename TimeDependentMapOptions<IsCylindrical>::template MapType<
    Frame::Grid, Frame::Inertial>
TimeDependentMapOptions<IsCylindrical>::grid_to_inertial_map(
    const IncludeDistortedMapType& include_distorted_map) const {
  bool block_has_shape_map = false;

  if constexpr (IsCylindrical) {
    block_has_shape_map = include_distorted_map;
  } else {
    const bool transition_ends_at_cube =
        Object == domain::ObjectLabel::A
            ? shape_options_A_->transition_ends_at_cube
            : shape_options_B_->transition_ends_at_cube;
    block_has_shape_map =
        include_distorted_map.has_value() and
        (transition_ends_at_cube or include_distorted_map.value() < 6);
  }

  if (block_has_shape_map) {
    const size_t index = get_index(Object);
    const std::optional<Shape>* shape{};
    if constexpr (IsCylindrical) {
      shape = &gsl::at(shape_maps_, index);
    } else {
      if (include_distorted_map.value() >= 12) {
        ERROR(
            "Invalid 'include_distorted_map' argument. Max value allowed is "
            "11, but it is "
            << include_distorted_map.value());
      }
      shape = &gsl::at(gsl::at(shape_maps_, index),
                       include_distorted_map.value() % 6);
    }
    if (not shape->has_value()) {
      ERROR(
          "Requesting grid to inertial map with distorted frame but shape map "
          "options were not specified.");
    }
    if (expansion_map_.has_value() and rotation_map_.has_value()) {
      return std::make_unique<detail::gi_map<Shape, Expansion, Rotation>>(
          shape->value(), expansion_map_.value(), rotation_map_.value());
    } else if (expansion_map_.has_value()) {
      return std::make_unique<detail::gi_map<Shape, Expansion>>(
          shape->value(), expansion_map_.value());
    } else if (rotation_map_.has_value()) {
      return std::make_unique<detail::gi_map<Shape, Rotation>>(
          shape->value(), rotation_map_.value());
    } else {
      return std::make_unique<detail::gi_map<Shape>>(shape->value());
    }
  } else {
    if (expansion_map_.has_value() and rotation_map_.has_value()) {
      return std::make_unique<detail::gi_map<Expansion, Rotation>>(
          expansion_map_.value(), rotation_map_.value());
    } else if (expansion_map_.has_value()) {
      return std::make_unique<detail::gi_map<Expansion>>(
          expansion_map_.value());
    } else if (rotation_map_.has_value()) {
      return std::make_unique<detail::gi_map<Rotation>>(rotation_map_.value());
    } else {
      ERROR(
          "Requesting grid to inertial map without a distorted frame and "
          "without a Rotation or Expansion map for object "
          << Object
          << ". This means there are no time dependent maps. If you don't want "
             "time dependent maps, specify 'None' for TimeDependentMapOptions. "
             "Otherwise specify at least one time dependent map.");
    }
  }
}

template <bool IsCylindrical>
size_t TimeDependentMapOptions<IsCylindrical>::get_index(
    const domain::ObjectLabel object) {
  ASSERT(object == domain::ObjectLabel::A or object == domain::ObjectLabel::B,
         "object label for TimeDependentMapOptions must be either A or B, not"
             << object);
  return object == domain::ObjectLabel::A ? 0_st : 1_st;
}

template class TimeDependentMapOptions<true>;
template class TimeDependentMapOptions<false>;

#define ISCYL(data) BOOST_PP_TUPLE_ELEM(0, data)
#define OBJECT(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data)                                                 \
  template TimeDependentMapOptions<ISCYL(data)>::MapType<Frame::Distorted,   \
                                                         Frame::Inertial>    \
  TimeDependentMapOptions<ISCYL(data)>::distorted_to_inertial_map<OBJECT(    \
      data)>(                                                                \
      const TimeDependentMapOptions<ISCYL(data)>::IncludeDistortedMapType&)  \
      const;                                                                 \
  template TimeDependentMapOptions<ISCYL(data)>::MapType<Frame::Grid,        \
                                                         Frame::Distorted>   \
  TimeDependentMapOptions<ISCYL(data)>::grid_to_distorted_map<OBJECT(data)>( \
      const TimeDependentMapOptions<ISCYL(data)>::IncludeDistortedMapType&)  \
      const;                                                                 \
  template TimeDependentMapOptions<ISCYL(data)>::MapType<Frame::Grid,        \
                                                         Frame::Inertial>    \
  TimeDependentMapOptions<ISCYL(data)>::grid_to_inertial_map<OBJECT(data)>(  \
      const TimeDependentMapOptions<ISCYL(data)>::IncludeDistortedMapType&)  \
      const;

GENERATE_INSTANTIATIONS(INSTANTIATE, (true, false),
                        (domain::ObjectLabel::A, domain::ObjectLabel::B,
                         domain::ObjectLabel::None))

#undef OBJECT
#undef ISCYL
#undef INSTANTIATE
}  // namespace domain::creators::bco
