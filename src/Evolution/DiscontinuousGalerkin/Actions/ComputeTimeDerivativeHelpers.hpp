// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <type_traits>

#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/CreateHasTypeAlias.hpp"

namespace evolution::dg::Actions::detail {
CREATE_HAS_TYPE_ALIAS(boundary_correction_base)
CREATE_HAS_TYPE_ALIAS_V(boundary_correction_base)

CREATE_HAS_TYPE_ALIAS(boundary_conditions_base)
CREATE_HAS_TYPE_ALIAS_V(boundary_conditions_base)

CREATE_HAS_TYPE_ALIAS(inverse_spatial_metric_tag)
CREATE_HAS_TYPE_ALIAS_V(inverse_spatial_metric_tag)

template <bool HasInverseSpatialMetricTag = false>
struct inverse_spatial_metric_tag_impl {
  template <typename System>
  using f = tmpl::list<>;
};

template <>
struct inverse_spatial_metric_tag_impl<true> {
  template <typename System>
  using f = tmpl::list<typename System::inverse_spatial_metric_tag>;
};

template <typename System>
using inverse_spatial_metric_tag = typename inverse_spatial_metric_tag_impl<
    has_inverse_spatial_metric_tag_v<System>>::template f<System>;

template <bool HasPrimitiveVars = false>
struct get_primitive_vars {
  template <typename BoundaryCorrection>
  using f = tmpl::list<>;

  template <typename BoundaryCondition>
  using boundary_condition_interior_tags = tmpl::list<>;
};

template <>
struct get_primitive_vars<true> {
  template <typename BoundaryCorrection>
  using f = typename BoundaryCorrection::dg_package_data_primitive_tags;

  template <typename BoundaryCondition>
  using boundary_condition_interior_tags =
      typename BoundaryCondition::dg_interior_primitive_variables_tags;
};

template <bool HasPrimitiveAndConservativeVars, typename BoundaryCorrection>
using boundary_correction_primitive_tags = typename get_primitive_vars<
    HasPrimitiveAndConservativeVars>::template f<BoundaryCorrection>;

template <bool HasPrimitiveAndConservativeVars, typename BoundaryCondition>
using boundary_condition_primitive_tags =
    typename get_primitive_vars<HasPrimitiveAndConservativeVars>::
        template boundary_condition_interior_tags<BoundaryCondition>;

template <typename BoundaryCorrection, typename = std::void_t<>>
struct interior_tags_for_boundary_correction {
  using type = tmpl::list<>;
};

template <typename BoundaryCorrection>
struct interior_tags_for_boundary_correction<
    BoundaryCorrection,
    std::void_t<typename BoundaryCorrection::
                    dg_project_from_interior_for_boundary_condition>> {
  using type = typename BoundaryCorrection::
      dg_project_from_interior_for_boundary_condition;
};

template <typename BoundaryCondition, typename = std::void_t<>>
struct derivative_tags_for_boundary_condition {
  using type = tmpl::list<>;
};

template <typename BoundaryCondition>
struct derivative_tags_for_boundary_condition<
    BoundaryCondition,
    std::void_t<typename BoundaryCondition::dg_interior_derivative_tags>> {
  using type = typename BoundaryCondition::dg_interior_derivative_tags;
};

template <typename System, bool = System::has_primitive_and_conservative_vars>
struct get_primitive_vars_tags_from_system_impl {
  using type = typename System::primitive_variables_tag::tags_list;
};

template <typename System>
struct get_primitive_vars_tags_from_system_impl<System, false> {
  using type = tmpl::list<>;
};

/// Returns a `tmpl::list` of the primitive tags. The list is empty if the
/// system does not have primitive tags.
template <typename System>
using get_primitive_vars_tags_from_system =
    typename get_primitive_vars_tags_from_system_impl<System>::type;

template <typename BoundaryCondition, typename = std::void_t<>>
struct get_dt_vars_from_boundary_condition_impl {
  using type = tmpl::list<>;
};

template <typename BoundaryCondition>
struct get_dt_vars_from_boundary_condition_impl<
    BoundaryCondition,
    std::void_t<typename BoundaryCondition::dg_interior_dt_vars_tags>> {
  using type = typename BoundaryCondition::dg_interior_dt_vars_tags;
};

/// Returns the `dg_interior_dt_vars_tags` if the boundary condition specifies
/// them, otherwise returns an empty list.
template <typename BoundaryCondition>
using get_dt_vars_from_boundary_condition =
    typename get_dt_vars_from_boundary_condition_impl<BoundaryCondition>::type;

template <typename BoundaryCondition, typename = std::void_t<>>
struct get_deriv_vars_from_boundary_condition_impl {
  using type = tmpl::list<>;
};

template <typename BoundaryCondition>
struct get_deriv_vars_from_boundary_condition_impl<
    BoundaryCondition,
    std::void_t<typename BoundaryCondition::dg_interior_deriv_vars_tags>> {
  using type = typename BoundaryCondition::dg_interior_deriv_vars_tags;
};

/// Returns the `dg_interior_deriv_vars_tags` if the boundary condition
/// specifies them, otherwise returns an empty list.
template <typename BoundaryCondition>
using get_deriv_vars_from_boundary_condition =
    typename get_deriv_vars_from_boundary_condition_impl<
        BoundaryCondition>::type;

template <typename BoundaryCorrection, typename = std::void_t<>>
struct get_sarah_list_impl {
  using type = tmpl::list<>;
};

template <typename BoundaryCorrection>
struct get_sarah_list_impl<
    BoundaryCorrection, std::void_t<typename BoundaryCorrection::sarah_list>> {
  using type = typename BoundaryCorrection::sarah_list;
};

template <typename BoundaryCorrection>
using get_sarah_list = typename get_sarah_list_impl<BoundaryCorrection>::type;
}  // namespace evolution::dg::Actions::detail
