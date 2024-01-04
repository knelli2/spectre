// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticSolutions/Xcts/CommonVariables.hpp"

#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/AnalyticData/Xcts/CommonVariables.tpp"
#include "PointwiseFunctions/Elasticity/Strain.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Xcts/LongitudinalOperator.hpp"
#include "Utilities/Gsl.hpp"

namespace Xcts::Solutions {

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<Scalar<DataType>*> conformal_factor,
    const gsl::not_null<Cache*> cache,
    Tags::ConformalFactor<DataType> /*meta*/) const {
  const auto& conformal_factor_minus_one =
      cache->get_var(*this, Tags::ConformalFactorMinusOne<DataType>{});
  get(*conformal_factor) = get(conformal_factor_minus_one) + 1.;
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<Scalar<DataType>*> lapse_times_conformal_factor,
    const gsl::not_null<Cache*> cache,
    Tags::LapseTimesConformalFactor<DataType> /*meta*/) const {
  const auto& lapse_times_conformal_factor_minus_one = cache->get_var(
      *this, Tags::LapseTimesConformalFactorMinusOne<DataType>{});
  get(*lapse_times_conformal_factor) =
      get(lapse_times_conformal_factor_minus_one) + 1.;
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<tnsr::I<DataType, Dim>*> conformal_factor_flux,
    const gsl::not_null<Cache*> cache,
    ::Tags::Flux<Tags::ConformalFactorMinusOne<DataType>, tmpl::size_t<Dim>,
                 Frame::Inertial> /*meta*/) const {
  const auto& conformal_factor_gradient = cache->get_var(
      *this, ::Tags::deriv<Tags::ConformalFactorMinusOne<DataType>,
                           tmpl::size_t<Dim>, Frame::Inertial>{});
  const auto& inv_conformal_metric = cache->get_var(
      *this, Tags::InverseConformalMetric<DataType, Dim, Frame::Inertial>{});
  raise_or_lower_index(conformal_factor_flux, conformal_factor_gradient,
                       inv_conformal_metric);
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<tnsr::ii<DataType, Dim>*> shift_strain,
    const gsl::not_null<Cache*> cache,
    Tags::ShiftStrain<DataType, Dim, Frame::Inertial> /*meta*/) const {
  const auto& shift_excess = cache->get_var(
      *this, Tags::ShiftExcess<DataType, Dim, Frame::Inertial>{});
  const auto& deriv_shift_excess = cache->get_var(
      *this, ::Tags::deriv<Tags::ShiftExcess<DataType, Dim, Frame::Inertial>,
                           tmpl::size_t<Dim>, Frame::Inertial>{});
  const auto& conformal_metric = cache->get_var(
      *this, Xcts::Tags::ConformalMetric<DataType, Dim, Frame::Inertial>{});
  const auto& deriv_conformal_metric = cache->get_var(
      *this,
      ::Tags::deriv<Xcts::Tags::ConformalMetric<DataType, Dim, Frame::Inertial>,
                    tmpl::size_t<Dim>, Frame::Inertial>{});
  const auto& conformal_christoffel_first_kind = cache->get_var(
      *this, Xcts::Tags::ConformalChristoffelFirstKind<DataType, Dim,
                                                       Frame::Inertial>{});
  Elasticity::strain(shift_strain, deriv_shift_excess, conformal_metric,
                     deriv_conformal_metric, conformal_christoffel_first_kind,
                     shift_excess);
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<tnsr::I<DataType, Dim>*>
        lapse_times_conformal_factor_flux,
    const gsl::not_null<Cache*> cache,
    ::Tags::Flux<Tags::LapseTimesConformalFactorMinusOne<DataType>,
                 tmpl::size_t<Dim>, Frame::Inertial> /*meta*/) const {
  const auto& lapse_times_conformal_factor_gradient = cache->get_var(
      *this, ::Tags::deriv<Tags::LapseTimesConformalFactorMinusOne<DataType>,
                           tmpl::size_t<Dim>, Frame::Inertial>{});
  const auto& inv_conformal_metric = cache->get_var(
      *this, Tags::InverseConformalMetric<DataType, Dim, Frame::Inertial>{});
  raise_or_lower_index(lapse_times_conformal_factor_flux,
                       lapse_times_conformal_factor_gradient,
                       inv_conformal_metric);
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<tnsr::II<DataType, Dim>*> longitudinal_shift_excess,
    const gsl::not_null<Cache*> cache,
    Tags::LongitudinalShiftExcess<DataType, Dim, Frame::Inertial> /*meta*/)
    const {
  const auto& shift_excess = cache->get_var(
      *this, Tags::ShiftExcess<DataType, Dim, Frame::Inertial>{});
  const auto& deriv_shift_excess = cache->get_var(
      *this, ::Tags::deriv<Tags::ShiftExcess<DataType, Dim, Frame::Inertial>,
                           tmpl::size_t<Dim>, Frame::Inertial>{});
  const auto& inv_conformal_metric = cache->get_var(
      *this, Tags::InverseConformalMetric<DataType, Dim, Frame::Inertial>{});
  const auto& conformal_christoffel_second_kind = cache->get_var(
      *this,
      Tags::ConformalChristoffelSecondKind<DataType, Dim, Frame::Inertial>{});
  Xcts::longitudinal_operator(longitudinal_shift_excess, shift_excess,
                              deriv_shift_excess, inv_conformal_metric,
                              conformal_christoffel_second_kind);
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<tnsr::I<DataType, Dim>*> shift,
    const gsl::not_null<Cache*> cache,
    gr::Tags::Shift<DataType, Dim> /*meta*/) const {
  *shift = cache->get_var(*this,
                          Tags::ShiftExcess<DataType, Dim, Frame::Inertial>{});
  const auto& shift_background = cache->get_var(
      *this, Tags::ShiftBackground<DataType, Dim, Frame::Inertial>{});
  for (size_t d = 0; d < Dim; ++d) {
    shift->get(d) += shift_background.get(d);
  }
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<tnsr::ii<DataType, Dim>*> spatial_metric,
    const gsl::not_null<Cache*> cache,
    gr::Tags::SpatialMetric<DataType, Dim> /*meta*/) const {
  *spatial_metric = cache->get_var(
      *this, Tags::ConformalMetric<DataType, Dim, Frame::Inertial>{});
  const auto& conformal_factor =
      cache->get_var(*this, Tags::ConformalFactor<DataType>{});
  for (size_t i = 0; i < spatial_metric->size(); ++i) {
    (*spatial_metric)[i] *= pow<4>(get(conformal_factor));
  }
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<tnsr::II<DataType, Dim>*> inv_spatial_metric,
    const gsl::not_null<Cache*> cache,
    gr::Tags::InverseSpatialMetric<DataType, Dim> /*meta*/)
    const {
  const auto& conformal_factor =
      cache->get_var(*this, Tags::ConformalFactor<DataType>{});
  *inv_spatial_metric = cache->get_var(
      *this, Tags::InverseConformalMetric<DataType, Dim, Frame::Inertial>{});
  for (size_t i = 0; i < inv_spatial_metric->size(); ++i) {
    (*inv_spatial_metric)[i] /= pow<4>(get(conformal_factor));
  }
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<tnsr::ijj<DataType, Dim>*> deriv_spatial_metric,
    const gsl::not_null<Cache*> cache,
    ::Tags::deriv<gr::Tags::SpatialMetric<DataType, Dim>,
                  tmpl::size_t<Dim>, Frame::Inertial> /*meta*/) const {
  const auto& conformal_metric = cache->get_var(
      *this, Tags::ConformalMetric<DataType, Dim, Frame::Inertial>{});
  const auto& conformal_factor =
      cache->get_var(*this, Tags::ConformalFactor<DataType>{});
  const auto& deriv_conformal_factor = cache->get_var(
      *this, ::Tags::deriv<Tags::ConformalFactorMinusOne<DataType>,
                           tmpl::size_t<Dim>, Frame::Inertial>{});
  *deriv_spatial_metric = cache->get_var(
      *this,
      ::Tags::deriv<Tags::ConformalMetric<DataType, Dim, Frame::Inertial>,
                    tmpl::size_t<3>, Frame::Inertial>{});
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k <= j; ++k) {
        deriv_spatial_metric->get(i, j, k) *= pow<4>(get(conformal_factor));
        deriv_spatial_metric->get(i, j, k) +=
            4. * pow<3>(get(conformal_factor)) * deriv_conformal_factor.get(i) *
            conformal_metric.get(j, k);
      }
    }
  }
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<Scalar<DataType>*>
        longitudinal_shift_minus_dt_conformal_metric_square,
    const gsl::not_null<Cache*> cache,
    Tags::LongitudinalShiftMinusDtConformalMetricSquare<DataType> /*meta*/)
    const {
  const auto& longitudinal_shift_background = cache->get_var(
      *this, Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
                 DataType, Dim, Frame::Inertial>{});
  const auto& longitudinal_shift_excess = cache->get_var(
      *this, Tags::LongitudinalShiftExcess<DataType, Dim, Frame::Inertial>{});
  const auto& conformal_metric = cache->get_var(
      *this, Tags::ConformalMetric<DataType, Dim, Frame::Inertial>{});
  get(*longitudinal_shift_minus_dt_conformal_metric_square) = 0.;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 3; ++k) {
        for (size_t l = 0; l < 3; ++l) {
          get(*longitudinal_shift_minus_dt_conformal_metric_square) +=
              conformal_metric.get(i, k) * conformal_metric.get(j, l) *
              (longitudinal_shift_background.get(i, j) +
               longitudinal_shift_excess.get(i, j)) *
              (longitudinal_shift_background.get(k, l) +
               longitudinal_shift_excess.get(k, l));
        }
      }
    }
  }
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<Scalar<DataType>*>
        longitudinal_shift_minus_dt_conformal_metric_over_lapse_square,
    const gsl::not_null<Cache*> cache,
    Tags::LongitudinalShiftMinusDtConformalMetricOverLapseSquare<
        DataType> /*meta*/) const {
  *longitudinal_shift_minus_dt_conformal_metric_over_lapse_square =
      cache->get_var(
          *this,
          Tags::LongitudinalShiftMinusDtConformalMetricSquare<DataType>{});
  const auto& lapse = cache->get_var(*this, gr::Tags::Lapse<DataType>{});
  get(*longitudinal_shift_minus_dt_conformal_metric_over_lapse_square) /=
      square(get(lapse));
}

template <typename DataType, typename Cache>
void CommonVariables<DataType, Cache>::operator()(
    const gsl::not_null<Scalar<DataType>*>
        shift_dot_deriv_extrinsic_curvature_trace,
    const gsl::not_null<Cache*> cache,
    Tags::ShiftDotDerivExtrinsicCurvatureTrace<DataType> /*meta*/) const {
  const auto& shift =
      cache->get_var(*this, gr::Tags::Shift<DataType, Dim>{});
  const auto& deriv_extrinsic_curvature_trace = cache->get_var(
      *this, ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataType>,
                           tmpl::size_t<Dim>, Frame::Inertial>{});
  dot_product(shift_dot_deriv_extrinsic_curvature_trace, shift,
              deriv_extrinsic_curvature_trace);
}

}  // namespace Xcts::Solutions
