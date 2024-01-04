// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "PointwiseFunctions/AnalyticData/Xcts/CommonVariables.hpp"
#include "PointwiseFunctions/GeneralRelativity/TagsDeclarations.hpp"
#include "Utilities/Gsl.hpp"

namespace Xcts::Solutions {

/// Tags for variables that solutions can share
template <typename DataType>
using common_tags = tmpl::push_back<
    AnalyticData::common_tags<DataType>,
    // Solved variables
    Tags::ConformalFactorMinusOne<DataType>, Tags::ConformalFactor<DataType>,
    Tags::LapseTimesConformalFactorMinusOne<DataType>,
    Tags::LapseTimesConformalFactor<DataType>,
    Tags::ShiftExcess<DataType, 3, Frame::Inertial>,
    // ADM variables
    gr::Tags::Lapse<DataType>, gr::Tags::Shift<DataType, 3>,
    gr::Tags::SpatialMetric<DataType, 3>,
    gr::Tags::InverseSpatialMetric<DataType, 3>,
    ::Tags::deriv<gr::Tags::SpatialMetric<DataType, 3>, tmpl::size_t<3>,
                  Frame::Inertial>,
    gr::Tags::ExtrinsicCurvature<DataType, 3>,
    // Derivatives of solved variables
    ::Tags::deriv<Tags::ConformalFactorMinusOne<DataType>, tmpl::size_t<3>,
                  Frame::Inertial>,
    ::Tags::deriv<Tags::LapseTimesConformalFactorMinusOne<DataType>,
                  tmpl::size_t<3>, Frame::Inertial>,
    ::Tags::deriv<Tags::ShiftExcess<DataType, 3, Frame::Inertial>,
                  tmpl::size_t<3>, Frame::Inertial>,
    Tags::ShiftStrain<DataType, 3, Frame::Inertial>,
    // Fluxes
    ::Tags::Flux<Tags::ConformalFactorMinusOne<DataType>, tmpl::size_t<3>,
                 Frame::Inertial>,
    ::Tags::Flux<Tags::LapseTimesConformalFactorMinusOne<DataType>,
                 tmpl::size_t<3>, Frame::Inertial>,
    Tags::LongitudinalShiftExcess<DataType, 3, Frame::Inertial>,
    // Background quantities for subsets of the XCTS equations
    Tags::LongitudinalShiftMinusDtConformalMetricSquare<DataType>,
    Tags::LongitudinalShiftMinusDtConformalMetricOverLapseSquare<DataType>,
    Tags::ShiftDotDerivExtrinsicCurvatureTrace<DataType>>;

/// Tags for hydro variables that are typically retrieved from a hydro solution
template <typename DataType>
using hydro_tags = AnalyticData::hydro_tags<DataType>;

/// Implementations for variables that solutions can share
template <typename DataType, typename Cache>
struct CommonVariables : AnalyticData::CommonVariables<DataType, Cache> {
  static constexpr size_t Dim = 3;
  using Base = AnalyticData::CommonVariables<DataType, Cache>;
  using Base::Base;
  using Base::operator();

  virtual void operator()(
      gsl::not_null<Scalar<DataType>*> conformal_factor_minus_one,
      gsl::not_null<Cache*> cache,
      Tags::ConformalFactorMinusOne<DataType> /*meta*/) const = 0;
  virtual void operator()(gsl::not_null<Scalar<DataType>*> conformal_factor,
                          gsl::not_null<Cache*> cache,
                          Tags::ConformalFactor<DataType> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<Scalar<DataType>*> lapse_times_conformal_factor_minus_one,
      gsl::not_null<Cache*> cache,
      Tags::LapseTimesConformalFactorMinusOne<DataType> /*meta*/) const = 0;
  virtual void operator()(
      gsl::not_null<Scalar<DataType>*> lapse_times_conformal_factor,
      gsl::not_null<Cache*> cache,
      Tags::LapseTimesConformalFactor<DataType> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<tnsr::I<DataType, Dim>*> shift_excess,
      gsl::not_null<Cache*> cache,
      Tags::ShiftExcess<DataType, Dim, Frame::Inertial> /*meta*/) const = 0;
  virtual void operator()(gsl::not_null<Scalar<DataType>*> lapse,
                          gsl::not_null<Cache*> cache,
                          gr::Tags::Lapse<DataType> /*meta*/) const = 0;
  virtual void operator()(
      gsl::not_null<tnsr::ii<DataType, Dim>*> spatial_metric,
      gsl::not_null<Cache*> cache,
      gr::Tags::SpatialMetric<DataType, Dim> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<tnsr::II<DataType, Dim>*> inv_spatial_metric,
      gsl::not_null<Cache*> cache,
      gr::Tags::InverseSpatialMetric<DataType, Dim> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<tnsr::ijj<DataType, Dim>*> deriv_spatial_metric,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<gr::Tags::SpatialMetric<DataType, Dim>, tmpl::size_t<Dim>,
                    Frame::Inertial> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<tnsr::i<DataType, Dim>*> deriv_conformal_factor,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Tags::ConformalFactorMinusOne<DataType>, tmpl::size_t<Dim>,
                    Frame::Inertial> /*meta*/) const = 0;
  virtual void operator()(
      gsl::not_null<tnsr::i<DataType, Dim>*> deriv_lapse_times_conformal_factor,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Tags::LapseTimesConformalFactorMinusOne<DataType>,
                    tmpl::size_t<Dim>, Frame::Inertial> /*meta*/) const = 0;
  virtual void operator()(
      gsl::not_null<tnsr::iJ<DataType, Dim>*> deriv_shift_excess,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Tags::ShiftExcess<DataType, 3, Frame::Inertial>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const = 0;
  virtual void operator()(
      gsl::not_null<tnsr::ii<DataType, Dim>*> shift_strain,
      gsl::not_null<Cache*> cache,
      Tags::ShiftStrain<DataType, Dim, Frame::Inertial> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<tnsr::I<DataType, Dim>*> conformal_factor_flux,
      gsl::not_null<Cache*> cache,
      ::Tags::Flux<Tags::ConformalFactorMinusOne<DataType>, tmpl::size_t<Dim>,
                   Frame::Inertial> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<tnsr::I<DataType, Dim>*> lapse_times_conformal_factor_flux,
      gsl::not_null<Cache*> cache,
      ::Tags::Flux<Tags::LapseTimesConformalFactorMinusOne<DataType>,
                   tmpl::size_t<Dim>, Frame::Inertial> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<tnsr::II<DataType, Dim>*> longitudinal_shift_excess,
      gsl::not_null<Cache*> cache,
      Tags::LongitudinalShiftExcess<DataType, Dim, Frame::Inertial> /*meta*/)
      const;
  virtual void operator()(gsl::not_null<tnsr::I<DataType, Dim>*> shift,
                          gsl::not_null<Cache*> cache,
                          gr::Tags::Shift<DataType, Dim> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<tnsr::ii<DataType, Dim>*> extrinsic_curvature,
      gsl::not_null<Cache*> cache,
      gr::Tags::ExtrinsicCurvature<DataType, Dim> /*meta*/) const = 0;
  virtual void operator()(
      gsl::not_null<Scalar<DataType>*>
          longitudinal_shift_minus_dt_conformal_metric_square,
      gsl::not_null<Cache*> cache,
      Tags::LongitudinalShiftMinusDtConformalMetricSquare<DataType> /*meta*/)
      const;
  virtual void operator()(
      gsl::not_null<Scalar<DataType>*>
          longitudinal_shift_minus_dt_conformal_metric_over_lapse_square,
      gsl::not_null<Cache*> cache,
      Tags::LongitudinalShiftMinusDtConformalMetricOverLapseSquare<
          DataType> /*meta*/) const;
  virtual void operator()(
      gsl::not_null<Scalar<DataType>*>
          shift_dot_deriv_extrinsic_curvature_trace,
      gsl::not_null<Cache*> cache,
      Tags::ShiftDotDerivExtrinsicCurvatureTrace<DataType> /*meta*/) const;
};

}  // namespace Xcts::Solutions
