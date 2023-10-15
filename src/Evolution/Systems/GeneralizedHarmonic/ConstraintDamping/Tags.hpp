// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/DampingFunction.hpp"
#include "Options/String.hpp"

/// \cond
namespace gh::OptionTags {
struct Group;
}  // namespace gh::OptionTags
/// \endcond

namespace gh::ConstraintDamping {
namespace OptionTags {
template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma0 {
  using type =
      std::unique_ptr<::gh::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  inline const static std::string help{
      "DampingFunction for damping parameter gamma0"};
  using group = gh::OptionTags::Group;
};

template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma1 {
  using type =
      std::unique_ptr<::gh::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  inline const static std::string help{
      "DampingFunction for damping parameter gamma1"};
  using group = gh::OptionTags::Group;
};

template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma2 {
  using type =
      std::unique_ptr<::gh::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  inline const static std::string help{
      "DampingFunction for damping parameter gamma2"};
  using group = gh::OptionTags::Group;
};
}  // namespace OptionTags

namespace Tags {
/*!
 * \brief Constraint dammping parameter \f$\gamma_0\f$ for the generalized
 * harmonic system (cf. \cite Lindblom2005qh).
 */
struct ConstraintGamma0 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

/*!
 * \brief Constraint dammping parameter \f$\gamma_1\f$ for the generalized
 * harmonic system (cf. \cite Lindblom2005qh).
 */
struct ConstraintGamma1 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

/*!
 * \brief Constraint dammping parameter \f$\gamma_2\f$ for the generalized
 * harmonic system (cf. \cite Lindblom2005qh).
 */
struct ConstraintGamma2 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

/*!
 * \brief A DampingFunction to compute the constraint damping parameter
 * \f$\gamma_0\f$.
 */
template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma0 : db::SimpleTag {
  using DampingFunctionType =
      ::gh::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::gh::ConstraintDamping::OptionTags::DampingFunctionGamma0<
          VolumeDim, Fr>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& damping_function) {
    return damping_function->get_clone();
  }
};

/*!
 * \brief A DampingFunction to compute the constraint damping parameter
 * \f$\gamma_0\f$.
 */
template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma1 : db::SimpleTag {
  using DampingFunctionType =
      ::gh::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::gh::ConstraintDamping::OptionTags::DampingFunctionGamma1<
          VolumeDim, Fr>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& damping_function) {
    return damping_function->get_clone();
  }
};

/*!
 * \brief A DampingFunction to compute the constraint damping parameter
 * \f$\gamma_0\f$.
 */
template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma2 : db::SimpleTag {
  using DampingFunctionType =
      ::gh::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::gh::ConstraintDamping::OptionTags::DampingFunctionGamma2<
          VolumeDim, Fr>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& damping_function) {
    return damping_function->get_clone();
  }
};
}  // namespace Tags
}  // namespace gh::ConstraintDamping
