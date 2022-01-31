// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>

#include "DataStructures/DataBox/Tag.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags/OptionsGroup.hpp"
#include "Options/Options.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace ScalarAdvection::subcell {
struct TciOptions {
  /*!
   * \brief The cutoff where the absolute value of the scalar field \f$U\f$ in
   * an element is below this we do not apply the Persson TCI.
   */
  struct UCutoff {
    using type = double;
    static type lower_bound() { return 0.0; }
    static constexpr Options::String help = {
        "The cutoff where the absolute value of the scalar field U in an "
        "element is below this we do not use subcell."};
  };

  using options = tmpl::list<UCutoff>;
  static constexpr Options::String help = {
      "Options for the troubled-cell indicator"};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/);

  double u_cutoff{std::numeric_limits<double>::signaling_NaN()};
};

namespace OptionTags {
struct TciOptions {
  using type = subcell::TciOptions;
  static constexpr Options::String help =
      "ScalarAdvection-specific options for the TCI";
  using group = ::dg::OptionTags::DiscontinuousGalerkinGroup;
};
}  // namespace OptionTags

namespace Tags {
struct TciOptions : db::SimpleTag {
  using type = subcell::TciOptions;
  using option_tags = tmpl::list<typename OptionTags::TciOptions>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& tci_options) {
    return tci_options;
  }
};
}  // namespace Tags
}  // namespace ScalarAdvection::subcell
