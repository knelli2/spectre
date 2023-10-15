// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/DgSubcell/Tags/Reconstructor.hpp"
#include "Evolution/DgSubcell/Tags/SubcellSolver.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/FiniteDifference/FilterOptions.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/FiniteDifference/Reconstructor.hpp"
#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

namespace grmhd::GhValenciaDivClean::fd {
/// Option tags for finite difference solver
namespace OptionTags {
/// \brief Option tag for the reconstructor
struct Reconstructor {
  using type = std::unique_ptr<fd::Reconstructor>;

  inline const static std::string help {"The reconstruction scheme to use."};
  using group = evolution::dg::subcell::OptionTags::SubcellSolverGroup;
};

/// \brief Option tag for the filter/dissipation options.
struct FilterOptions {
  using type = fd::FilterOptions;

  inline const static std::string help =type::help;
  using group = evolution::dg::subcell::OptionTags::SubcellSolverGroup;
};
}  // namespace OptionTags

/// %Tags for finite difference solver
namespace Tags {
/// \brief Tag for the reconstructor
struct Reconstructor : db::SimpleTag,
                       evolution::dg::subcell::Tags::Reconstructor {
  using type = std::unique_ptr<fd::Reconstructor>;
  using option_tags = tmpl::list<OptionTags::Reconstructor>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& reconstructor) {
    return reconstructor->get_clone();
  }
};

/// \brief Tag for filter/dissipation options.
struct FilterOptions : db::SimpleTag {
  using type = fd::FilterOptions;
  using option_tags = tmpl::list<OptionTags::FilterOptions>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& filter_options) {
    return filter_options;
  }
};
}  // namespace Tags
}  // namespace grmhd::GhValenciaDivClean::fd
