// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/DgSubcell/Tags/Reconstructor.hpp"
#include "Evolution/DgSubcell/Tags/SubcellSolver.hpp"
#include "Evolution/Systems/Burgers/FiniteDifference/Reconstructor.hpp"
#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

namespace Burgers::fd {
namespace OptionTags {
/*!
 * \brief Holds the subcell reconstructor in the input file
 */
struct Reconstructor {
  using type = std::unique_ptr<fd::Reconstructor>;

  inline const static std::string help {"The reconstruction scheme to use."};
  using group = evolution::dg::subcell::OptionTags::SubcellSolverGroup;
};
}  // namespace OptionTags

namespace Tags {
/*!
 * \brief Tag for the reconstructor
 */
struct Reconstructor : db::SimpleTag,
                       evolution::dg::subcell::Tags::Reconstructor {
  using type = std::unique_ptr<fd::Reconstructor>;
  using option_tags = tmpl::list<OptionTags::Reconstructor>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& reconstructor) {
    return reconstructor->get_clone();
  }
};
}  // namespace Tags
}  // namespace Burgers::fd
