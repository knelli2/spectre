// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/DgSubcell/Tags/Reconstructor.hpp"
#include "Evolution/DgSubcell/Tags/SubcellSolver.hpp"
#include "Evolution/Systems/NewtonianEuler/FiniteDifference/Reconstructor.hpp"
#include "Options/String.hpp"

namespace NewtonianEuler::fd {
/// Option tags for reconstruction
namespace OptionTags {
/// \brief Option tag for the reconstructor
template <size_t Dim>
struct Reconstructor {
  using type = std::unique_ptr<fd::Reconstructor<Dim>>;

  inline const static std::string help {"The reconstruction scheme to use."};
  using group = evolution::dg::subcell::OptionTags::SubcellSolverGroup;
};
}  // namespace OptionTags

/// %Tags for reconstruction
namespace Tags {
/// \brief Tag for the reconstructor
template <size_t Dim>
struct Reconstructor : db::SimpleTag,
                       evolution::dg::subcell::Tags::Reconstructor {
  using type = std::unique_ptr<fd::Reconstructor<Dim>>;
  using option_tags = tmpl::list<OptionTags::Reconstructor<Dim>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& reconstructor) {
    return reconstructor->get_clone();
  }
};
}  // namespace Tags
}  // namespace NewtonianEuler::fd
