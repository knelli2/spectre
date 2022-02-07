// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <type_traits>

#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/Algorithms/AlgorithmGroup.hpp"
#include "Parallel/Algorithms/AlgorithmNodegroup.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/TypeTraits.hpp"

namespace PUP {
class er;
}  // namespace PUP

namespace {
class PupableClass {
 public:
  void pup(PUP::er&) {}  // NOLINT
};

class NonpupableClass {};

struct Metavariables {
  enum class Phase { Initialization, Exit };
};

struct SingletonParallelComponent {
  using chare_type = Parallel::Algorithms::Singleton;
  using metavariables = Metavariables;
  using initialization_tags = tmpl::list<>;
};
struct ArrayParallelComponent {
  using chare_type = Parallel::Algorithms::Array;
  using metavariables = Metavariables;
  using initialization_tags = tmpl::list<>;
};
struct GroupParallelComponent {
  using chare_type = Parallel::Algorithms::Group;
  using metavariables = Metavariables;
  using initialization_tags = tmpl::list<>;
};
struct NodegroupParallelComponent {
  using chare_type = Parallel::Algorithms::Nodegroup;
  using metavariables = Metavariables;
  using initialization_tags = tmpl::list<>;
};

using singleton = SingletonParallelComponent;
using array = ArrayParallelComponent;
using group = GroupParallelComponent;
using nodegroup = NodegroupParallelComponent;

// If passing proxy to a full array, we expect it to be from a
// Parallel::Algorithms::Singleton. If passing proxy to an element, we expect it
// to be from a Parallel::Algorithms::Array.
using array_proxy = CProxy_AlgorithmArray<SingletonParallelComponent, int>;
using array_element_proxy =
    CProxyElement_AlgorithmArray<ArrayParallelComponent, int>;
using group_proxy = CProxy_AlgorithmGroup<GroupParallelComponent, int>;
using nodegroup_proxy =
    CProxy_AlgorithmNodegroup<NodegroupParallelComponent, int>;
}  // namespace

static_assert(Parallel::is_array_proxy<array_proxy>::value,
              "Failed testing type trait is_array_proxy");
static_assert(not Parallel::is_array_proxy<array_element_proxy>::value,
              "Failed testing type trait is_array_proxy");
static_assert(not Parallel::is_array_proxy<group_proxy>::value,
              "Failed testing type trait is_array_proxy");
static_assert(not Parallel::is_array_proxy<nodegroup_proxy>::value,
              "Failed testing type trait is_array_proxy");

static_assert(not Parallel::is_array_element_proxy<array_proxy>::value,
              "Failed testing type trait is_array_element_proxy");
static_assert(Parallel::is_array_element_proxy<array_element_proxy>::value,
              "Failed testing type trait is_array_element_proxy");
static_assert(not Parallel::is_array_element_proxy<group_proxy>::value,
              "Failed testing type trait is_array_element_proxy");
static_assert(not Parallel::is_array_element_proxy<nodegroup_proxy>::value,
              "Failed testing type trait is_array_element_proxy");

static_assert(not Parallel::is_group_proxy<array_proxy>::value,
              "Failed testing type trait is_group_proxy");
static_assert(not Parallel::is_group_proxy<array_element_proxy>::value,
              "Failed testing type trait is_group_proxy");
static_assert(Parallel::is_group_proxy<group_proxy>::value,
              "Failed testing type trait is_group_proxy");
static_assert(not Parallel::is_group_proxy<nodegroup_proxy>::value,
              "Failed testing type trait is_group_proxy");

static_assert(not Parallel::is_node_group_proxy<array_proxy>::value,
              "Failed testing type trait is_node_group_proxy");
static_assert(not Parallel::is_node_group_proxy<array_element_proxy>::value,
              "Failed testing type trait is_node_group_proxy");
static_assert(not Parallel::is_node_group_proxy<group_proxy>::value,
              "Failed testing type trait is_node_group_proxy");
static_assert(Parallel::is_node_group_proxy<nodegroup_proxy>::value,
              "Failed testing type trait is_node_group_proxy");

// [has_pup_member_example]
static_assert(Parallel::has_pup_member<PupableClass>::value,
              "Failed testing type trait has_pup_member");
static_assert(Parallel::has_pup_member_t<PupableClass>::value,
              "Failed testing type trait has_pup_member");
static_assert(Parallel::has_pup_member_v<PupableClass>,
              "Failed testing type trait has_pup_member");
static_assert(not Parallel::has_pup_member<NonpupableClass>::value,
              "Failed testing type trait has_pup_member");
// [has_pup_member_example]

// [is_pupable_example]
static_assert(Parallel::is_pupable<PupableClass>::value,
              "Failed testing type trait is_pupable");
static_assert(Parallel::is_pupable_t<PupableClass>::value,
              "Failed testing type trait is_pupable");
static_assert(Parallel::is_pupable_v<PupableClass>,
              "Failed testing type trait is_pupable");
static_assert(not Parallel::is_pupable<NonpupableClass>::value,
              "Failed testing type trait is_pupable");
// [is_pupable_example]

static_assert(
    Parallel::is_chare_type<Parallel::Algorithms::Singleton, singleton>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Singleton.");
static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Singleton, array>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Singleton.");
static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Singleton, group>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Singleton.");
static_assert(not Parallel::is_chare_type<Parallel::Algorithms::Singleton,
                                          nodegroup>::value,
              "Failed testing type trait is_chare_type for "
              "Parallel::Algorithms::Singleton.");

static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Array, singleton>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Array.");
static_assert(
    Parallel::is_chare_type<Parallel::Algorithms::Array, array>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Array.");
static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Array, group>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Array.");
static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Array, nodegroup>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Array.");

static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Group, singleton>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Group.");
static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Group, array>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Group.");
static_assert(
    Parallel::is_chare_type<Parallel::Algorithms::Group, group>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Group.");
static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Group, nodegroup>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Group.");

static_assert(not Parallel::is_chare_type<Parallel::Algorithms::Nodegroup,
                                          singleton>::value,
              "Failed testing type trait is_chare_type for "
              "Parallel::Algorithms::Nodegroup.");
static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Nodegroup, array>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Nodegroup.");
static_assert(
    not Parallel::is_chare_type<Parallel::Algorithms::Nodegroup, group>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Nodegroup.");
static_assert(
    Parallel::is_chare_type<Parallel::Algorithms::Nodegroup, nodegroup>::value,
    "Failed testing type trait is_chare_type for "
    "Parallel::Algorithms::Nodegroup.");

static_assert(Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Singleton, array_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Singleton.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Singleton, array_element_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Singleton.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Singleton, group_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Singleton.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Singleton, nodegroup_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Singleton.");

static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Array, array_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Array.");
static_assert(Parallel::is_chare_type_from_proxy<Parallel::Algorithms::Array,
                                                 array_element_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Array.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Array, group_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Array.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Array, nodegroup_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Array.");

static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Group, array_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Group.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Group, array_element_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Group.");
static_assert(Parallel::is_chare_type_from_proxy<Parallel::Algorithms::Group,
                                                 group_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Group.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Group, nodegroup_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Group.");

static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Nodegroup, array_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Nodegroup.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Nodegroup, array_element_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Nodegroup.");
static_assert(not Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Nodegroup, group_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Nodegroup.");
static_assert(Parallel::is_chare_type_from_proxy<
                  Parallel::Algorithms::Nodegroup, nodegroup_proxy>::value,
              "Failed testing type trait is_chare_type_from_proxy for "
              "Parallel::Algorithms::Nodegroup.");
