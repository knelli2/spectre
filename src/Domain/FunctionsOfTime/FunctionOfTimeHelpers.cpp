// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/FunctionsOfTime/FunctionOfTimeHelpers.hpp"

#include <atomic>
#include <list>
#include <pup.h>
#include <pup_stl.h>

#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace domain::FunctionsOfTime::FunctionOfTimeHelpers {
template <size_t MaxDerivPlusOne, bool StoreCoefs>
StoredInfo<MaxDerivPlusOne, StoreCoefs>::StoredInfo(
    const double t, std::array<DataVector, MaxDerivPlusOne> init_quantities)
    : time(t), stored_quantities{std::move(init_quantities)} {
  // Convert derivs to coefficients for polynomial evaluation
  // The coefficient of x^N is the Nth deriv rescaled by 1/factorial(N)
  double fact = 1.0;
  if (MaxDerivPlusOne > 2) {
    for (size_t j = 2; j < MaxDerivPlusOne; j++) {
      fact *= j;
      gsl::at(stored_quantities, j) /= fact;
    }
  }
}

template <size_t MaxDerivPlusOne>
StoredInfo<MaxDerivPlusOne, false>::StoredInfo(
    const double t, std::array<DataVector, MaxDerivPlusOne> init_quantities)
    : time(t), stored_quantities{std::move(init_quantities)} {}

template <size_t MaxDerivPlusOne, bool StoreCoefs>
void StoredInfo<MaxDerivPlusOne, StoreCoefs>::pup(PUP::er& p) {
  p | time;
  p | stored_quantities;
}

template <size_t MaxDerivPlusOne>
void StoredInfo<MaxDerivPlusOne, false>::pup(PUP::er& p) {
  p | time;
  p | stored_quantities;
}

void reset_expiration_time(
    const gsl::not_null<std::atomic<double>*> prev_expiration_time,
    const double next_expiration_time) {
  const double loaded_prev_expr_time =
      prev_expiration_time->load(std::memory_order_acquire);
  if (next_expiration_time < loaded_prev_expr_time) {
    ERROR("Attempted to change expiration time to "
          << next_expiration_time
          << ", which precedes the previous expiration time of "
          << loaded_prev_expr_time << ".");
  }
  prev_expiration_time->store(next_expiration_time, std::memory_order_release);
}

template <size_t MaxDerivPlusOne, bool StoreCoefs>
const StoredInfo<MaxDerivPlusOne, StoreCoefs>& stored_info_from_upper_bound(
    const double t,
    const std::list<StoredInfo<MaxDerivPlusOne, StoreCoefs>>& all_stored_infos,
    const size_t all_info_size) {
  // this function assumes that the times in stored_info_at_update_times is
  // sorted, which is enforced by the update function of a piecewise polynomial.

  // We check the passed in parameter rather than the actual size of
  // all_stored_infos because another thread may have added an entry to
  // all_stored_infos during a call to this function.
  ASSERT(all_info_size != 0,
         "List of StoredInfos you are trying to access is empty. Was it "
         "constructed properly?");

  auto upper_bound_stored_info = all_stored_infos.begin();
  // Linear search is faster with std::list because of expensive
  // std::advance/distance inside std::lower_bound
  for (size_t i = 0; i < all_info_size; i++, upper_bound_stored_info++) {
    if (not(upper_bound_stored_info->time < t)) {
      break;
    }
  }

  if (upper_bound_stored_info == all_stored_infos.begin()) {
    // all elements of times are greater than t
    // check if t is just less than the min element by roundoff
    if (not equal_within_roundoff(upper_bound_stored_info->time, t)) {
      ERROR("requested time " << t << " precedes earliest time "
                              << upper_bound_stored_info->time << " of times.");
    }
    return *upper_bound_stored_info;
  }
  // t is either greater than all elements of times
  // or t is within the range of times.
  // In both cases, 'upper_bound_stored_info' currently points to one index past
  // the desired index.
  return *std::prev(upper_bound_stored_info);
}

template <size_t MaxDerivPlusOne, bool StoreCoefs>
bool operator==(
    const domain::FunctionsOfTime::FunctionOfTimeHelpers::StoredInfo<
        MaxDerivPlusOne, StoreCoefs>& lhs,
    const domain::FunctionsOfTime::FunctionOfTimeHelpers::StoredInfo<
        MaxDerivPlusOne, StoreCoefs>& rhs) {
  return lhs.time == rhs.time and
         lhs.stored_quantities == rhs.stored_quantities;
}

template <size_t MaxDerivPlusOne, bool StoreCoefs>
bool operator!=(
    const domain::FunctionsOfTime::FunctionOfTimeHelpers::StoredInfo<
        MaxDerivPlusOne, StoreCoefs>& lhs,
    const domain::FunctionsOfTime::FunctionOfTimeHelpers::StoredInfo<
        MaxDerivPlusOne, StoreCoefs>& rhs) {
  return not(lhs == rhs);
}

template <size_t MaxDerivPlusOne, bool StoreCoefs>
std::ostream& operator<<(std::ostream& os,
                         const StoredInfo<MaxDerivPlusOne, StoreCoefs>& info) {
  os << "t=" << info.time << ": ";
  for (size_t i = 0; i < MaxDerivPlusOne - 1; ++i) {
    os << gsl::at(info.stored_quantities, i) << " ";
  }
  os << info.stored_quantities[MaxDerivPlusOne - 1];
  return os;
}

// explicit instantiation of StoredInfo class and stored_info_from_upper_bound
// function for MaxDerivPlusOne = {1,2,3,4,5}
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define STORECOEF(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data)                             \
  template class StoredInfo<DIM(data), STORECOEF(data)>; \
  template bool operator==<DIM(data), STORECOEF(data)>(  \
      const StoredInfo<DIM(data), STORECOEF(data)>&,     \
      const StoredInfo<DIM(data), STORECOEF(data)>&);    \
  template bool operator!=<DIM(data), STORECOEF(data)>(  \
      const StoredInfo<DIM(data), STORECOEF(data)>&,     \
      const StoredInfo<DIM(data), STORECOEF(data)>&);    \
  template std::ostream& operator<<(                     \
      std::ostream&, const StoredInfo<DIM(data), STORECOEF(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3, 4, 5), (true, false))

#undef INSTANTIATE

#define INSTANTIATE(_, data)                                                  \
  template const StoredInfo<DIM(data), STORECOEF(data)>&                      \
  stored_info_from_upper_bound(                                               \
      const double, const std::list<StoredInfo<DIM(data), STORECOEF(data)>>&, \
      const size_t);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3, 4, 5), (true, false))

#undef INSTANTIATE
#undef STORECOEF
#undef DIM
}  // namespace domain::FunctionsOfTime::FunctionOfTimeHelpers
