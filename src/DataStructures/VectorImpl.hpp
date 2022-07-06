// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>  // IWYU pragma: keep  // for std::fill
#include <array>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/CustomVector.h>
#include <blaze/math/DenseVector.h>
#include <blaze/math/GroupTag.h>
#include <blaze/math/PaddingFlag.h>
#include <blaze/math/TransposeFlag.h>
#include <cstddef>
#include <cstring>
#include <functional>  // IWYU pragma: keep  // for std::plus, etc.
#include <initializer_list>
#include <limits>
#include <memory>
#include <ostream>
#include <pup.h>
#include <type_traits>

#include "DataStructures/Blaze/StepFunction.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ForceInline.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/MakeWithValue.hpp"  // IWYU pragma: keep
#include "Utilities/MemoryHelpers.hpp"
#include "Utilities/PrintHelpers.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TypeTraits/IsComplexOfFundamental.hpp"

/*!
 * \ingroup DataStructuresGroup
 * \brief Base class template for various DataVector and related types
 *
 * \details The `VectorImpl` class is the generic parent class for vectors
 * representing collections of related function values, such as `DataVector`s
 * for contiguous data over a computational domain.
 *
 * The `VectorImpl` does not itself define any particular mathematical
 * operations on the contained values. The `VectorImpl` template class and the
 * macros defined in `VectorImpl.hpp` assist in the construction of various
 * derived classes supporting a chosen set of mathematical operations.
 *
 * In addition, the equivalence operator `==` is inherited from the underlying
 * `blaze::CustomVector` type, and returns true if and only if the size and
 * contents of the two compared vectors are equivalent.
 *
 * Template parameters:
 * - `T` is the underlying stored type, e.g. `double`, `std::complex<double>`,
 *   `float`, etc.
 * - `VectorType` is the type that should be associated with the VectorImpl
 *    during mathematical computations. In most cases, inherited types should
 *    have themselves as the second template argument, e.g.
 *  ```
 *  class DataVector : VectorImpl<double, DataVector> {
 *  ```
 *  The second template parameter communicates arithmetic type restrictions to
 *  the underlying Blaze framework. For example, if `VectorType` is
 *  `DataVector`, then the underlying architecture will prevent addition with a
 *  vector type whose `ResultType` (which is aliased to its `VectorType`) is
 *  `ModalVector`.  Since `DataVector`s and `ModalVector`s represent data in
 *  different spaces, we wish to forbid several operations between them. This
 *  vector-type-tracking through an expression prevents accidental mixing of
 *  vector types in math expressions.
 *
 * \note
 * - If created with size 0, then  `data()` will return `nullptr`
 * - If either `SPECTRE_DEBUG` or `SPECTRE_NAN_INIT` are defined, then the
 *   `VectorImpl` is default initialized to `signaling_NaN()`. Otherwise, the
 *   vector is filled with uninitialized memory for performance.
 */
template <typename T, typename VectorType>
class VectorImpl
    : public blaze::CustomVector<
          T, blaze::AlignmentFlag::unaligned, blaze::PaddingFlag::unpadded,
          blaze::defaultTransposeFlag, blaze::GroupTag<0>, VectorType> {
 public:
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using BaseType = blaze::CustomVector<
      T, blaze::AlignmentFlag::unaligned, blaze::PaddingFlag::unpadded,
      blaze::defaultTransposeFlag, blaze::GroupTag<0>, VectorType>;
  static constexpr bool transpose_flag = blaze::defaultTransposeFlag;

  using ElementType = T;
  using TransposeType = VectorImpl<T, VectorType>;
  using CompositeType = const VectorImpl<T, VectorType>&;
  using iterator = typename BaseType::Iterator;
  using const_iterator = typename BaseType::ConstIterator;

  using BaseType::operator[];
  using BaseType::begin;
  using BaseType::cbegin;
  using BaseType::cend;
  using BaseType::data;
  using BaseType::end;
  using BaseType::size;

  /// @{
  /// Upcast to `BaseType`
  /// \attention
  /// upcast should only be used when implementing a derived vector type, not in
  /// calling code
  const BaseType& operator*() const {
    return static_cast<const BaseType&>(*this);
  }
  BaseType& operator*() { return static_cast<BaseType&>(*this); }
  /// @}

  /// Create with the given size. In debug mode, the vector is initialized to
  /// 'NaN' by default. If not initialized to 'NaN', the memory is allocated but
  /// not initialized.
  ///
  /// - `set_size` number of values
  explicit VectorImpl(size_t set_size)
      : owned_data_(
            set_size > 0
                ? cpp20::make_unique_for_overwrite<value_type[]>(set_size)
                : nullptr) {
#if defined(SPECTRE_DEBUG) || defined(SPECTRE_NAN_INIT)
    std::fill(owned_data_.get(), owned_data_.get() + set_size,
              std::numeric_limits<value_type>::signaling_NaN());
#endif  // SPECTRE_DEBUG
    reset_pointer_vector(set_size);
  }

  /// Create with the given size and value.
  ///
  /// - `set_size` number of values
  /// - `value` the value to initialize each element
  VectorImpl(size_t set_size, T value)
      : owned_data_(
            set_size > 0
                ? cpp20::make_unique_for_overwrite<value_type[]>(set_size)
                : nullptr) {
    std::fill(owned_data_.get(), owned_data_.get() + set_size, value);
    reset_pointer_vector(set_size);
  }

  /// Create a non-owning VectorImpl that points to `start`
  VectorImpl(T* start, size_t set_size)
      : BaseType(start, set_size), owning_(false) {}

  /// Create from an initializer list of `T`.
  template <class U, Requires<std::is_same_v<U, T>> = nullptr>
  VectorImpl(std::initializer_list<U> list)
      : owned_data_(
            list.size() > 0
                ? cpp20::make_unique_for_overwrite<value_type[]>(list.size())
                : nullptr) {
    // Note: can't use memcpy with an initializer list.
    std::copy(list.begin(), list.end(), owned_data_.get());
    reset_pointer_vector(list.size());
  }

  /// Empty VectorImpl
  VectorImpl() = default;
  /// \cond HIDDEN_SYMBOLS
  ~VectorImpl() = default;

  VectorImpl(const VectorImpl<T, VectorType>& rhs);
  VectorImpl& operator=(const VectorImpl<T, VectorType>& rhs);
  VectorImpl(VectorImpl<T, VectorType>&& rhs);
  VectorImpl& operator=(VectorImpl<T, VectorType>&& rhs);

  // This is a converting constructor. clang-tidy complains that it's not
  // explicit, but we want it to allow conversion.
  // clang-tidy: mark as explicit (we want conversion to VectorImpl type)
  template <
      typename VT, bool VF,
      Requires<std::is_same_v<typename VT::ResultType, VectorType>> = nullptr>
  VectorImpl(const blaze::DenseVector<VT, VF>& expression);  // NOLINT

  template <typename VT, bool VF>
  VectorImpl& operator=(const blaze::DenseVector<VT, VF>& expression);
  /// \endcond

  VectorImpl& operator=(const T& rhs);

  decltype(auto) SPECTRE_ALWAYS_INLINE operator[](const size_t index) {
    ASSERT(index < size(), "Out-of-range access to element "
                               << index << " of a size " << size()
                               << " Blaze vector.");
    return BaseType::operator[](index);
  }

  decltype(auto) SPECTRE_ALWAYS_INLINE operator[](const size_t index) const {
    ASSERT(index < size(), "Out-of-range access to element "
                               << index << " of a size " << size()
                               << " Blaze vector.");
    return BaseType::operator[](index);
  }

  /// @{
  /// Set the VectorImpl to be a reference to another VectorImpl object
  void set_data_ref(gsl::not_null<VectorType*> rhs) {
    set_data_ref(rhs->data(), rhs->size());
  }

  void set_data_ref(T* const start, const size_t set_size) {
    owned_data_.reset();
    if (start == nullptr) {
      (**this).reset();
    } else {
      (**this).reset(start, set_size);
    }
    owning_ = false;
  }
  /// @}

  /*!
   * \brief A common operation for checking the size and resizing a memory
   * buffer if needed to ensure that it has the desired size. This operation is
   * not permitted on a non-owning vector.
   *
   * \note This utility should NOT be used when it is anticipated that the
   *   supplied buffer will typically be the wrong size (in that case, suggest
   *   either manual checking or restructuring so that resizing is less common).
   *   This uses `UNLIKELY` to perform the check most quickly when the buffer
   *   needs no resizing, but will be slower when resizing is common.
   */
  void SPECTRE_ALWAYS_INLINE destructive_resize(const size_t new_size) {
    if (UNLIKELY(size() != new_size)) {
      ASSERT(owning_,
             MakeString{}
                 << "Attempting to resize a non-owning vector from size: "
                 << size() << " to size: " << new_size
                 << " but we may not destructively resize a non-owning vector");
      owned_data_ = cpp20::make_unique_for_overwrite<value_type[]>(new_size);
      reset_pointer_vector(new_size);
    }
  }

  /// Returns true if the class owns the data
  bool is_owning() const { return owning_; }

  /// Serialization for Charm++
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 protected:
  std::unique_ptr<value_type[]> owned_data_{};
  bool owning_{true};

  SPECTRE_ALWAYS_INLINE void reset_pointer_vector(const size_t set_size) {
    if (set_size == 0) {
      return;
    }
    if (owned_data_ == nullptr) {
      ERROR(
          "VectorImpl::reset_pointer_vector cannot be called when owned_data_ "
          "is nullptr.");
    }
    this->reset(owned_data_.get(), set_size);
  }
};

template <typename T, typename VectorType>
VectorImpl<T, VectorType>::VectorImpl(const VectorImpl<T, VectorType>& rhs)
    : BaseType{rhs},
      owned_data_(
          rhs.size() > 0
              ? cpp20::make_unique_for_overwrite<value_type[]>(rhs.size())
              : nullptr) {
  reset_pointer_vector(rhs.size());
  std::memcpy(data(), rhs.data(), size() * sizeof(value_type));
}

template <typename T, typename VectorType>
VectorImpl<T, VectorType>& VectorImpl<T, VectorType>::operator=(
    const VectorImpl<T, VectorType>& rhs) {
  if (this != &rhs) {
    if (owning_) {
      if (size() != rhs.size()) {
        owned_data_.reset();
        if (rhs.size() > 0) {
          owned_data_ =
              cpp20::make_unique_for_overwrite<value_type[]>(rhs.size());
        }
      }
      reset_pointer_vector(rhs.size());
    } else {
      ASSERT(rhs.size() == size(), "Must copy into same size, not "
                                       << rhs.size() << " into " << size());
    }
    std::memcpy(data(), rhs.data(), size() * sizeof(value_type));
  }
  return *this;
}

template <typename T, typename VectorType>
VectorImpl<T, VectorType>::VectorImpl(VectorImpl<T, VectorType>&& rhs) {
  owned_data_ = std::move(rhs.owned_data_);
  **this = std::move(*rhs);
  owning_ = rhs.owning_;
  rhs.owning_ = true;
  rhs.reset();
}

template <typename T, typename VectorType>
VectorImpl<T, VectorType>& VectorImpl<T, VectorType>::operator=(
    VectorImpl<T, VectorType>&& rhs) {
  if (this != &rhs) {
    if (owning_) {
      owned_data_ = std::move(rhs.owned_data_);
      **this = std::move(*rhs);
      owning_ = rhs.owning_;
    } else {
      ASSERT(rhs.size() == size(), "Must copy into same size, not "
                                       << rhs.size() << " into " << size());
      std::memcpy(data(), rhs.data(), size() * sizeof(value_type));
    }
    rhs.owning_ = true;
    rhs.reset();
  }
  return *this;
}

/// \cond HIDDEN_SYMBOLS
// This is a converting constructor. clang-tidy complains that it's not
// explicit, but we want it to allow conversion.
// clang-tidy: mark as explicit (we want conversion to VectorImpl)
template <typename T, typename VectorType>
template <typename VT, bool VF,
          Requires<std::is_same_v<typename VT::ResultType, VectorType>>>
VectorImpl<T, VectorType>::VectorImpl(
    const blaze::DenseVector<VT, VF>& expression)  // NOLINT
    : owned_data_(cpp20::make_unique_for_overwrite<value_type[]>(
          (*expression).size())) {
  static_assert(std::is_same_v<typename VT::ResultType, VectorType>,
                "You are attempting to assign the result of an expression "
                "that is not consistent with the VectorImpl type you are "
                "assigning to.");
  reset_pointer_vector((*expression).size());
  **this = expression;
}

template <typename T, typename VectorType>
template <typename VT, bool VF>
VectorImpl<T, VectorType>& VectorImpl<T, VectorType>::operator=(
    const blaze::DenseVector<VT, VF>& expression) {
  static_assert(std::is_same_v<typename VT::ResultType, VectorType>,
                "You are attempting to assign the result of an expression "
                "that is not consistent with the VectorImpl type you are "
                "assigning to.");
  if (owning_ and (*expression).size() != size()) {
    owned_data_ =
        cpp20::make_unique_for_overwrite<value_type[]>((*expression).size());
    reset_pointer_vector((*expression).size());
  } else if (not owning_) {
    ASSERT((*expression).size() == size(), "Must copy into same size, not "
                                               << (*expression).size()
                                               << " into " << size());
  }
  **this = expression;
  return *this;
}
/// \endcond

// The case of assigning a type apart from the same VectorImpl or a
// `blaze::DenseVector` forwards the assignment to the `blaze::CustomVector`
// base type. In the case of a single compatible value, this fills the vector
// with that value.
template <typename T, typename VectorType>
VectorImpl<T, VectorType>& VectorImpl<T, VectorType>::operator=(const T& rhs) {
  **this = rhs;
  return *this;
}

template <typename T, typename VectorType>
void VectorImpl<T, VectorType>::pup(PUP::er& p) {  // NOLINT
  ASSERT(owning_ or p.isSizing(),
         "Cannot pup a non-owning vector unless you are sizing!");
  auto my_size = size();
  p | my_size;
  if (my_size > 0) {
    if (p.isUnpacking()) {
      owning_ = true;
      owned_data_ = cpp20::make_unique_for_overwrite<value_type[]>(my_size);
      reset_pointer_vector(my_size);
    }
    PUParray(p, data(), size());
  }
}

/// Output operator for VectorImpl
template <typename T, typename VectorType>
std::ostream& operator<<(std::ostream& os, const VectorImpl<T, VectorType>& d) {
  sequence_print_helper(os, d.begin(), d.end());
  return os;
}

#define DECLARE_GENERAL_VECTOR_BLAZE_TRAITS(VECTOR_TYPE)         \
  template <>                                                    \
  struct IsDenseVector<VECTOR_TYPE> : public blaze::TrueType {}; \
                                                                 \
  template <>                                                    \
  struct IsVector<VECTOR_TYPE> : public blaze::TrueType {};      \
                                                                 \
  template <>                                                    \
  struct CustomTransposeType<VECTOR_TYPE> {                      \
    using Type = VECTOR_TYPE;                                    \
  }

/*!
 * \ingroup DataStructuresGroup
 * \brief Instructs Blaze to provide the appropriate vector result type after
 * math operations. This is accomplished by specializing Blaze's type traits
 * that are used for handling return type deduction and specifying the `using
 * Type =` nested type alias in the traits.
 *
 * \param VECTOR_TYPE The vector type, which matches the type of the operation
 * result (e.g. `DataVector`)
 *
 * \param BLAZE_MATH_TRAIT The blaze trait/expression for which you want to
 * specify the return type (e.g. `AddTrait`).
 */
#define BLAZE_TRAIT_SPECIALIZE_BINARY_TRAIT(VECTOR_TYPE, BLAZE_MATH_TRAIT) \
  template <>                                                              \
  struct BLAZE_MATH_TRAIT<VECTOR_TYPE, VECTOR_TYPE> {                      \
    using Type = VECTOR_TYPE;                                              \
  };                                                                       \
  template <>                                                              \
  struct BLAZE_MATH_TRAIT<VECTOR_TYPE, VECTOR_TYPE::value_type> {          \
    using Type = VECTOR_TYPE;                                              \
  };                                                                       \
  template <>                                                              \
  struct BLAZE_MATH_TRAIT<VECTOR_TYPE::value_type, VECTOR_TYPE> {          \
    using Type = VECTOR_TYPE;                                              \
  }

/*!
 * \ingroup DataStructuresGroup
 * \brief Instructs Blaze to provide the appropriate vector result type of an
 * operator between `VECTOR_TYPE` and `COMPATIBLE`, where the operation is
 * represented by `BLAZE_MATH_TRAIT`
 *
 * \param VECTOR_TYPE The vector type, which matches the type of the operation
 * result (e.g. `ComplexDataVector`)
 *
 * \param COMPATIBLE the type for which you want math operations to work with
 * `VECTOR_TYPE` smoothly (e.g. `DataVector`)
 *
 * \param BLAZE_MATH_TRAIT The blaze trait for which you want declare the Type
 * field (e.g. `AddTrait`)
 *
 * \param RESULT_TYPE The type which should be used as the 'return' type for the
 * binary operation
 */
#define BLAZE_TRAIT_SPECIALIZE_COMPATIBLE_BINARY_TRAIT(     \
    VECTOR_TYPE, COMPATIBLE, BLAZE_MATH_TRAIT, RESULT_TYPE) \
  template <>                                               \
  struct BLAZE_MATH_TRAIT<VECTOR_TYPE, COMPATIBLE> {        \
    using Type = RESULT_TYPE;                               \
  };                                                        \
  template <>                                               \
  struct BLAZE_MATH_TRAIT<COMPATIBLE, VECTOR_TYPE> {        \
    using Type = RESULT_TYPE;                               \
  }

/*!
 * \ingroup DataStructuresGroup
 * \brief Instructs Blaze to provide the appropriate vector result type of
 * arithmetic operations for `VECTOR_TYPE`. This is accomplished by specializing
 * Blaze's type traits that are used for handling return type deduction.
 *
 * \details Type definitions here are suitable for contiguous data
 * (e.g. `DataVector`), but this macro might need to be tweaked for other types
 * of data, for instance Fourier coefficients.
 *
 * \param VECTOR_TYPE The vector type, which for the arithmetic operations is
 * the type of the operation result (e.g. `DataVector`)
 */
#define VECTOR_BLAZE_TRAIT_SPECIALIZE_ARITHMETIC_TRAITS(VECTOR_TYPE) \
  template <>                                                        \
  struct TransposeFlag<VECTOR_TYPE>                                  \
      : BoolConstant<VECTOR_TYPE::transpose_flag> {};                \
  BLAZE_TRAIT_SPECIALIZE_BINARY_TRAIT(VECTOR_TYPE, AddTrait);        \
  BLAZE_TRAIT_SPECIALIZE_BINARY_TRAIT(VECTOR_TYPE, SubTrait);        \
  BLAZE_TRAIT_SPECIALIZE_BINARY_TRAIT(VECTOR_TYPE, MultTrait);       \
  BLAZE_TRAIT_SPECIALIZE_BINARY_TRAIT(VECTOR_TYPE, DivTrait)

/*!
 * \ingroup DataStructuresGroup
 * \brief Instructs Blaze to provide the appropriate vector result type of `Map`
 * operations (unary and binary) acting on `VECTOR_TYPE`. This is accomplished
 * by specializing Blaze's type traits that are used for handling return type
 * deduction.
 *
 * \details Type declarations here are suitable for contiguous data (e.g.
 * `DataVector`), but this macro might need to be tweaked for other types of
 * data, for instance Fourier coefficients.
 *
 * \param VECTOR_TYPE The vector type, which for the `Map` operations is
 * the type of the operation result (e.g. `DataVector`)
 */
#define VECTOR_BLAZE_TRAIT_SPECIALIZE_ALL_MAP_TRAITS(VECTOR_TYPE) \
  template <typename Operator>                                    \
  struct MapTrait<VECTOR_TYPE, Operator> {                        \
    using Type = VECTOR_TYPE;                                     \
  };                                                              \
  template <typename Operator>                                    \
  struct MapTrait<VECTOR_TYPE, VECTOR_TYPE, Operator> {           \
    using Type = VECTOR_TYPE;                                     \
  }

/*!
 * \ingroup DataStructuresGroup
 * \brief Defines the set of binary operations often supported for
 * `std::array<VECTOR_TYPE, size>`, for arbitrary `size`.
 *
 *  \param VECTOR_TYPE The vector type (e.g. `DataVector`)
 */
#define MAKE_STD_ARRAY_VECTOR_BINOPS(VECTOR_TYPE)                            \
  DEFINE_STD_ARRAY_BINOP(VECTOR_TYPE, VECTOR_TYPE::value_type,               \
                         VECTOR_TYPE, operator+, std::plus<>())              \
  DEFINE_STD_ARRAY_BINOP(VECTOR_TYPE, VECTOR_TYPE,                           \
                         VECTOR_TYPE::value_type, operator+, std::plus<>())  \
  DEFINE_STD_ARRAY_BINOP(VECTOR_TYPE, VECTOR_TYPE, VECTOR_TYPE, operator+,   \
                         std::plus<>())                                      \
                                                                             \
  DEFINE_STD_ARRAY_BINOP(VECTOR_TYPE, VECTOR_TYPE::value_type,               \
                         VECTOR_TYPE, operator-, std::minus<>())             \
  DEFINE_STD_ARRAY_BINOP(VECTOR_TYPE, VECTOR_TYPE,                           \
                         VECTOR_TYPE::value_type, operator-, std::minus<>()) \
  DEFINE_STD_ARRAY_BINOP(VECTOR_TYPE, VECTOR_TYPE, VECTOR_TYPE, operator-,   \
                         std::minus<>())                                     \
                                                                             \
  DEFINE_STD_ARRAY_INPLACE_BINOP(VECTOR_TYPE, VECTOR_TYPE, operator-=,       \
                                 std::minus<>())                             \
  DEFINE_STD_ARRAY_INPLACE_BINOP(                                            \
      VECTOR_TYPE, VECTOR_TYPE::value_type, operator-=, std::minus<>())      \
  DEFINE_STD_ARRAY_INPLACE_BINOP(VECTOR_TYPE, VECTOR_TYPE, operator+=,       \
                                 std::plus<>())                              \
  DEFINE_STD_ARRAY_INPLACE_BINOP(                                            \
      VECTOR_TYPE, VECTOR_TYPE::value_type, operator+=, std::plus<>())

/*!
 * \ingroup DataStructuresGroup
 * \brief Defines the `MakeWithValueImpl` `apply` specialization
 *
 * \details The `MakeWithValueImpl<VECTOR_TYPE, VECTOR_TYPE>` member
 * `apply(VECTOR_TYPE, VECTOR_TYPE::value_type)` specialization defined by this
 * macro produces an object with the same size as the `input` argument,
 * initialized with the `value` argument in every entry.
 *
 * \param VECTOR_TYPE The vector type (e.g. `DataVector`)
 */
#define MAKE_WITH_VALUE_IMPL_DEFINITION_FOR(VECTOR_TYPE)                  \
  namespace MakeWithValueImpls {                                          \
  template <>                                                             \
  struct NumberOfPoints<VECTOR_TYPE> {                                    \
    static SPECTRE_ALWAYS_INLINE size_t apply(const VECTOR_TYPE& input) { \
      return input.size();                                                \
    }                                                                     \
  };                                                                      \
  template <>                                                             \
  struct MakeWithSize<VECTOR_TYPE> {                                      \
    static SPECTRE_ALWAYS_INLINE VECTOR_TYPE                              \
    apply(const size_t size, const VECTOR_TYPE::value_type value) {       \
      return VECTOR_TYPE(size, value);                                    \
    }                                                                     \
  };                                                                      \
  }  // namespace MakeWithValueImpls

/// @{
/*!
 * \ingroup DataStructuresGroup
 * \ingroup TypeTraitsGroup
 * \brief Helper struct to determine the element type of a VectorImpl or
 * container of VectorImpl
 *
 * \details Extracts the element type of a `VectorImpl`, a std::array of
 * `VectorImpl`, or a reference or pointer to a `VectorImpl`. In any of these
 * cases, the `type` member is defined as the `ElementType` of the `VectorImpl`
 * in question. If, instead, `get_vector_element_type` is passed an arithmetic
 * or complex arithemetic type, the `type` member is defined as the passed type.
 *
 * \snippet DataStructures/Test_VectorImpl.cpp get_vector_element_type_example
 */
// cast to bool needed to avoid the compiler mistaking the type to be determined
// by T
template <typename T,
          bool = static_cast<bool>(tt::is_complex_of_fundamental_v<T> or
                                   std::is_fundamental_v<T>)>
struct get_vector_element_type;
template <typename T>
struct get_vector_element_type<T, true> {
  using type = T;
};
template <typename T>
struct get_vector_element_type<const T, false> {
  using type = typename get_vector_element_type<T>::type;
};
template <typename T>
struct get_vector_element_type<T, false> {
  using type = typename get_vector_element_type<
      typename T::ResultType::ElementType>::type;
};
template <typename T>
struct get_vector_element_type<T*, false> {
  using type = typename get_vector_element_type<T>::type;
};
template <typename T>
struct get_vector_element_type<T&, false> {
  using type = typename get_vector_element_type<T>::type;
};
template <typename T, size_t S>
struct get_vector_element_type<std::array<T, S>, false> {
  using type = typename get_vector_element_type<T>::type;
};
/// @}

template <typename T>
using get_vector_element_type_t = typename get_vector_element_type<T>::type;

namespace detail {
template <typename... VectorImplTemplateArgs>
std::true_type is_derived_of_vector_impl_impl(
    const VectorImpl<VectorImplTemplateArgs...>*);

std::false_type is_derived_of_vector_impl_impl(...);
}  // namespace detail

/// \ingroup TypeTraitsGroup
/// This is `std::true_type` if the provided type possesses an implicit
/// conversion to any `VectorImpl`, which is the primary feature of SpECTRE
/// vectors generally. Otherwise, it is `std::false_type`.
template <typename T>
using is_derived_of_vector_impl =
    decltype(detail::is_derived_of_vector_impl_impl(std::declval<T*>()));

template <typename T>
constexpr bool is_derived_of_vector_impl_v =
    is_derived_of_vector_impl<T>::value;

// impose strict equality for derived classes of VectorImpl; note that this
// overrides intended behavior in blaze for comparison operators to use
// approximate equality in favor of equality between containers being
// appropriately recursive. This form primarily works by using templates to
// ensure that our comparison operator is resolved with higher priority than the
// blaze form as of blaze 3.8
template <
    typename Lhs, typename Rhs,
    Requires<(is_derived_of_vector_impl_v<Lhs> or
              is_derived_of_vector_impl_v<
                  Rhs>)and not(std::is_base_of_v<blaze::Computation, Lhs> or
                               std::is_base_of_v<blaze::Computation, Rhs>) and
             not(std::is_same_v<Rhs, typename Lhs::ElementType> or
                 std::is_same_v<Lhs, typename Rhs::ElementType>)> = nullptr>
bool operator==(const Lhs& lhs, const Rhs& rhs) {
  return blaze::equal<blaze::strict>(lhs, rhs);
}

template <
    typename Lhs, typename Rhs,
    Requires<(is_derived_of_vector_impl_v<Lhs> or
              is_derived_of_vector_impl_v<
                  Rhs>)and not(std::is_base_of_v<blaze::Computation, Lhs> or
                               std::is_base_of_v<blaze::Computation, Rhs>) and
             not(std::is_same_v<Rhs, typename Lhs::ElementType> or
                 std::is_same_v<Lhs, typename Lhs::ElementType>)> = nullptr>
bool operator!=(const Lhs& lhs, const Rhs& rhs) {
  return not(lhs == rhs);
}

// Impose strict equality for any expression templates; note that
// this overrides intended behavior in blaze for comparison
// operators to use approximate equality in favor of equality
// between containers being appropriately recursive. This form
// primarily works by using templates to ensure that our
// comparison operator is resolved with higher priority than the
// blaze form as of blaze 3.8
template <typename Lhs, typename Rhs,
          Requires<std::is_base_of_v<blaze::Computation, Lhs> or
                   std::is_base_of_v<blaze::Computation, Rhs>> = nullptr>
bool operator==(const Lhs& lhs, const Rhs& rhs) {
  return blaze::equal<blaze::strict>(lhs, rhs);
}

template <typename Lhs, typename Rhs,
          Requires<std::is_base_of_v<blaze::Computation, Lhs> or
                   std::is_base_of_v<blaze::Computation, Rhs>> = nullptr>
bool operator!=(const Lhs& lhs, const Rhs& rhs) {
  return not(lhs == rhs);
}

template <typename Lhs, Requires<is_derived_of_vector_impl_v<Lhs>> = nullptr>
bool operator==(const Lhs& lhs, const typename Lhs::ElementType& rhs) {
  for (const auto& element : lhs) {
    if (element != rhs) {
      return false;
    }
  }
  return true;
}

template <typename Lhs, Requires<is_derived_of_vector_impl_v<Lhs>> = nullptr>
bool operator!=(const Lhs& lhs, const typename Lhs::ElementType& rhs) {
  return not(lhs == rhs);
}

template <typename Rhs, Requires<is_derived_of_vector_impl_v<Rhs>> = nullptr>
bool operator==(const typename Rhs::ElementType& lhs, const Rhs& rhs) {
  return rhs == lhs;
}

template <typename Rhs, Requires<is_derived_of_vector_impl_v<Rhs>> = nullptr>
bool operator!=(const typename Rhs::ElementType& lhs, const Rhs& rhs) {
  return not(lhs == rhs);
}

/// \ingroup DataStructuresGroup
/// Make the input `view` a `const` view of the const data `vector`, at
/// offset `offset` and length `extent`.
///
/// \warning This DOES modify the (const) input `view`. The reason `view` is
/// taken by const pointer is to try to insist that the object to be a `const`
/// view is actually const. Of course, there are ways of subverting this
/// intended functionality and editing the data pointed into by `view` after
/// this function is called; doing so is highly discouraged and results in
/// undefined behavior.
template <typename VectorType,
          Requires<is_derived_of_vector_impl_v<VectorType>> = nullptr>
void make_const_view(const gsl::not_null<const VectorType*> view,
                     const VectorType& vector, const size_t offset,
                     const size_t extent) {
  const_cast<VectorType*>(view.get())  // NOLINT
      ->set_data_ref(
          const_cast<typename VectorType::value_type*>(vector.data())  // NOLINT
              + offset,                                                // NOLINT
          extent);
}
