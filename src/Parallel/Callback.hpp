// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines Parallel::Callback.

#pragma once

#include <iomanip>
#include <ios>
#include <optional>
#include <pup.h>
#include <sstream>
#include <tuple>
#include <utility>

#include "Parallel/Invoke.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/StlStreamDeclarations.hpp"
#include "Utilities/TypeTraits/IsStreamable.hpp"

namespace Parallel {
namespace detail {
template <typename... Args>
std::string stream_args(const std::tuple<Args...>& args) {
  std::stringstream ss{};
  ss << std::setprecision(std::numeric_limits<double>::digits10 + 4)
     << std::scientific;
  ss << "(";

  tmpl::for_each<tmpl::make_sequence<tmpl::size_t<0>,
                                     tmpl::size<tmpl::list<Args...>>::value>>(
      [&](const auto index_v) {
        constexpr size_t index = tmpl::type_from<decltype(index_v)>::value;

        if constexpr (tt::is_streamable_v<std::ostream,
                                          decltype(std::get<index>(args))>) {
          using ::operator<<;
          ss << std::get<index>(args) << ", ";
        }
      });

  ss << ")";

  return ss.str();
}
}  // namespace detail

/// An abstract base class, whose derived class holds a function that
/// can be invoked at a later time.  The function is intended to be
/// invoked only once.
class Callback : public PUP::able {
 public:
  WRAPPED_PUPable_abstract(Callback);  // NOLINT
  Callback() = default;
  Callback(const Callback&) = default;
  Callback& operator=(const Callback&) = default;
  Callback(Callback&&) = default;
  Callback& operator=(Callback&&) = default;
  ~Callback() override = default;
  explicit Callback(CkMigrateMessage* msg) : PUP::able(msg) {}
  virtual void invoke() = 0;
  virtual void register_with_charm() = 0;
  virtual std::optional<bool> is_equal_to(const Callback& rhs) const = 0;
  virtual std::string name() const = 0;
};

/// Wraps a call to a simple action and its arguments.
/// Can be invoked only once.
template <typename SimpleAction, typename Proxy, typename... Args>
class SimpleActionCallback : public Callback {
 public:
  WRAPPED_PUPable_decl_template(SimpleActionCallback);  // NOLINT
  SimpleActionCallback() = default;
  // NOLINTNEXTLINE(google-explicit-constructor)
  SimpleActionCallback(Proxy proxy, std::decay_t<Args>... args)
      : proxy_(proxy), args_(std::move(args)...) {}
  explicit SimpleActionCallback(CkMigrateMessage* msg) : Callback(msg) {}
  using PUP::able::register_constructor;
  void invoke() override {
    std::apply(
        [this](auto&&... args) {
          Parallel::simple_action<SimpleAction>(proxy_, args...);
        },
        std::move(args_));
  }
  void pup(PUP::er& p) override {
    p | proxy_;
    p | args_;
  }

  void register_with_charm() override {
    static bool done_registration{false};
    if (done_registration) {
      return;
    }
    done_registration = true;
    register_classes_with_charm<SimpleActionCallback>();
  }

  std::optional<bool> is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr = dynamic_cast<const SimpleActionCallback*>(&rhs);
    if (downcast_ptr == nullptr) {
      return std::nullopt;
    }
    return  // proxy_ == downcast_ptr->proxy_ and
        {args_ == downcast_ptr->args_};
  }

  std::string name() const override {
    return "SimpleActionCallback(" + pretty_type::get_name<SimpleAction>() +
           "," + pretty_type::name<Proxy>() +
           ",Args=" + detail::stream_args(args_) + ")";
  }

 private:
  std::decay_t<Proxy> proxy_{};
  std::tuple<std::decay_t<Args>...> args_{};
};

/// Wraps a call to a simple action without arguments.
template <typename SimpleAction, typename Proxy>
class SimpleActionCallback<SimpleAction, Proxy> : public Callback {
 public:
  WRAPPED_PUPable_decl_template(SimpleActionCallback);  // NOLINT
  SimpleActionCallback() = default;
  // NOLINTNEXTLINE(google-explicit-constructor)
  SimpleActionCallback(Proxy proxy) : proxy_(proxy) {}
  explicit SimpleActionCallback(CkMigrateMessage* msg) : Callback(msg) {}
  using PUP::able::register_constructor;
  void invoke() override { Parallel::simple_action<SimpleAction>(proxy_); }

  void pup(PUP::er& p) override { p | proxy_; }

  void register_with_charm() override {
    static bool done_registration{false};
    if (done_registration) {
      return;
    }
    done_registration = true;
    register_classes_with_charm<SimpleActionCallback>();
  }

  std::optional<bool> is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr = dynamic_cast<const SimpleActionCallback*>(&rhs);
    return downcast_ptr == nullptr ? std::nullopt : std::optional{true};
  }

  std::string name() const override {
    return "SimpleActionCallback(" + pretty_type::get_name<SimpleAction>() +
           "," + pretty_type::name<Proxy>() + ")";
  }

 private:
  std::decay_t<Proxy> proxy_{};
};

/// Wraps a call to a threaded action and its arguments.
/// Can be invoked only once.
template <typename ThreadedAction, typename Proxy, typename... Args>
class ThreadedActionCallback : public Callback {
 public:
  WRAPPED_PUPable_decl_template(ThreadedActionCallback);  // NOLINT
  ThreadedActionCallback() = default;
  // NOLINTNEXTLINE(google-explicit-constructor)
  ThreadedActionCallback(Proxy proxy, std::decay_t<Args>... args)
      : proxy_(proxy), args_(std::move(args)...) {}
  explicit ThreadedActionCallback(CkMigrateMessage* msg) : Callback(msg) {}
  using PUP::able::register_constructor;
  void invoke() override {
    std::apply(
        [this](auto&&... args) {
          Parallel::threaded_action<ThreadedAction>(proxy_, args...);
        },
        std::move(args_));
  }
  void pup(PUP::er& p) override {
    p | proxy_;
    p | args_;
  }

  void register_with_charm() override {
    static bool done_registration{false};
    if (done_registration) {
      return;
    }
    done_registration = true;
    register_classes_with_charm<ThreadedActionCallback>();
  }

  std::optional<bool> is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr =
        dynamic_cast<const ThreadedActionCallback*>(&rhs);
    if (downcast_ptr == nullptr) {
      return std::nullopt;
    }
    return  // proxy_ == downcast_ptr->proxy_ and
        {args_ == downcast_ptr->args_};
  }

  std::string name() const override {
    return "ThreadedActionCallback(" + pretty_type::get_name<ThreadedAction>() +
           "," + pretty_type::name<Proxy>() +
           ",Args=" + detail::stream_args(args_) + ")";
  }

 private:
  std::decay_t<Proxy> proxy_{};
  std::tuple<std::decay_t<Args>...> args_{};
};

/// Wraps a call to a threaded action without arguments.
template <typename ThreadedAction, typename Proxy>
class ThreadedActionCallback<ThreadedAction, Proxy> : public Callback {
 public:
  WRAPPED_PUPable_decl_template(ThreadedActionCallback);  // NOLINT
  ThreadedActionCallback() = default;
  // NOLINTNEXTLINE(google-explicit-constructor)
  ThreadedActionCallback(Proxy proxy) : proxy_(proxy) {}
  explicit ThreadedActionCallback(CkMigrateMessage* msg) : Callback(msg) {}
  using PUP::able::register_constructor;
  void invoke() override { Parallel::threaded_action<ThreadedAction>(proxy_); }

  void pup(PUP::er& p) override { p | proxy_; }

  void register_with_charm() override {
    static bool done_registration{false};
    if (done_registration) {
      return;
    }
    done_registration = true;
    register_classes_with_charm<ThreadedActionCallback>();
  }

  std::optional<bool> is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr =
        dynamic_cast<const ThreadedActionCallback*>(&rhs);
    return downcast_ptr == nullptr ? std::nullopt : std::optional{true};
  }

  std::string name() const override {
    return "ThreadedActionCallback(" + pretty_type::get_name<ThreadedAction>() +
           "," + pretty_type::name<Proxy>() + ")";
  }

 private:
  std::decay_t<Proxy> proxy_{};
};

/// Wraps a call to perform_algorithm.
template <typename Proxy>
class PerformAlgorithmCallback : public Callback {
 public:
  WRAPPED_PUPable_decl_template(PerformAlgorithmCallback);  // NOLINT
  PerformAlgorithmCallback() = default;
  // NOLINTNEXTLINE(google-explicit-constructor)
  PerformAlgorithmCallback(Proxy proxy) : proxy_(proxy) {}
  explicit PerformAlgorithmCallback(CkMigrateMessage* msg) : Callback(msg) {}
  using PUP::able::register_constructor;
  void invoke() override { proxy_.perform_algorithm(); }
  void pup(PUP::er& p) override { p | proxy_; }

  void register_with_charm() override {
    static bool done_registration{false};
    if (done_registration) {
      return;
    }
    done_registration = true;
    register_classes_with_charm<PerformAlgorithmCallback>();
  }

  std::optional<bool> is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr =
        dynamic_cast<const PerformAlgorithmCallback*>(&rhs);
    return downcast_ptr == nullptr ? std::nullopt : std::optional{true};
  }

  std::string name() const override {
    return "PerformAlgorithmCallback(" + pretty_type::name<Proxy>() + ")";
  }

 private:
  std::decay_t<Proxy> proxy_{};
};

/// \cond
template <typename Proxy>
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
PUP::able::PUP_ID PerformAlgorithmCallback<Proxy>::my_PUP_ID = 0;
template <typename SimpleAction, typename Proxy, typename... Args>
PUP::able::PUP_ID
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    SimpleActionCallback<SimpleAction, Proxy, Args...>::my_PUP_ID =
        0;  // NOLINT
template <typename SimpleAction, typename Proxy>
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
PUP::able::PUP_ID SimpleActionCallback<SimpleAction, Proxy>::my_PUP_ID =
    0;  // NOLINT
template <typename ThreadedAction, typename Proxy, typename... Args>
PUP::able::PUP_ID
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    ThreadedActionCallback<ThreadedAction, Proxy, Args...>::my_PUP_ID =
        0;  // NOLINT
template <typename ThreadedAction, typename Proxy>
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
PUP::able::PUP_ID ThreadedActionCallback<ThreadedAction, Proxy>::my_PUP_ID =
    0;  // NOLINT
/// \endcond

}  // namespace Parallel
