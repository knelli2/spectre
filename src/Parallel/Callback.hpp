// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines Parallel::Callback.

#pragma once

#include <pup.h>
#include <tuple>
#include <utility>

#include "Parallel/Invoke.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace Parallel {
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
  virtual bool is_equal_to(const Callback& rhs) const = 0;
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

  bool is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr =
        dynamic_cast<const SimpleActionCallback*>(&rhs);
    if (downcast_ptr == nullptr) {
      return false;
    }
    return // proxy_ == downcast_ptr->proxy_ and
        args_ == downcast_ptr->args_;
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

  bool is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr =
        dynamic_cast<const SimpleActionCallback*>(&rhs);
    return downcast_ptr != nullptr;
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

  bool is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr =
        dynamic_cast<const ThreadedActionCallback*>(&rhs);
    if (downcast_ptr == nullptr) {
      return false;
    }
    return  // proxy_ == downcast_ptr->proxy_ and
        args_ == downcast_ptr->args_;
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

  bool is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr =
        dynamic_cast<const ThreadedActionCallback*>(&rhs);
    return downcast_ptr != nullptr;
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

  bool is_equal_to(const Callback& rhs) const override {
    const auto* downcast_ptr =
        dynamic_cast<const PerformAlgorithmCallback*>(&rhs);
    return downcast_ptr != nullptr;
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
