// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/ControlErrors/Size/StateHistory.hpp"

#include <deque>
#include <memory>
#include <pup.h>
#include <pup_stl.h>
#include <unordered_map>
#include <utility>

#include "ControlSystem/ControlErrors/Size/AhSpeed.hpp"
#include "ControlSystem/ControlErrors/Size/DeltaR.hpp"
#include "ControlSystem/ControlErrors/Size/Info.hpp"
#include "ControlSystem/ControlErrors/Size/Initial.hpp"
#include "ControlSystem/ControlErrors/Size/State.hpp"
#include "DataStructures/DataVector.hpp"

namespace control_system::size {
StateHistory::StateHistory() { initialize_stored_control_errors(); }

StateHistory::StateHistory(const size_t num_times_to_store)
    : num_times_to_store_(num_times_to_store) {
  initialize_stored_control_errors();
}

StateHistory::StateHistory(const StateHistory& history) {
  num_times_to_store_ = history.num_times_to_store_;
  stored_control_errors_ = history.stored_control_errors_;
  // Don't need to copy the state
}

StateHistory& StateHistory::operator=(const StateHistory& history) {
  num_times_to_store_ = history.num_times_to_store_;
  stored_control_errors_ = history.stored_control_errors_;
  // Don't need to copy the state
  return *this;
}

void StateHistory::initialize_stored_control_errors() {
  const auto initialize_state = [this](auto state) {
    using StateType = std::decay_t<decltype(state)>;
    current_state_ = std::make_unique<StateType>();
    stored_control_errors_[current_state_->number()];
  };

  initialize_state(States::Initial{});
  initialize_state(States::DeltaR{});
  initialize_state(States::AhSpeed{});
}

void StateHistory::store(double time, const Info& info,
                         const ControlErrorArgs& control_error_args) {
  const auto store_state = [this, &time, &info,
                            &control_error_args](auto state) {
    using StateType = std::decay_t<decltype(state)>;
    current_state_ = std::make_unique<StateType>();
    double control_error =
        current_state_->control_error(info, control_error_args);
    size_t state_number = current_state_->number();
    stored_control_errors_.at(state_number)
        .emplace_back(std::make_pair(time, control_error));
    while (stored_control_errors_.at(state_number).size() >
           num_times_to_store_) {
      stored_control_errors_.at(state_number).pop_front();
    }
  };

  store_state(States::Initial{});
  store_state(States::DeltaR{});
  store_state(States::AhSpeed{});
}

std::deque<std::pair<double, double>> StateHistory::state_history(
    size_t state_number) const {
  return stored_control_errors_.at(state_number);
}

void StateHistory::pup(PUP::er& p) {
  p | num_times_to_store_;
  p | current_state_;
  p | stored_control_errors_;
}
}  // namespace control_system::size
