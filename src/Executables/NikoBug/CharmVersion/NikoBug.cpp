// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NikoBug.hpp"

#include <charm++.h>
#include <charm.h>
#include <ckcallback.h>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <iomanip>
#include <pup.h>
#include <sstream>
#include <string>

namespace {
constexpr uint64_t expected_number_of_messages = 4294967296;  // 2^32
constexpr uint64_t batch_size = 100000000;

std::string pretty_wall_time() {
  // Subseconds don't really matter so just ignore them. This gives nice round
  // numbers.
  int total = static_cast<int>(CkWallTimer());
  const int day = total / (24 * 3600);

  total %= (24 * 3600);
  const int hour = total / 3600;

  total %= 3600;
  const int minutes = total / 60;

  total %= 60;
  const int seconds = total;

  std::stringstream ss{};
  ss << std::setfill('0');
  if (day > 0) {
    ss << std::setw(2) << day << "-";
  }

  // std::setw() isn't sticky so it has to be used for every insertion
  ss << std::setw(2) << hour << ":";
  ss << std::setw(2) << minutes << ":";
  ss << std::setw(2) << seconds;
  return ss.str();
}
}  // namespace

NikoBug::NikoBug(CkArgMsg* /*msg*/) {
  const int global_num_procs = CkNumPes();
  const int num_nodes = CkNumNodes();
  // NOLINTNEXTLINE
  printf("Running on %u nodes and %u cores\n", num_nodes, global_num_procs);
  CProxy_Receiver receiver_proxy = CProxy_Receiver::ckNew();
  CProxy_Sender sender_proxy = CProxy_Sender::ckNew();

  const int receiver_core = 0;
  printf("Allocating Receiver on global proc %u\n", receiver_core);  // NOLINT
  receiver_proxy[0].insert(sender_proxy, receiver_core);
  receiver_proxy.doneInserting();

  int sender_core = 0;
  if (num_nodes > 1) {
    sender_core = CkNodeFirst(1);
  } else if (global_num_procs > 1) {
    sender_core = 1;
  }
  printf("Allocating Sender on global proc %u\n", sender_core);  // NOLINT
  sender_proxy[0].insert(receiver_proxy, sender_core);
  sender_proxy.doneInserting();

  sender_proxy[0].send_messages_to_receiver();

  CkStartQD(CkCallback(CkIndex_Receiver::print_results(), receiver_proxy));
}

Sender::Sender(CProxy_Receiver receiver_proxy)
    : receiver_proxy_(std::move(receiver_proxy)) {
  // NOLINTNEXTLINE
  printf("Sender: My node = %u, my core = %u\n", CkMyNode(), CkMyPe());
}

void Sender::send_messages_to_receiver() {
  for (uint64_t i = 0; i < batch_size; i++) {
    if (messages_sent_ >= expected_number_of_messages) {
      break;
    }

    receiver_proxy_.receive_messages_from_sender();
    ++messages_sent_;
  }

  // NOLINTNEXTLINE
  printf("Sent %lu messages at %s\n", messages_sent_,
         pretty_wall_time().c_str());
}

void Sender::pup(PUP::er& p) { p | messages_sent_; }

Receiver::Receiver(CProxy_Sender sender_proxy)
    : sender_proxy_(std::move(sender_proxy)) {
  // NOLINTNEXTLINE
  printf("Receiver: My node = %u, my core = %u\n", CkMyNode(), CkMyPe());
}

void Receiver::receive_messages_from_sender() {
  ++messages_received_;
  if (messages_received_ % batch_size == 0) {
    // NOLINTNEXTLINE
    printf("Received %lu messages at %s\n", messages_received_,
           pretty_wall_time().c_str());
    sender_proxy_.send_messages_to_receiver();
  }
}

void Receiver::print_results() const {
  // NOLINTNEXTLINE
  printf(
      "At: %s\n"
      " Expected number of messages: %lu\n"
      " Received number of messages: %lu\n",
      pretty_wall_time().c_str(), expected_number_of_messages,
      messages_received_);

  CkExit(0);
  // the following call is never reached, but suppresses the warning that
  // a 'noreturn' function does return
  std::terminate();
}

void Receiver::pup(PUP::er& p) { p | messages_received_; }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "NikoBug.def.h"
#pragma GCC diagnostic pop
