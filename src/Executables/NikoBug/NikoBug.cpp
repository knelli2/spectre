// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Executables/NikoBug/NikoBug.hpp"

#include <Executables/NikoBug/NikoBug.decl.h>
#include <charm++.h>
#include <ckcallback.h>
#include <cstdint>
#include <cstdio>
#include <pup.h>

#include "Utilities/System/Exit.hpp"
#include "Utilities/System/ParallelInfo.hpp"

namespace {
constexpr uint64_t expected_number_of_messages = 4294967296;  // 2^32
}  // namespace

NikoBug::NikoBug(CkArgMsg* /*msg*/) {
  const int num_procs = sys::number_of_procs();
  printf("Running on %u cores\n", num_procs);  // NOLINT
  CProxy_Receiver receiver_proxy = CProxy_Receiver::ckNew();
  receiver_proxy[0].insert(0_st);
  receiver_proxy.doneInserting();

  CProxy_Sender sender_proxy = CProxy_Sender::ckNew();
  sender_proxy[0].insert(receiver_proxy, num_procs > 1 ? 1_st : 0_st);
  sender_proxy.doneInserting();

  sender_proxy[0].send_messages_to_receiver();
}

Sender::Sender(CProxy_Receiver receiver_proxy)
    : receiver_proxy_(std::move(receiver_proxy)) {}

void Sender::send_messages_to_receiver() {
  for (uint64_t i = 0; i < expected_number_of_messages; i++) {
    receiver_proxy_.receive_messages_from_sender();
    if (i % 100000000 == 0) {
      printf("Sent %lu messages already\n", i);  // NOLINT
    }
  }

  CkStartQD(CkCallback(CkIndex_Receiver::print_results(), receiver_proxy_));
}

Receiver::Receiver() : messages_received(0) {}

void Receiver::receive_messages_from_sender() { ++messages_received; }

void Receiver::print_results() const {
  // NOLINTNEXTLINE
  printf(
      "Expected number of messages: %lu\n"
      "Received number of messages: %lu\n",
      expected_number_of_messages, messages_received);

  sys::exit();
}

void Receiver::pup(PUP::er& p) { p | messages_received; }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "Executables/NikoBug/NikoBug.def.h"
#pragma GCC diagnostic pop
