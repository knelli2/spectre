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
  const int global_num_procs = sys::number_of_procs();
  const int num_nodes = sys::number_of_nodes();
  // NOLINTNEXTLINE
  printf("Running on %u nodes and %u cores\n", num_nodes, global_num_procs);
  CProxy_Receiver receiver_proxy = CProxy_Receiver::ckNew();
  const int receiver_core = 0;
  printf("Allocating Receiver on global proc %u\n", receiver_core);  // NOLINT
  receiver_proxy[0].insert(receiver_core);
  receiver_proxy.doneInserting();

  CProxy_Sender sender_proxy = CProxy_Sender::ckNew();
  int sender_core = 0;
  if (num_nodes > 1) {
    sender_core = sys::first_proc_on_node(1);
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
  printf("Sender: My node = %u, my core = %u\n", sys::my_node(),
         sys::my_proc());
}

void Sender::send_messages_to_receiver() {
  for (uint64_t i = 0; i < expected_number_of_messages; i++) {
    receiver_proxy_.receive_messages_from_sender();
    if (i % 100000000 == 0) {
      printf("Sent %lu messages already\n", i);  // NOLINT
    }
  }
}

Receiver::Receiver() {
  // NOLINTNEXTLINE
  printf("Receiver: My node = %u, my core = %u\n", sys::my_node(),
         sys::my_proc());
}

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
