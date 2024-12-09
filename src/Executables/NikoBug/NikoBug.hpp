// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Declares the chares for NikoBug

#pragma once

#include <cstdint>
#include <pup.h>
#include "Executables/NikoBug/NikoBug.decl.h"
#include "Parallel/Reduction.hpp"

/// \cond
class CkArgMsg;
class CkCallback;
class Sender;
class Receiver;
class CProxy_Receiver;
/// \endcond

class NikoBug : public CBase_NikoBug {
 public:
  explicit NikoBug(CkArgMsg* msg);
};

class Sender : public CBase_Sender {
 public:
  Sender(CProxy_Receiver receiver_proxy);

  void send_messages_to_receiver();

 private:
  CProxy_Receiver receiver_proxy_;
};

class Receiver : public CBase_Receiver {
 public:
  Receiver();

  void receive_messages_from_sender();
  [[noreturn]] void print_results() const;

  void pup(PUP::er& p) override;

 private:
  uint64_t messages_received;
};
/// \endcond
