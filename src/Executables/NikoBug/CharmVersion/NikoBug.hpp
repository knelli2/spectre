// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Declares the chares for NikoBug

#pragma once

#include <cstdint>
#include <pup.h>
#include "NikoBug.decl.h"

/// \cond
class CkArgMsg;
class CkCallback;
class Sender;
class Receiver;
class CProxy_Receiver;
class CProxy_Sender;
/// \endcond

class NikoBug : public CBase_NikoBug {
 public:
  explicit NikoBug(CkArgMsg* msg);
};

class Sender : public CBase_Sender {
 public:
  Sender(CProxy_Receiver receiver_proxy);

  void send_messages_to_receiver();

  void pup(PUP::er& p) override;

 private:
  CProxy_Receiver receiver_proxy_;
  uint64_t messages_sent_{0};
};

class Receiver : public CBase_Receiver {
 public:
  Receiver(CProxy_Sender sender_proxy);

  void receive_messages_from_sender();
  [[noreturn]] void print_results() const;

  void pup(PUP::er& p) override;

 private:
  CProxy_Sender sender_proxy_;
  uint64_t messages_received_{0};
};
/// \endcond
