// Copyright 2026 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "research/oram_okvs/oram/lat/latoram.h"

#include <memory>
#include <utility>

#include "research/oram_okvs/oram/lat/latoram_client.h"
#include "research/oram_okvs/oram/lat/latoram_server.h"
#include "research/oram_okvs/oram/lat/latoram_storage.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

LatOram::LatOram(const LatOramConfig& config,
                 std::unique_ptr<LatOramStorage> storage) {
  if (!storage) {
    storage = std::make_unique<LatOramStorage>();
  }
  server_ = std::make_shared<LatOramServer>(std::move(storage));
  client_ = LatOramClient::Create(config, server_);
}

LatOram::LatOram(const LatOramConfig& config,
                 std::shared_ptr<LatOramServerInterface> server)
    : server_(std::move(server)),
      client_(LatOramClient::Create(config, server_)) {}

}  // namespace oram_okvs
