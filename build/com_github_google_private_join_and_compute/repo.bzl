# Copyright 2021 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Repository rules/macros for com_github_google_private_join_and_compute.
"""

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def com_github_google_private_join_and_compute_repo():
    if "com_github_google_private_join_and_compute" not in native.existing_rules():
        http_archive(
            name = "com_github_google_private_join_and_compute",
            sha256 = "b1a83e1bc778fe902b782ae6d06fdf590a1f74684954c05592463ddad75f8ddb",
            strip_prefix = "private-join-and-compute-505ba981d66c9e5e73e18cfa647b4685f74784cb",
            url = "https://github.com/google/private-join-and-compute/archive/505ba981d66c9e5e73e18cfa647b4685f74784cb.tar.gz",
        )
