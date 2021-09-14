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
Repository rules/macros for com_github_gflags_gflags.
"""

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def com_github_google_shell_repo():
    if "com_github_google_shell" not in native.existing_rules():
        http_archive(
            name = "com_github_google_shell",
            sha256 = "807d42caed3056cea63b9048a2fd69122c071740d43c9de546cc0fabded87a5c",
            strip_prefix = "shell-encryption-507781e129a03f8178c9716d79163fae23d34b6a",
            url = "https://github.com/google/shell-encryption/archive/507781e129a03f8178c9716d79163fae23d34b6a.tar.gz",
        )
