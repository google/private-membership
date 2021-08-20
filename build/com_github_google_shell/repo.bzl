# Copyright 2021 The Cross-Media Measurement Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
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
            sha256 = "a6404a1fdc794c8e5424ed7233a925872d2abd7b08462a366f4a6485ea0747e7",
            strip_prefix = "shell-encryption-4e2c0dded993ec96eb319453e35e7d95ce689b45",
            url = "https://github.com/google/shell-encryption/archive/4e2c0dded993ec96eb319453e35e7d95ce689b45.tar.gz",
        )
