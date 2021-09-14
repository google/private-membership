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
Adds external repos necessary for private-membership.
"""

load("//build/tink_base:repo.bzl", "tink_base_repo")
load("//build/tink_cc:repo.bzl", "tink_cc_repo")
load("//build/rules_cc:repo.bzl", "rules_cc_repo")
load("//build/rules_proto:repo.bzl", "rules_proto_repo")
load("//build/com_github_google_googletest:repo.bzl", "com_github_google_googletest_repo")
load("//build/com_google_absl:repo.bzl", "com_google_absl_repo")
load("//build/com_github_google_glog:repo.bzl", "com_github_google_glog_repo")
load("//build/com_github_gflags_gflags:repo.bzl", "com_github_gflags_gflags_repo")
load("//build/com_github_google_shell:repo.bzl", "com_github_google_shell_repo")
load("//build/boringssl:repo.bzl", "boringssl_repo")

def private_membership_repositories():
    """
    Adds all external repos necessary for private_membership.
    """
    tink_base_repo()
    tink_cc_repo()
    rules_cc_repo()
    rules_proto_repo()
    com_github_google_googletest_repo()
    com_google_absl_repo()
    com_github_google_glog_repo()
    com_github_gflags_gflags_repo()
    com_github_google_shell_repo()
    boringssl_repo()
