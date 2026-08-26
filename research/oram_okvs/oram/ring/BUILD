# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_cc//cc:defs.bzl", "cc_test")

package(
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])

cc_library(
    name = "ringoram_storage",
    srcs = ["ringoram_storage.cc"],
    hdrs = ["ringoram_storage.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        "//research/oram_okvs/oram:crypto_utils",
        "//research/oram_okvs/oram:oblivious_container",
        "//research/oram_okvs/oram:oram_cc_proto",
        "//research/oram_okvs/oram:types",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings:str_format",
    ],
)

cc_library(
    name = "ringoram_server",
    srcs = ["ringoram_server.cc"],
    hdrs = ["ringoram_server.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":ringoram_storage",
        "//research/oram_okvs/oram:block_retrieval",
        "//research/oram_okvs/oram:oram_cc_proto",
        "//research/oram_okvs/pir:onion_pir_params",
        "//research/oram_okvs/pir:onion_pir_server",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)

cc_library(
    name = "ringoram_client",
    srcs = [
        "ok_ringoram_client.cc",
        "okp_ringoram_client.cc",
        "ringoram_client.cc",
        "tree_ringoram_client.cc",
    ],
    hdrs = [
        "ok_ringoram_client.h",
        "okp_ringoram_client.h",
        "ringoram_client.h",
        "tree_ringoram_client.h",
    ],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":ringoram_server",
        ":ringoram_storage",
        "//research/oram_okvs/okvs:polynomial_okvs",
        "//research/oram_okvs/okvs:random_band_okvs",
        "//research/oram_okvs/oram",
        "//research/oram_okvs/oram:block_retrieval",
        "//research/oram_okvs/oram:crypto_utils",
        "//research/oram_okvs/oram:oram_cc_proto",
        "//research/oram_okvs/oram:types",
        "//research/oram_okvs/pir:onion_pir_params",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "ringoram",
    srcs = ["ringoram.cc"],
    hdrs = ["ringoram.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":ringoram_client",
        ":ringoram_server",
        ":ringoram_storage",
        "//research/oram_okvs/oram",
        "//research/oram_okvs/oram:types",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "ringoram_test",
    srcs = ["ringoram_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":ringoram",
        ":ringoram_client",
        ":ringoram_server",
        ":ringoram_storage",
        "//research/oram_okvs/oram",
        "//research/oram_okvs/oram:types",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)
