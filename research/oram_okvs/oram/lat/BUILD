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

load("@rules_cc//cc:defs.bzl", "cc_binary")
load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_cc//cc:defs.bzl", "cc_test")

package(
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])

cc_library(
    name = "multi_perm",
    hdrs = ["multi_perm.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "latoram_storage",
    srcs = ["latoram_storage.cc"],
    hdrs = ["latoram_storage.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        "//research/oram_okvs/okvs:polynomial_okvs",
        "//research/oram_okvs/okvs:random_band_okvs",
        "//research/oram_okvs/oram:oblivious_container",
        "//research/oram_okvs/oram:types",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "latoram_server",
    srcs = ["latoram_server.cc"],
    hdrs = ["latoram_server.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":latoram_storage",
        "//research/oram_okvs/oram:block_retrieval",
        "//research/oram_okvs/oram:crypto_utils",
        "//research/oram_okvs/oram:oram_cc_proto",
        "//research/oram_okvs/oram:types",
        "//research/oram_okvs/pir:onion_pir_params",
        "//research/oram_okvs/pir:onion_pir_server",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "latoram_client",
    srcs = [
        "latoram_client.cc",
        "map_latoram_client.cc",
        "ok_latoram_client.cc",
        "okp_latoram_client.cc",
    ],
    hdrs = [
        "latoram_client.h",
        "map_latoram_client.h",
        "ok_latoram_client.h",
        "okp_latoram_client.h",
    ],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":latoram_server",
        ":latoram_storage",
        ":multi_perm",
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
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "latoram",
    srcs = ["latoram.cc"],
    hdrs = ["latoram.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":latoram_client",
        ":latoram_server",
        ":latoram_storage",
        "//research/oram_okvs/oram",
        "//research/oram_okvs/oram:types",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "latoram_test",
    srcs = ["latoram_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":latoram",
        ":latoram_client",
        ":latoram_server",
        ":latoram_storage",
        ":multi_perm",
        "//research/oram_okvs/oram",
        "//research/oram_okvs/oram:crypto_utils",
        "//research/oram_okvs/oram:oram_cc_proto",
        "//research/oram_okvs/oram:types",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_binary(
    name = "latoram_benchmark",
    srcs = ["latoram_benchmark.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":latoram",
        ":multi_perm",
        "//research/oram_okvs/oram:types",
        "@com_github_google_benchmark//:benchmark",
        "@com_github_google_benchmark//:benchmark_main",
    ],
)
