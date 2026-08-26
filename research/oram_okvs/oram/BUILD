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
load("@rules_proto//proto:defs.bzl", "proto_library")
load("@rules_cc//cc:defs.bzl", "cc_library")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])

proto_library(
    name = "oram_proto",
    srcs = ["oram.proto"],
    deps = [
        "//research/oram_okvs/pir:pir_proto",
    ],
)

cc_proto_library(
    name = "oram_cc_proto",
    deps = [":oram_proto"],
)

cc_library(
    name = "types",
    hdrs = ["types.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        "//research/oram_okvs/pir:onion_pir_params",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "oram",
    hdrs = ["oram.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":types",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "crypto_utils",
    srcs = ["crypto_utils.cc"],
    hdrs = ["crypto_utils.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":types",
        "@boringssl//:ssl",
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "block_retrieval",
    srcs = ["block_retrieval.cc"],
    hdrs = ["block_retrieval.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":oram_cc_proto",
        ":types",
        "//research/oram_okvs/pir:onion_pir_client",
        "//research/oram_okvs/pir:onion_pir_params",
        "//research/oram_okvs/pir:onion_pir_server",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "block_retrieval_test",
    srcs = ["block_retrieval_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":block_retrieval",
        ":oram_cc_proto",
        "//research/oram_okvs/pir:onion_pir_params",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)

cc_library(
    name = "oblivious_container",
    hdrs = ["oblivious_container.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":crypto_utils",
        ":oram_cc_proto",
        ":types",
        "//research/oram_okvs/okvs:polynomial_okvs",
        "//research/oram_okvs/okvs:random_band_okvs",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "oblivious_container_test",
    srcs = ["oblivious_container_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":crypto_utils",
        ":oblivious_container",
        ":types",
        "//research/oram_okvs/okvs:polynomial_okvs",
        "//research/oram_okvs/okvs:random_band_okvs",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)
