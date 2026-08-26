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
load("@rules_proto//proto:defs.bzl", "proto_library")
load("@rules_cc//cc:defs.bzl", "cc_library")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])

proto_library(
    name = "pir_proto",
    srcs = ["pir.proto"],
    deps = [
        "@com_github_google_shell//shell_encryption/rns:rns_serialization_proto",
    ],
)

cc_proto_library(
    name = "pir_cc_proto",
    deps = [":pir_proto"],
)

cc_library(
    name = "pir",
    hdrs = ["pir.h"],
    deps = [
        ":pir_cc_proto",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "fast_ntt",
    srcs = ["fast_ntt.cc"],
    hdrs = ["fast_ntt.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        "@com_github_google_highway//:hwy",
        "@com_github_google_shell//shell_encryption:dft_transformations",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_github_google_shell//shell_encryption/rns:rns_polynomial",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "fast_ntt_test",
    srcs = ["fast_ntt_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":fast_ntt",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption:dft_transformations",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption/rns:rns_context",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_github_google_shell//shell_encryption/rns:rns_polynomial",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
        "@com_google_absl//absl/status",
    ],
)

cc_library(
    name = "onion_pir_params",
    hdrs = ["onion_pir_params.h"],
    deps = [
        ":pir_cc_proto",
        "@com_github_google_shell//shell_encryption:constants",
        "@com_github_google_shell//shell_encryption:integral_types",
    ],
)

cc_library(
    name = "onion_pir_types",
    hdrs = ["onion_pir_types.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":fast_ntt",
        ":onion_pir_params",
        ":pir_cc_proto",
        "@com_google_absl//absl/strings:string_view",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_github_google_shell//shell_encryption/rns:rns_bfv_ciphertext",
        "@com_github_google_shell//shell_encryption/rns:rns_ciphertext",
        "@com_github_google_shell//shell_encryption/rns:rns_context",
        "@com_github_google_shell//shell_encryption/rns:rns_error_params",
        "@com_github_google_shell//shell_encryption/rns:rns_gadget",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_github_google_shell//shell_encryption/rns:rns_polynomial",
        "@com_github_google_shell//shell_encryption/rns:rns_serialization_cc_proto",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "onion_pir_client",
    srcs = ["onion_pir_client.cc"],
    hdrs = ["onion_pir_client.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":fast_ntt",
        ":onion_pir_params",
        ":onion_pir_types",
        ":pir",
        ":pir_cc_proto",
        "@com_github_google_shell//shell_encryption:constants",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_github_google_shell//shell_encryption/prng:hkdf_prng",
        "@com_github_google_shell//shell_encryption/rns:coefficient_encoder",
        "@com_github_google_shell//shell_encryption/rns:error_distribution",
        "@com_github_google_shell//shell_encryption/rns:rns_bfv_ciphertext",
        "@com_github_google_shell//shell_encryption/rns:rns_context",
        "@com_github_google_shell//shell_encryption/rns:rns_error_params",
        "@com_github_google_shell//shell_encryption/rns:rns_gadget",
        "@com_github_google_shell//shell_encryption/rns:rns_galois_key",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_github_google_shell//shell_encryption/rns:rns_polynomial",
        "@com_github_google_shell//shell_encryption/rns:rns_secret_key",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "onion_pir_server",
    srcs = ["onion_pir_server.cc"],
    hdrs = ["onion_pir_server.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":fast_ntt",
        ":onion_pir_params",
        ":onion_pir_types",
        ":pir",
        ":pir_cc_proto",
        "@com_github_google_shell//shell_encryption:constants",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_github_google_shell//shell_encryption/prng:hkdf_prng",
        "@com_github_google_shell//shell_encryption/rns:rns_bfv_ciphertext",
        "@com_github_google_shell//shell_encryption/rns:rns_context",
        "@com_github_google_shell//shell_encryption/rns:rns_error_params",
        "@com_github_google_shell//shell_encryption/rns:rns_gadget",
        "@com_github_google_shell//shell_encryption/rns:rns_galois_key",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_github_google_shell//shell_encryption/rns:rns_polynomial",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "onion_pir_test",
    srcs = ["onion_pir_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":fast_ntt",
        ":onion_pir_client",
        ":onion_pir_params",
        ":onion_pir_server",
        ":onion_pir_types",
        ":pir",
        ":pir_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption:constants",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_github_google_shell//shell_encryption/prng:hkdf_prng",
        "@com_github_google_shell//shell_encryption/rns:rns_bfv_ciphertext",
        "@com_github_google_shell//shell_encryption/rns:rns_context",
        "@com_github_google_shell//shell_encryption/rns:rns_error_params",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_github_google_shell//shell_encryption/rns:rns_polynomial",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_binary(
    name = "onion_pir_benchmark",
    srcs = ["onion_pir_benchmark.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":fast_ntt",
        ":onion_pir_client",
        ":onion_pir_params",
        ":onion_pir_server",
        ":onion_pir_types",
        "@com_github_google_benchmark//:benchmark",
        "@com_github_google_benchmark//:benchmark_main",
        "@com_github_google_shell//shell_encryption:dft_transformations",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption/rns:rns_context",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_google_absl//absl/strings",
    ],
)

cc_binary(
    name = "ntt_benchmark",
    srcs = ["ntt_benchmark.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":fast_ntt",
        "@com_github_google_benchmark//:benchmark",
        "@com_github_google_benchmark//:benchmark_main",
        "@com_github_google_shell//shell_encryption:dft_transformations",
        "@com_github_google_shell//shell_encryption:dft_transformations_hwy",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption/rns:rns_context",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_google_absl//absl/strings",
    ],
)
