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
    name = "gf128",
    srcs = ["gf128.cc"],
    hdrs = ["gf128.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
    ],
)

cc_test(
    name = "gf128_test",
    srcs = ["gf128_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/numeric:int128",
    ],
)

cc_binary(
    name = "gf128_benchmark",
    srcs = ["gf128_benchmark.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        "@com_github_google_benchmark//:benchmark",
        "@com_github_google_benchmark//:benchmark_main",
        "@com_google_absl//absl/base:core_headers",
    ],
)

cc_library(
    name = "okvs",
    hdrs = ["okvs.h"],
    deps = [
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "random_band_okvs",
    hdrs = ["random_band_okvs.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        ":okvs",
        "@boringssl//:ssl",
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "random_band_okvs_test",
    srcs = ["random_band_okvs_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        ":random_band_okvs",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)

cc_binary(
    name = "random_band_okvs_benchmark",
    srcs = ["random_band_okvs_benchmark.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        ":random_band_okvs",
        "@com_github_google_benchmark//:benchmark",
        "@com_github_google_benchmark//:benchmark_main",
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "polynomial_okvs",
    hdrs = ["polynomial_okvs.h"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        ":okvs",
        "@boringssl//:ssl",
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "polynomial_okvs_test",
    srcs = ["polynomial_okvs_test.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        ":polynomial_okvs",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
    ],
)

cc_binary(
    name = "polynomial_okvs_benchmark",
    srcs = ["polynomial_okvs_benchmark.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        ":polynomial_okvs",
        "@com_github_google_benchmark//:benchmark",
        "@com_github_google_benchmark//:benchmark_main",
        "@com_google_absl//absl/strings",
    ],
)

cc_binary(
    name = "rb_okvs_benchmark",
    srcs = ["rb_okvs_benchmark.cc"],
    copts = [
        "-msse4.2",
        "-mpclmul",
        "-mavx",
        "-mavx2",
    ],
    deps = [
        ":gf128",
        ":random_band_okvs",
        "@com_google_absl//absl/flags:flag",
        "@com_google_absl//absl/flags:parse",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/time",
    ],
)
