load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

licenses(["notice"])

package(default_visibility = ["//:__subpackages__"])

proto_library(
    name = "client_proto",
    srcs = ["client.proto"],
    visibility = ["//visibility:public"],
    deps = [
        "//private_membership/rlwe/batch/cpp:shared_proto",
        "@com_github_google_shell//:serialization_proto",
    ],
)

cc_proto_library(
    name = "client_cc_proto",
    visibility = ["//visibility:public"],
    deps = [
        ":client_proto",
    ],
)

cc_library(
    name = "client",
    srcs = ["client.cc"],
    hdrs = ["client.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":client_cc_proto",
        ":client_helper",
        "//private_membership/rlwe/batch/cpp:shared",
        "//private_membership/rlwe/batch/cpp:shared_cc_proto",
        "@com_github_google_shell//:context",
        "@com_github_google_shell//:galois_key",
        "@com_github_google_shell//:integral_types",
        "@com_github_google_shell//:oblivious_expand",
        "@com_github_google_shell//:polynomial",
        "@com_github_google_shell//:serialization_cc_proto",
        "@com_github_google_shell//:symmetric_encryption",
        "@com_github_google_shell//:symmetric_encryption_with_prng",
        "@com_github_google_shell//:transcription",
        "@com_github_google_shell//prng",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)

# This library contains functions that are shared with tests.
cc_library(
    name = "client_helper",
    srcs = ["client_helper.cc"],
    hdrs = ["client_helper.h"],
    deps = [
        ":client_cc_proto",
        "//private_membership/rlwe/batch/cpp:shared_cc_proto",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "client_test",
    srcs = ["client_test.cc"],
    deps = [
        ":client",
        ":client_cc_proto",
        ":client_helper",
        "//private_membership/rlwe/batch/cpp:shared",
        "//private_membership/rlwe/batch/cpp:shared_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//:galois_key",
        "@com_github_google_shell//:oblivious_expand",
        "@com_github_google_shell//:polynomial",
        "@com_github_google_shell//:serialization_cc_proto",
        "@com_github_google_shell//:symmetric_encryption",
        "@com_github_google_shell//:symmetric_encryption_with_prng",
        "@com_github_google_shell//:transcription",
        "@com_google_absl//absl/strings",
    ],
)
