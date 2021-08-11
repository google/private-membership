load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//:__subpackages__"])

licenses(["notice"])

proto_library(
    name = "shared_proto",
    srcs = ["shared.proto"],
    visibility = ["//visibility:public"],
    deps = [
        "@com_github_google_shell//:serialization_proto",
    ],
)

cc_proto_library(
    name = "shared_cc_proto",
    visibility = ["//visibility:public"],
    deps = [
        ":shared_proto",
    ],
)

cc_library(
    name = "shared",
    srcs = ["shared.cc"],
    hdrs = ["shared.h"],
    deps = [
        ":shared_cc_proto",
        "@com_github_google_shell//:context",
        "@com_github_google_shell//:integral_types",
        "@com_github_google_shell//:montgomery",
        "@com_github_google_shell//prng",
        "@com_github_google_shell//prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "shared_test",
    srcs = ["shared_test.cc"],
    deps = [
        ":shared",
        "@com_github_google_googletest//:gtest_main",
    ],
)
