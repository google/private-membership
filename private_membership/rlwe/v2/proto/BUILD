load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//visibility:public"])

proto_library(
    name = "private_membership_rlwe_v2_proto",
    srcs = ["private_membership_rlwe_v2.proto"],
)

cc_proto_library(
    name = "private_membership_rlwe_v2_cc_proto",
    deps = [":private_membership_rlwe_v2_proto"],
)

proto_library(
    name = "rlwe_proto",
    srcs = ["rlwe.proto"],
)

cc_proto_library(
    name = "rlwe_cc_proto",
    deps = [":rlwe_proto"],
)

proto_library(
    name = "pir_proto",
    srcs = ["pir.proto"],
    deps = [
        ":private_membership_rlwe_v2_proto",
        ":rlwe_proto",
    ],
)

cc_proto_library(
    name = "pir_cc_proto",
    deps = [":pir_proto"],
)
