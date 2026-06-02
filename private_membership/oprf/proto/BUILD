load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//visibility:public"])

proto_library(
    name = "oprf_proto",
    srcs = ["oprf.proto"],
    option_deps = [
        "@com_google_protobuf//:cpp_features_proto",
    ],
    deps = [
        "//private_membership/rlwe/client/proto:private_membership_proto",
        "//private_membership/rlwe/client/proto:private_membership_rlwe_proto",
    ],
)

cc_proto_library(
    name = "oprf_cc_proto",
    deps = [":oprf_proto"],
)
