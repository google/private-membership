load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//:__subpackages__"])

proto_library(
    name = "regression_test_data_proto",
    srcs = ["regression_test_data.proto"],
    deps = ["//private_membership/rlwe/client/proto:private_membership_rlwe_proto"],
)

cc_proto_library(
    name = "regression_test_data_cc_proto",
    deps = [":regression_test_data_proto"],
)

filegroup(
    name = "regression_test_data",
    testonly = 1,
    srcs = [
        "test_data.binarypb",
    ],
)
