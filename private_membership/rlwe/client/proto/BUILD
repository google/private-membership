load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//visibility:public"])

proto_library(
    name = "private_membership_proto",
    srcs = ["private_membership.proto"],
)

cc_proto_library(
    name = "private_membership_cc_proto",
    deps = [":private_membership_proto"],
)

proto_library(
    name = "private_membership_rlwe_proto",
    srcs = ["private_membership_rlwe.proto"],
    deps = [
        ":private_membership_proto",
        "@com_github_google_shell//shell_encryption:serialization_proto",
    ],
)

cc_proto_library(
    name = "private_membership_rlwe_cc_proto",
    deps = [":private_membership_rlwe_proto"],
)
