licenses(["notice"])

package(default_visibility = ["//:__subpackages__"])

cc_library(
    name = "client",
    srcs = ["client.cc"],
    hdrs = ["client.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":client_helper",
        "//private_membership/rlwe/batch/cpp:constants",
        "//private_membership/rlwe/batch/cpp:context",
        "//private_membership/rlwe/batch/cpp:padding",
        "//private_membership/rlwe/batch/cpp:parameters",
        "//private_membership/rlwe/batch/cpp:prng",
        "//private_membership/rlwe/batch/proto:client_cc_proto",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_shell//shell_encryption:context",
        "@com_github_google_shell//shell_encryption:galois_key",
        "@com_github_google_shell//shell_encryption:integral_types",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:oblivious_expand",
        "@com_github_google_shell//shell_encryption:polynomial",
        "@com_github_google_shell//shell_encryption:serialization_cc_proto",
        "@com_github_google_shell//shell_encryption:symmetric_encryption",
        "@com_github_google_shell//shell_encryption:symmetric_encryption_with_prng",
        "@com_github_google_shell//shell_encryption:transcription",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

# This library contains functions that are shared with tests.
cc_library(
    name = "client_helper",
    srcs = ["client_helper.cc"],
    hdrs = ["client_helper.h"],
    deps = [
        "//private_membership/rlwe/batch/cpp:parameters",
        "//private_membership/rlwe/batch/proto:client_cc_proto",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
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
        ":client_helper",
        "//private_membership/rlwe/batch/cpp:constants",
        "//private_membership/rlwe/batch/cpp:context",
        "//private_membership/rlwe/batch/cpp:encoding",
        "//private_membership/rlwe/batch/cpp:padding",
        "//private_membership/rlwe/batch/cpp:prng",
        "//private_membership/rlwe/batch/cpp:test_helper",
        "//private_membership/rlwe/batch/proto:client_cc_proto",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption:context",
        "@com_github_google_shell//shell_encryption:galois_key",
        "@com_github_google_shell//shell_encryption:int256",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:oblivious_expand",
        "@com_github_google_shell//shell_encryption:polynomial",
        "@com_github_google_shell//shell_encryption:serialization_cc_proto",
        "@com_github_google_shell//shell_encryption:symmetric_encryption",
        "@com_github_google_shell//shell_encryption:symmetric_encryption_with_prng",
        "@com_github_google_shell//shell_encryption:transcription",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)
