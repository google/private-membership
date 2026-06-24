load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "pir_client",
    srcs = ["pir_client.cc"],
    hdrs = ["pir_client.h"],
    deps = [
        "//private_membership/rlwe/v2/crypto:encryption",
        "//private_membership/rlwe/v2/crypto:polynomial",
        "//private_membership/rlwe/v2/crypto/pir:pir_params",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)

cc_test(
    name = "pir_client_test",
    srcs = ["pir_client_test.cc"],
    deps = [
        ":pir_client",
        "//private_membership/rlwe/v2/crypto:encryption",
        "//private_membership/rlwe/v2/crypto/pir:pir_params",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_github_google_shell//shell_encryption/prng:chacha_prng",
        "@com_github_google_shell//shell_encryption/testing:status_is_fork",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)
