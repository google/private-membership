package(default_visibility = ["//:__subpackages__"])

cc_library(
    name = "private_membership_rlwe_client",
    srcs = ["private_membership_rlwe_client.cc"],
    hdrs = [
        "private_membership_rlwe_client.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//private_membership/rlwe/client/cpp/internal:constants",
        "//private_membership/rlwe/client/cpp/internal:crypto_utils",
        "//private_membership/rlwe/client/cpp/internal:encrypted_bucket_id",
        "//private_membership/rlwe/client/cpp/internal:hashed_bucket_id",
        "//private_membership/rlwe/client/cpp/internal:rlwe_id_utils",
        "//private_membership/rlwe/client/cpp/internal:rlwe_params",
        "//private_membership/rlwe/client/cpp/internal:utils",
        "//private_membership/rlwe/client/proto:private_membership_cc_proto",
        "//private_membership/rlwe/client/proto:private_membership_rlwe_cc_proto",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:polynomial",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_github_google_shell//shell_encryption:symmetric_encryption",
        "@com_github_google_shell//shell_encryption:symmetric_encryption_with_prng",
        "@com_github_google_shell//shell_encryption:transcription",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
    ],
)

cc_test(
    name = "private_membership_rlwe_client_regression_test",
    srcs = ["private_membership_rlwe_client_regression_test.cc"],
    data = ["//private_membership/rlwe/client/cpp/internal/testing/regression_test_data"],
    deps = [
        ":private_membership_rlwe_client",
        "//private_membership/rlwe/client/cpp/internal/testing/regression_test_data:regression_test_data_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
        "@com_google_absl//absl/status",
    ],
)
