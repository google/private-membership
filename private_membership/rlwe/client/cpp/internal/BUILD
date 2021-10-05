package(default_visibility = ["//:__subpackages__"])

cc_library(
    name = "aes_ctr_256_with_fixed_iv",
    srcs = [
        "aes_ctr_256_with_fixed_iv.cc",
    ],
    hdrs = [
        "aes_ctr_256_with_fixed_iv.h",
    ],
    deps = [
        "@boringssl//:ssl",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "aes_ctr_256_with_fixed_iv_test",
    srcs = [
        "aes_ctr_256_with_fixed_iv_test.cc",
    ],
    deps = [
        ":aes_ctr_256_with_fixed_iv",
        "@boringssl//:ssl",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/testing:status_is_fork",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
        "@com_google_absl//absl/memory",
    ],
)

cc_library(
    name = "crypto_utils",
    srcs = [
        "crypto_utils.cc",
    ],
    hdrs = [
        "crypto_utils.h",
    ],
    deps = [
        ":aes_ctr_256_with_fixed_iv",
        ":id_utils",
        "//private_membership/rlwe/client/proto:private_membership_cc_proto",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:bn_util",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "crypto_utils_test",
    srcs = ["crypto_utils_test.cc"],
    deps = [
        ":aes_ctr_256_with_fixed_iv",
        ":crypto_utils",
        "//private_membership/rlwe/client/proto:private_membership_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/testing:status_is_fork",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
    ],
)

cc_library(
    name = "id_utils",
    srcs = [
        "id_utils.cc",
    ],
    hdrs = [
        "id_utils.h",
    ],
    deps = [
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "id_utils_test",
    srcs = [
        "id_utils_test.cc",
    ],
    deps = [
        ":id_utils",
        "@com_github_google_googletest//:gtest_main",
    ],
)

cc_library(
    name = "oprf_utils",
    srcs = [
        "oprf_utils.cc",
    ],
    hdrs = [
        "oprf_utils.h",
    ],
    deps = [
        "//private_membership/rlwe/client/proto:private_membership_cc_proto",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_shell//shell_encryption:statusor_fork",
    ],
)

cc_test(
    name = "oprf_utils_test",
    srcs = [
        "oprf_utils_test.cc",
    ],
    deps = [
        ":oprf_utils",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_shell//shell_encryption/testing:status_is_fork",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
    ],
)

cc_library(
    name = "constants",
    hdrs = [
        "constants.h",
    ],
    deps = [
        "@boringssl//:ssl",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption/prng:chacha_prng",
        "@com_github_google_shell//shell_encryption/prng:single_thread_chacha_prng",
    ],
)

cc_library(
    name = "encrypted_bucket_id",
    srcs = [
        "encrypted_bucket_id.cc",
    ],
    hdrs = [
        "encrypted_bucket_id.h",
    ],
    deps = [
        ":rlwe_id_utils",
        ":utils",
        "//private_membership/rlwe/client/cpp/internal:crypto_utils",
        "//private_membership/rlwe/client/proto:private_membership_rlwe_cc_proto",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "encrypted_bucket_id_test",
    srcs = [
        "encrypted_bucket_id_test.cc",
    ],
    deps = [
        ":encrypted_bucket_id",
        ":rlwe_id_utils",
        "//private_membership/rlwe/client/cpp/internal/testing:constants",
        "//private_membership/rlwe/client/proto:private_membership_rlwe_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/testing:status_is_fork",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
    ],
)

cc_library(
    name = "hashed_bucket_id",
    srcs = [
        "hashed_bucket_id.cc",
    ],
    hdrs = [
        "hashed_bucket_id.h",
    ],
    deps = [
        ":rlwe_id_utils",
        ":utils",
        "//private_membership/rlwe/client/proto:private_membership_cc_proto",
        "//private_membership/rlwe/client/proto:private_membership_rlwe_cc_proto",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "hashed_bucket_id_test",
    srcs = [
        "hashed_bucket_id_test.cc",
    ],
    deps = [
        ":hashed_bucket_id",
        "//private_membership/rlwe/client/proto:private_membership_cc_proto",
        "//private_membership/rlwe/client/proto:private_membership_rlwe_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/testing:status_is_fork",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
    ],
)

cc_library(
    name = "rlwe_id_utils",
    srcs = [
        "rlwe_id_utils.cc",
    ],
    hdrs = [
        "rlwe_id_utils.h",
    ],
    deps = [
        ":constants",
        "//private_membership/rlwe/client/cpp/internal:crypto_utils",
        "//private_membership/rlwe/client/proto:private_membership_cc_proto",
        "//private_membership/rlwe/client/proto:private_membership_rlwe_cc_proto",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "rlwe_id_utils_test",
    srcs = ["rlwe_id_utils_test.cc"],
    deps = [
        ":constants",
        ":rlwe_id_utils",
        "//private_membership/rlwe/client/cpp/internal/testing:constants",
        "//private_membership/rlwe/client/proto:private_membership_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/testing:status_is_fork",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
        "@com_google_absl//absl/status",
    ],
)

cc_library(
    name = "rlwe_params",
    srcs = [
        "rlwe_params.cc",
    ],
    hdrs = [
        "rlwe_params.h",
    ],
    deps = [
        ":constants",
        "//private_membership/rlwe/client/proto:private_membership_rlwe_cc_proto",
        "@com_github_google_shell//shell_encryption:context",
        "@com_github_google_shell//shell_encryption:error_params",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption:statusor_fork",
    ],
)

cc_library(
    name = "utils",
    srcs = [
        "utils.cc",
    ],
    hdrs = [
        "utils.h",
    ],
    deps = [
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "utils_test",
    srcs = [
        "utils_test.cc",
    ],
    deps = [
        ":utils",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/testing:status_is_fork",
        "@com_github_google_shell//shell_encryption/testing:status_testing",
    ],
)
