load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "oprf_constants",
    hdrs = ["oprf_constants.h"],
    deps = [
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "oprf_utils",
    srcs = ["oprf_utils.cc"],
    hdrs = ["oprf_utils.h"],
    deps = [
        ":oprf_constants",
        "//private_membership/oprf/proto:oprf_cc_proto",
        "//private_membership/rlwe/client/cpp/internal:aes_ctr_256_no_auth",
        "//private_membership/rlwe/client/cpp/internal:aes_ctr_256_with_fixed_iv",
        "//private_membership/rlwe/client/cpp/internal:aes_gcm_256",
        "//private_membership/rlwe/client/cpp/internal:cipher",
        "//private_membership/rlwe/client/cpp/internal:crypto_utils",
        "@boringssl//:ssl",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:bn_util",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_util",
        "@com_github_google_private_join_and_compute//private_join_and_compute/util:status_includes",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@tink_cc//:deterministic_aead",
        "@tink_cc//subtle",
        "@tink_cc//subtle:aes_siv_boringssl",
        "@tink_cc//subtle:common_enums",
        "@tink_cc//subtle:hkdf",
        "@tink_cc//util:secret_data",
    ],
)
