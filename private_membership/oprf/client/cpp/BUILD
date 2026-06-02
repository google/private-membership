load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "oprf_client",
    srcs = ["oprf_client.cc"],
    hdrs = ["oprf_client.h"],
    deps = [
        "//private_membership/oprf/proto:oprf_cc_proto",
        "//private_membership/oprf/shared:oprf_utils",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:bn_util",
        "@com_github_google_private_join_and_compute//private_join_and_compute/crypto:ec_commutative_cipher",
        "@com_github_google_private_join_and_compute//private_join_and_compute/util:status_includes",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)
