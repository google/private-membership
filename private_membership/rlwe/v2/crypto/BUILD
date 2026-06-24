load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "polynomial",
    srcs = ["polynomial.cc"],
    hdrs = ["polynomial.h"],
    deps = [
        "//private_membership/rlwe/v2/proto:rlwe_cc_proto",
        "@com_github_google_shell//shell_encryption:constants",
        "@com_github_google_shell//shell_encryption:dft_transformations",
        "@com_github_google_shell//shell_encryption:int256",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:ntt_parameters",
        "@com_github_google_shell//shell_encryption:statusor_fork",
        "@com_github_google_shell//shell_encryption/rns:crt_interpolation",
        "@com_github_google_shell//shell_encryption/rns:rns_modulus",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings:string_view",
    ],
)

cc_library(
    name = "encryption",
    srcs = ["encryption.cc"],
    hdrs = ["encryption.h"],
    deps = [
        ":polynomial",
        "@com_github_google_shell//shell_encryption:integral_types",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:sample_error",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_github_google_shell//shell_encryption/sampler:uniform_ternary",
        "@com_google_absl//absl/numeric:bits",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)
