load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "pir_params",
    srcs = ["pir_params.cc"],
    hdrs = ["pir_params.h"],
    deps = [
        "//private_membership/rlwe/v2/crypto:encryption",
        "//private_membership/rlwe/v2/proto:pir_cc_proto",
        "@com_google_absl//absl/numeric:bits",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)
