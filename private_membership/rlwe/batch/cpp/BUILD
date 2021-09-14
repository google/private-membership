package(default_visibility = ["//:__subpackages__"])

licenses(["notice"])

cc_library(
    name = "constants",
    hdrs = [
        "constants.h",
    ],
    deps = [
        "@com_github_google_shell//shell_encryption:context",
        "@com_github_google_shell//shell_encryption:integral_types",
        "@com_github_google_shell//shell_encryption:montgomery",
    ],
)

cc_library(
    name = "context",
    srcs = [
        "context.cc",
    ],
    hdrs = [
        "context.h",
    ],
    deps = [
        ":constants",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_shell//shell_encryption:context",
        "@com_github_google_shell//shell_encryption:integral_types",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)

cc_library(
    name = "encoding",
    srcs = [
        "encoding.cc",
    ],
    hdrs = [
        "encoding.h",
    ],
    deps = [
        ":constants",
        ":padding",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_shell//shell_encryption:context",
        "@com_github_google_shell//shell_encryption:galois_key",
        "@com_github_google_shell//shell_encryption:integral_types",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:polynomial",
        "@com_github_google_shell//shell_encryption:serialization_cc_proto",
        "@com_github_google_shell//shell_encryption:transcription",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_github_google_shell//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "padding",
    srcs = [
        "padding.cc",
    ],
    hdrs = [
        "padding.h",
    ],
    deps = [
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_shell//shell_encryption:context",
        "@com_github_google_shell//shell_encryption:integral_types",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:transcription",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_github_google_shell//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "prng",
    srcs = [
        "prng.cc",
    ],
    hdrs = [
        "prng.h",
    ],
    deps = [
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_shell//shell_encryption:context",
        "@com_github_google_shell//shell_encryption:integral_types",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:transcription",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_github_google_shell//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "encoding_test",
    srcs = ["encoding_test.cc"],
    deps = [
        ":constants",
        ":context",
        ":encoding",
        ":padding",
        ":prng",
        ":test_helper",
        "//private_membership/rlwe/batch/proto:client_cc_proto",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption:galois_key",
        "@com_github_google_shell//shell_encryption:polynomial",
        "@com_github_google_shell//shell_encryption:serialization_cc_proto",
        "@com_google_absl//absl/status:statusor",
    ],
)

cc_test(
    name = "padding_test",
    srcs = ["padding_test.cc"],
    deps = [
        ":padding",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "prng_test",
    srcs = ["prng_test.cc"],
    deps = [
        ":prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_google_absl//absl/status:statusor",
    ],
)

cc_test(
    name = "context_test",
    srcs = ["context_test.cc"],
    deps = [
        ":context",
        ":test_helper",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status:statusor",
    ],
)

cc_library(
    name = "test_helper",
    testonly = 1,
    srcs = ["test_helper.cc"],
    hdrs = ["test_helper.h"],
    deps = [
        ":constants",
        ":prng",
        "//private_membership/rlwe/batch/cpp/client",
        "//private_membership/rlwe/batch/proto:client_cc_proto",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_github_google_shell//shell_encryption:montgomery",
        "@com_github_google_shell//shell_encryption:polynomial",
        "@com_github_google_shell//shell_encryption:serialization_cc_proto",
        "@com_github_google_shell//shell_encryption:symmetric_encryption",
        "@com_github_google_shell//shell_encryption/prng",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)
