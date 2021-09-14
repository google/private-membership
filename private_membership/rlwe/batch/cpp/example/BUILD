licenses(["notice"])

package(default_visibility = ["//:__subpackages__"])

cc_binary(
    name = "example",
    srcs = ["example.cc"],
    deps = [
        "//private_membership/rlwe/batch/cpp/client",
        "//private_membership/rlwe/batch/cpp/server",
        "//private_membership/rlwe/batch/proto:client_cc_proto",
        "//private_membership/rlwe/batch/proto:server_cc_proto",
        "//private_membership/rlwe/batch/proto:shared_cc_proto",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/time",
    ],
)
