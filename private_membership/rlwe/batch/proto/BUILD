load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//:__subpackages__"])

licenses(["notice"])

proto_library(
    name = "shared_proto",
    srcs = ["shared.proto"],
    visibility = ["//visibility:public"],
    deps = [
        "@com_github_google_shell//shell_encryption:serialization_proto",
    ],
)

cc_proto_library(
    name = "shared_cc_proto",
    visibility = ["//visibility:public"],
    deps = [
        ":shared_proto",
    ],
)

java_proto_library(
    name = "shared_java_proto",
    visibility = ["//visibility:public"],
    deps = [
        ":shared_proto",
    ],
)

proto_library(
    name = "client_proto",
    srcs = ["client.proto"],
    visibility = ["//visibility:public"],
    deps = [
        ":shared_proto",
        "@com_github_google_shell//shell_encryption:serialization_proto",
    ],
)

cc_proto_library(
    name = "client_cc_proto",
    visibility = ["//visibility:public"],
    deps = [
        ":client_proto",
    ],
)

java_proto_library(
    name = "client_java_proto",
    visibility = ["//visibility:public"],
    deps = [
        ":client_proto",
    ],
)

proto_library(
    name = "server_proto",
    srcs = ["server.proto"],
    visibility = ["//visibility:public"],
    deps = [
        ":shared_proto",
        "@com_github_google_shell//shell_encryption:serialization_proto",
    ],
)

cc_proto_library(
    name = "server_cc_proto",
    visibility = ["//visibility:public"],
    deps = [
        ":server_proto",
    ],
)

java_proto_library(
    name = "server_java_proto",
    visibility = ["//visibility:public"],
    deps = [
        ":server_proto",
    ],
)
