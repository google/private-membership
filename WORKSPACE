load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# Tink
http_archive(
    name = "tink_base",
    urls = ["https://github.com/google/tink/archive/master.zip"],
    strip_prefix = "tink-master/",
    sha256 = "694230dea9d064fa423a81bc21c0aa7bddafd1cc8b017594756cda301f078be3",
)

http_archive(
    name = "tink_cc",
    urls = ["https://github.com/google/tink/archive/master.zip"],
    strip_prefix = "tink-master/cc",
    sha256 = "694230dea9d064fa423a81bc21c0aa7bddafd1cc8b017594756cda301f078be3",
)

load("@tink_base//:tink_base_deps.bzl", "tink_base_deps")
tink_base_deps()

load("@tink_base//:tink_base_deps_init.bzl", "tink_base_deps_init")
tink_base_deps_init()

load("@tink_cc//:tink_cc_deps.bzl", "tink_cc_deps")
tink_cc_deps()

load("@tink_cc//:tink_cc_deps_init.bzl", "tink_cc_deps_init")
tink_cc_deps_init()

# rules_cc defines rules for generating C++ code from Protocol Buffers.
http_archive(
    name = "rules_cc",  # 2021-06-07T16:41:49Z
    sha256 = "b295cad8c5899e371dde175079c0a2cdc0151f5127acc92366a8c986beb95c76",
    strip_prefix = "rules_cc-daf6ace7cfeacd6a83e9ff2ed659f416537b6c74",
    urls = ["https://github.com/bazelbuild/rules_cc/archive/daf6ace7cfeacd6a83e9ff2ed659f416537b6c74.zip"],
)

load("@rules_cc//cc:repositories.bzl", "rules_cc_dependencies")
rules_cc_dependencies()

# rules_proto defines abstract rules for building Protocol Buffers.
# https://github.com/bazelbuild/rules_proto
http_archive(
    name = "rules_proto",
    sha256 = "e0cab008a9cdc2400a1d6572167bf9c5afc72e19ee2b862d18581051efab42c9",
    strip_prefix = "rules_proto-c0b62f2f46c85c16cb3b5e9e921f0d00e3101934",
    urls = [
        "https://mirror.bazel.build/github.com/bazelbuild/rules_proto/archive/c0b62f2f46c85c16cb3b5e9e921f0d00e3101934.tar.gz",
        "https://github.com/bazelbuild/rules_proto/archive/c0b62f2f46c85c16cb3b5e9e921f0d00e3101934.tar.gz",
    ],
)
load("@rules_proto//proto:repositories.bzl", "rules_proto_dependencies", "rules_proto_toolchains")
rules_proto_dependencies()
rules_proto_toolchains()

# Install gtest: tag = "release-1.10.0"
http_archive(
    name = "com_github_google_googletest",
    urls = [
        "https://github.com/google/googletest/archive/703bd9caab50b139428cea1aaff9974ebee5742e.tar.gz",
    ],
    sha256 = "d17b1b83a57b3933565a6d0616fe261107326d47de20288d0949ed038e1c342d",
    strip_prefix = "googletest-703bd9caab50b139428cea1aaff9974ebee5742e",
)

# Install abseil-cpp.
http_archive(
    name = "com_google_absl",
    sha256 = "ebe2ad1480d27383e4bf4211e2ca2ef312d5e6a09eba869fd2e8a5c5d553ded2",
    strip_prefix = "abseil-cpp-20200923.3",
    urls = [
         "https://github.com/abseil/abseil-cpp/archive/20200923.3.tar.gz",
    ],
)

# BoringSSL
git_repository(
    name = "boringssl",
    commit = "67ffb9606462a1897d3a5edf5c06d329878ba600",  # https://boringssl.googlesource.com/boringssl/+/refs/heads/master-with-bazel
    remote = "https://boringssl.googlesource.com/boringssl",
    shallow_since = "1585767053 +0000"
)

# Logging
http_archive(
    name = "com_github_google_glog",
    urls = ["https://github.com/google/glog/archive/96a2f23dca4cc7180821ca5f32e526314395d26a.zip"],
    strip_prefix = "glog-96a2f23dca4cc7180821ca5f32e526314395d26a",
    sha256 = "6281aa4eeecb9e932d7091f99872e7b26fa6aacece49c15ce5b14af2b7ec050f",
)

# gflags, needed for glog
http_archive(
    name = "com_github_gflags_gflags",
    urls = ["https://github.com/gflags/gflags/archive/v2.2.2.tar.gz"],
    sha256 = "34af2f15cf7367513b352bdcd2493ab14ce43692d2dcd9dfc499492966c64dcf",
    strip_prefix = "gflags-2.2.2",
)

# Install shell-encryption.
http_archive(
    name = "com_github_google_shell",
    urls = [
        "https://github.com/google/shell-encryption/archive/4e2c0dded993ec96eb319453e35e7d95ce689b45.tar.gz",
    ],
    strip_prefix = "shell-encryption-4e2c0dded993ec96eb319453e35e7d95ce689b45",
    sha256 = "a6404a1fdc794c8e5424ed7233a925872d2abd7b08462a366f4a6485ea0747e7",
)
