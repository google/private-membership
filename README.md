# Private Set Membership (PSM)

## Problem Statement

Private Set Membership (PSM) is a cryptographic protocol that allows clients to privately query whether the client's identifier is a member of a set of identifiers held by a server in a privacy-preserving manner.

At a high level, PSM provides the following privacy guarantees:

*   The server does not learn the client's queried identifier in the
    plaintext.
*   The server does not learn whether the client's query results in a
    membership or non-membership determination.
*   The querying client does not learn any information about the set of
    identifiers that are stored by the server beyond whether the querying
    client's identifier is a member or not of the server-held set of identifiers.
    In other words, the querying client learns the bare minimum amount of
    information which is only the answer of the membership query.

## Dependencies

The Private Set Membership library requires the following dependencies:

*   [Abseil](https://github.com/abseil/abseil-cpp) for C++ common libraries.

*   [Bazel](https://github.com/bazelbuild/bazel) for building the library.

*   [BoringSSL](https://github.com/google/boringssl) for underlying
    cryptographic operations.

*   [GFlag](https://github.com/gflags/gflags) for flags. Needed to use glog.

*   [GLog](https://github.com/google/glog) for logging.

*   [Google Test](https://github.com/google/googletest) for unit testing the
    library.

*   [Protocol Buffers](https://github.com/google/protobuf) for data
    serialization.

*   [Shell](https://github.com/google/shell-encryption) for fully homomorphic encryption.

*   [Tink](https://github.com/google/tink) for cryptographic PRNGs.

## How to build

In order to run this library, you need to install Bazel, if you don't have
it already.
[Follow the instructions for your platform on the Bazel website. Make sure you
 are installing version 4.2.1 or above.]
(https://docs.bazel.build/versions/master/install.html)

You also need to install Git, if you don't have it already.
[Follow the instructions for your platform on the Git website.](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

Once you've installed Bazel and Git, open a Terminal and clone the repository into a local folder.

Navigate into the `private-membership` folder you just created, and build the
library and dependencies using Bazel. Note, the library must be built using C++17.

```bash
cd private-membership
bazel build ... --cxxopt='-std=c++17'
```

You may also run all tests (recursively) using the following command:

```bash
bazel test ... --cxxopt='-std=c++17'
```

## Disclaimers

This is not an officially supported Google product. The software is provided as-is without any guarantees or warranties, express or implied.
