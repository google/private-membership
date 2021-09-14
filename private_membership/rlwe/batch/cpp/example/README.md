# Example for Batch RLWE-based Private Set Membership

This directory contains a binary showing a simple example of the entire batch
pipeline from client key and request generation, server response and client decryption.

To run the pipeline, you can run the following command in this directory:

```bash
bazel run -c opt :example --cxxopt='-std=c++17'
```
