# Batch RLWE-based Private Set Membership

## Batch Queries

In this directory, we present libraries to perform batch Private Set Membership (PSM)
queries. Batch queries are typically a large set of membership queries that
need to be performed at one time.

## Library

The libraries in this directory provide the core cryptographic operations
necessary for performing a batch of PSM queries. The libraries are purposedly
left a lower level so that they may be used in large-scale data processing
workflows (such as Apache Beam).
