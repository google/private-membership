# RLWE-based Private Set Membership

## Cryptographic Protocol

This protocol for Private Set Membership (PSM) utilizes two cryptographic primitives:

*  **Fully Homomorphic Encryption**: FHE enables computation over encrypted data without ever needing to decrypt the data.
*  **Oblivious Hashing**: This is a cryptographic protocol that enables a sender and a receiver to jointly compute a hash, H(K, x), where the sender holds the key K and the receiver holds the input x. The hash is given to the receiver. The sender's key K and receiver's input x remain private from the other party.

To prepare for PSM queries, the server does the following:

1. **Hash Set**: The server generates a hash key K. For all items s in the set, the server computes a hash of the item H(K, s).
2. **Bucketize**: A prefix of the hash is used to bucketize the set. The remainder of the hash is stored in the bucket's contents. For example, if bucket identifiers are 10 bits long, the first 10 bits of H(K, s) determine the bucket and the remaining bits are stored in the associated bucket.

When querying for an item x, the client does the following:

1. **Oblivious Hashing**: The client and server execute the oblivious hashing protocol so that the client receives H(K, x). The server does not learn the queried item x and the client does not learn the key K.
2. **Compute Bucket Identifier**: The client computes the bucket identifier and contents of the item in the bucket using H(K, x) in the same way as the server.
3. **Fully Homomorphic Encryption**:  The client encrypts the bucket identifier using fully homomorphic encryption.
4. **Bucket Retrieval**: The server uses the FHE encrypted bucket identifier to compute the associated bucket that is returned to the client. The server does this without learning the identity of the queried bucket.
5. **Client-side Matching**: The client decrypts the bucket contents. Then, the client matches the hash of the queried item with the bucket contents to generate the final membership response. This match occurs exclusively on the client-side and the server does not learn the membership result.
