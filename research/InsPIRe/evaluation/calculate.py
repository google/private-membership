from math import ceil, log2, floor, sqrt, pi
from matplotlib import pyplot as plt
from enum import Enum


# Other possible ways for future work:
# - in the medium payload case, using different plaintext moduli for the first and second layers -> would this require separate packing keys for the first and second layers?
# - in the medium payload case, using diffenent rlwe parameters for the first and second layers -> this would potentially require separate packing keys for the first and second layers

def ours_stats(ell_prime, s):

  log_d1 = 10
  log_d2 = 11

  log_q1 = 32
  log_q1_tilde = 28

  log_q2 = 56
  log_q2_1_tilde = 28
  log_q2_2_tilde = 20

  log_z = 19

  log_N =8
  log_p = 15

  db_size = ell_prime * s * log_N

  log_ell_1 = ceil(log2(ell_prime) / 2)
  log_ell_2 = floor(log2(ell_prime) / 2)

  ell_1 = 2**log_ell_1
  ell_2 = 2**log_ell_2

  d1 = 2**log_d1
  d2 = 2**log_d2

  kappa = ceil(log_q1_tilde / log_p)
  rho = ceil(kappa * (d1+1) / d2)

  key_size = log_d2 * ceil(log_q2 / log_z) * d2 * log_q2
  query_size = ell_1 * log_q1 + ell_2 * log_q2
  response_size = rho * (d2 * log_q2_1_tilde + s * log_q2_2_tilde)
  
  return {
      "db_size": db_size,
      "ell_prime": ell_prime,
      "s": s,
      "key_size": key_size,
      "query_size": query_size,
      "response_size": response_size,
      "total_size": key_size + query_size + response_size
  }

def ypir_stats(ell_prime, s):

  log_d1 = 10
  log_d2 = 11

  log_q1 = 32
  log_q1_tilde = 28

  log_q2 = 56
  log_q2_1_tilde = 28
  log_q2_2_tilde = 20

  log_z = 19

  log_N =8
  log_p = 15

  db_size = ell_prime * s * log_N

  log_ell_1 = ceil(log2(ell_prime) / 2)
  log_ell_2 = floor(log2(ell_prime) / 2)

  ell_1 = 2**log_ell_1
  ell_2 = 2**log_ell_2

  d1 = 2**log_d1
  d2 = 2**log_d2

  kappa = ceil(log_q1_tilde / log_p)
  rho = ceil(kappa * (d1+1) / d2)

  key_size = log_d2 * ceil(log_q2 / log_z) * d2 * log_q2
  query_size = ell_1 * log_q1 + ell_2 * log_q2
  response_size = s * rho * (d2 * log_q2_1_tilde + log_q2_2_tilde)
  
  return {
      "db_size": db_size,
      "ell_prime": ell_prime,
      "s": s,
      "key_size": key_size,
      "query_size": query_size,
      "response_size": response_size,
      "total_size": key_size + query_size + response_size
  }

def plotting():
  # plot total size as a function of s
  ell_prime = 10**9
  s_s = list(range(1, 100))
  ypir_sizes_KB = []
  ours_sizes_KB = []
  for s in s_s:
    stats = ypir_stats(ell_prime, s)
    ypir_sizes_KB.append(stats["total_size"]/8096)
    
    stats = ours_stats(ell_prime, s)
    ours_sizes_KB.append(stats["total_size"]/8096)
    
  plt.xlabel("s")
  plt.ylabel("total size (KB)")
  plt.plot(s_s, ypir_sizes_KB)
  plt.plot(s_s, ours_sizes_KB)
  plt.savefig("total_size.png")
  
def b_to_KB(x):
  return x//8096

def B_to_KB(x):
  return x/1024

class ProtocolType(Enum):
  def __str__(self):
    return self.name

  # first inner product pow2, second inner product over the ring then pack
  DOUBLE_YPIR_SMALL_PAYLOAD_WITHOUT_ONLINE_PACKING = 1

  # first inner product pow2, second inner product over the ring then pack
  DOUBLE_YPIR_SMALL_PAYLOAD_WITH_ONLINE_PACKING = 2

  # first inner product over the ring, pack,
  # second inner product over the ring then pack
  DOUBLE_YPIR_MEDIUM_PAYLOAD_WITHOUT_ONLINE_PACKING = 3

  # first inner product over the ring, pack,
  # second inner product over the ring then pack
  DOUBLE_YPIR_MEDIUM_PAYLOAD_WITH_ONLINE_PACKING = 4

  # inner product over the ring then pack
  SIMPLE_YPIR = 5 

def compute_costs_of_simple_ypir(N, payload_bits_input):

  valid_outputs = []

  for param_set in [
      {
          "max_log_q": 29,
          "log_d": 10,
          "sigma": 3.2 * sqrt(2 * pi),
      },
      {
          "max_log_q": 56,
          "log_d": 11,
          "sigma": 6.4 * sqrt(2 * pi),
      },
  ]:

    # RLWE parameters (for second layer)
    max_log_q = param_set["max_log_q"]
    log_d = param_set["log_d"]
    sigma = param_set["sigma"]
    d = 2**log_d

    for reduction in range(16):
      log_q = max_log_q - reduction

      for log_p in range(2, 17):
        s_input_based_on_params = ceil(payload_bits_input / log_p)
        log_s_input_based_on_params = ceil(log2((s_input_based_on_params)))
        if log_s_input_based_on_params > log_q:
          continue

        for log_max_packed in [log_d - 1, log_d]:
          max_packed = 2**log_max_packed

          log_q_prime = 2 + log_p + ceil(log2(d+1))

          s_working = ceil(sqrt(N*max_packed*log_q / ((d+max_packed)*log_q_prime)))
          s_working = s_input_based_on_params * ceil(s_working / s_input_based_on_params)
          log_s_working = ceil(log2(s_working))

          dim1 = ceil(N / s_working)

          def stats(B):
            number_of_decomp_digits = ceil(log_q / log2(B))
            if log_max_packed == log_d:
              key_size = 2 * number_of_decomp_digits * (d * log_q)
            else:
              key_size = number_of_decomp_digits * (d * log_q)

            query_size = dim1 * log_q
            
            ##### DoubleYPIR with medium payload with online packing #####
            protocol_type = ProtocolType.SIMPLE_YPIR
            response_size = ceil(s_working / max_packed) * (d + max_packed) * log_q_prime

            total_size = key_size + query_size + response_size
            return {
                "protocol_type": protocol_type.__str__(),
                "comm": {
                    "key_size": key_size,
                    "query_size": {
                        "total": query_size,
                    },
                    "response_size": {
                        "total": response_size,
                    },
                    "total_size": total_size,
                    "key_size_KB": b_to_KB(key_size),
                    "query_size_KB": b_to_KB(query_size),
                    "response_size_KB": b_to_KB(response_size),
                    "total_size_KB": b_to_KB(total_size),
                },
            }

          def is_valid_B(B):
            first_layer_rhs = (
                sqrt(d) * sqrt(dim1) * (2**log_p)
                + sqrt(max_packed) * ceil(log_q / log2(B)) * B * sqrt(d)
            ) * sigma
            first_layer_check = log_q - log_p - 1 - 1 > log2(first_layer_rhs)
            
            return first_layer_check

          # find largest B that satisfies the condition
          base = 2
          while is_valid_B(2*base):
            base *= 2

          if is_valid_B(base):
            digits = ceil(log_q / log2(base))
            valid_outputs.append({
                "params": {
                    "N": N,
                    "payload_bits_input": payload_bits_input,
                    "payload_bytes_input": floor(payload_bits_input / 8),
                    "log_d": log_d,
                    "log_q": log_q,
                    "log_q_prime": log_q_prime,
                    "log_p": log_p,
                    "dim1": dim1,
                    "log_s_working": log_s_working,
                    "s_working": s_working,
                    "N_working": dim1,
                    "log_max_packed": log_max_packed,
                    "max_packed": max_packed,
                    "log_B": ceil(log2(base)),
                    "digits": digits,
                },
                "stats": stats(base),
            })
  return valid_outputs

def compute_costs_of_double_ypir_medium_payload(N, payload_bits_input):

  valid_outputs = []

  for param_set in [
      # {
      #     "max_log_q": 29,
      #     "log_d": 10,
      #     "sigma": 3.2 * sqrt(2 * pi),
      # },
      {
          "max_log_q": 56,
          "log_d": 11,
          "sigma": 6.4 * sqrt(2 * pi),
      },
  ]:

    # RLWE parameters (for second layer)
    max_log_q = param_set["max_log_q"]
    log_d = param_set["log_d"]
    sigma = param_set["sigma"]
    d = 2**log_d

    max_reduction = 16
    for reduction in range(max_reduction):
      log_q = max_log_q - reduction

      for log_p in range(2, 17):
        s_input_based_on_params = ceil(payload_bits_input / log_p)
        log_s_input_based_on_params = ceil(log2((s_input_based_on_params)))
        if log_s_input_based_on_params > log_q:
          continue

        for log_resize_factor in range(log_d - log_s_input_based_on_params):
          resize_factor = 2**log_resize_factor
          log_s_working = log_s_input_based_on_params + log_resize_factor
          assert log_s_working <= log_d
          s_working = 2**log_s_working

          for log2_gamma_0 in range(1, log_d + 1):
            for log2_gamma_1 in range(1, log_d + 1):
              for log2_gamma_2 in range(1, log_d + 1):

                gamma_0 = 2**log2_gamma_0
                gamma_1 = 2**log2_gamma_1
                gamma_2 = 2**log2_gamma_2
                
                log_q_prime = 28 # 2 + log_p + ceil(log2(d + 1))
                tau = log_q_prime / log_p
                
                N_working = ceil(N / resize_factor)

                dim0 = ceil(sqrt(N_working))
                dim1 = floor(sqrt(N_working))

                for online_packing in [True, False]:
                  
                  def stats(B):
                    number_of_decomp_digits = ceil(log_q / log2(B))
                    if log2_gamma_0 == log_d:
                      key_size = 2 * number_of_decomp_digits * (d * log_q)
                    else:
                      key_size = number_of_decomp_digits * (d * log_q)

                    first_dim_query = dim0 * log_q
                    second_dim_query = dim1 * log_q
                    query_size = first_dim_query + second_dim_query
                    
                    ##### InsPIRe with online packing #####
                    if online_packing:
                      protocol_type = ProtocolType.DOUBLE_YPIR_MEDIUM_PAYLOAD_WITH_ONLINE_PACKING
                      packed_rlwes = (ceil(num_partitions * d * tau / gamma_0) \
                                      + ceil(num_partitions * s_working * tau / gamma_0)) * (d + gamma_0) * log_q_prime
                      unpacked_lwes = 0
                    
                    ##### InsPIRe without online packing ###
                    else:
                      protocol_type = ProtocolType.DOUBLE_YPIR_MEDIUM_PAYLOAD_WITHOUT_ONLINE_PACKING
                      packed_rlwes = ceil(num_partitions * d * tau / gamma_0) * (d + gamma_0) * log_q_prime 
                      unpacked_lwes = ceil(num_partitions * s_working * tau) * (d + 1) * log_q_prime

                    response_size = packed_rlwes + unpacked_lwes
                    total_size = key_size + query_size + response_size
                    return {
                        "protocol_type": protocol_type.__str__(),
                        "comm": {
                            "key_size": key_size,
                            "query_size": {
                                "total": query_size,
                                "first_dim_query": first_dim_query,
                                "second_dim_query": second_dim_query,
                            },
                            "response_size": {
                                "total": response_size,
                                "packed_lwes": packed_rlwes,
                                "unpacked_lwes": unpacked_lwes,
                            },
                            "total_size": total_size,
                            "key_size_KB": b_to_KB(key_size),
                            "query_size_KB": b_to_KB(query_size),
                            "response_size_KB": b_to_KB(response_size),
                            "total_size_KB": b_to_KB(total_size),
                        },
                    }

                  def is_valid_B(B):
                    first_layer_rhs = (
                        sqrt(d) * sqrt(dim0) * (2**log_p)
                        + sqrt(max_packed) * ceil(log_q / log2(B)) * B * sqrt(d)
                    ) * sigma
                    first_layer_check = log_q - log_p - 1 - 1 > log2(first_layer_rhs)

                    second_layer_rhs = (
                        sqrt(d) * sqrt(dim1) * (2**log_p)
                        + sqrt(max_packed) * ceil(log_q / log2(B)) * B * sqrt(d)
                    ) * sigma
                    
                    second_layer_check = log_q - log_p - 1 - 1 > log2(second_layer_rhs)

                    return first_layer_check and second_layer_check

                  # find largest B that satisfies the condition
                  base = 2
                  while is_valid_B(2*base):
                    base *= 2

                  if is_valid_B(base):
                    digits = ceil(log_q / log2(base))
                    valid_outputs.append({
                        "params": {
                            "N": N,
                            "payload_bits_input": payload_bits_input,
                            "payload_bytes_input": floor(payload_bits_input / 8),
                            "log_d": log_d,
                            "log_q": log_q,
                            "log_q_prime": log_q_prime,
                            "log_p": log_p,
                            "tau": tau,
                            "dim1": dim0,
                            "dim2": dim1,
                            "log_s_working": log_s_working,
                            "s_working": s_working,
                            "num_partitions": num_partitions,
                            "N_working": N_working,
                            "log2_gamma_0": log2_gamma_0,
                            "gamma_0": gamma_0,
                            "log_resize_factor": log_resize_factor,
                            "resize_factor": 2**log_resize_factor,
                            "log_B": ceil(log2(base)),
                            "digits": digits,
                        },
                        "stats": stats(base),
                    })
  return valid_outputs

def compute_costs_of_double_ypir_small_payload(N, payload_bits_input):

  valid_outputs = []

  for param_set in [
      {
          "max_log_q": 29,
          "log_d": 10,
          "sigma": 3.2 * sqrt(2 * pi),
      },
      {
          "max_log_q": 56,
          "log_d": 11,
          "sigma": 6.4 * sqrt(2 * pi),
      },
  ]:

    # RLWE parameters (for second layer)
    max_log_q = param_set["max_log_q"]
    log_d = param_set["log_d"]
    sigma = param_set["sigma"]
    d = 2**log_d

    # LWE parameters (for first layer)
    lwe_n = 1024
    lwe_log_q = 32
    lwe_log_p = 8
    lwe_noise_std = 11 * sqrt(2*pi)
  
    for reduction in range(16):
      log_q = max_log_q - reduction

      for log_p in range(2, 17):
        s_input_based_on_params = ceil(payload_bits_input / log_p)
        log_s_input_based_on_params = ceil(log2((s_input_based_on_params)))
        if log_s_input_based_on_params > log_q:
          continue

        for log_resize_factor in range(log_d - log_s_input_based_on_params):
          resize_factor = 2**log_resize_factor
          log_s_working = log_s_input_based_on_params + log_resize_factor
          assert log_s_working <= log_d
          s_working = 2**log_s_working

          num_partitions = s_working

          for log_max_packed in range(1, log_d + 1):
            max_packed = 2**log_max_packed

            log_q_prime = 2 + log_p + ceil(log2(d + 1))
            lwe_log_q_prime = 2 + lwe_log_p + ceil(log2(lwe_n+1))
            
            tau = lwe_log_q_prime / log_p
            
            N_working = ceil(N / resize_factor)

            dim0 = ceil(sqrt(N_working * log_q / lwe_log_q))
            dim1 = ceil(N_working / dim0)

            for online_packing in [True, False]:
              
              def stats(B):
                number_of_decomp_digits = ceil(log_q / log2(B))
                if log_max_packed == log_d:
                  key_size = 2 * number_of_decomp_digits * (d * log_q)
                else:
                  key_size = number_of_decomp_digits * (d * log_q)

                first_dim_query = dim0 * lwe_log_q
                second_dim_query = dim1 * log_q
                query_size = first_dim_query + second_dim_query
                
                ##### DoubleYPIR with small payload with online packing #####
                if online_packing:
                  protocol_type = ProtocolType.DOUBLE_YPIR_SMALL_PAYLOAD_WITH_ONLINE_PACKING
                  packed_rlwes = (
                      ceil(num_partitions * d * tau / max_packed) + \
                      ceil(num_partitions * tau / max_packed) \
                                      ) * (d + max_packed) * log_q_prime
                  unpacked_lwes = 0
                
                ##### DoubleYPIR with small payload without online packing ###
                else:
                  protocol_type = ProtocolType.DOUBLE_YPIR_SMALL_PAYLOAD_WITHOUT_ONLINE_PACKING
                  packed_rlwes = ceil(num_partitions * d * tau / max_packed) * (d + max_packed) * log_q_prime 
                  unpacked_lwes = ceil(num_partitions * s_working * tau) * (d + 1) * log_q_prime

                response_size = packed_rlwes + unpacked_lwes
                total_size = key_size + query_size + response_size
                return {
                    "protocol_type": protocol_type.__str__(),
                    "comm": {
                        "key_size": key_size,
                        "query_size": {
                            "total": query_size,
                            "first_dim_query": first_dim_query,
                            "second_dim_query": second_dim_query,
                        },
                        "response_size": {
                            "total": response_size,
                            "packed_lwes": packed_rlwes,
                            "unpacked_lwes": unpacked_lwes,
                        },
                        "total_size": total_size,
                        "key_size_KB": b_to_KB(key_size),
                        "query_size_KB": b_to_KB(query_size),
                        "response_size_KB": b_to_KB(response_size),
                        "total_size_KB": b_to_KB(total_size),
                    },
                }

              def is_valid_B(B):
                first_layer_rhs = (
                    sqrt(dim0) * (2**lwe_log_p)
                ) * lwe_noise_std
                first_layer_check = lwe_log_q - lwe_log_p - 1 - 1 > log2(first_layer_rhs)

                second_layer_rhs = (
                    sqrt(d) * sqrt(dim1) * (2**log_p)
                    + sqrt(max_packed) * ceil(log_q / log2(B)) * B * sqrt(d)
                ) * sigma
                
                second_layer_check = log_q - log_p - 1 - 1 > log2(second_layer_rhs)

                return first_layer_check and second_layer_check

              # find largest B that satisfies the condition
              base = 2
              while is_valid_B(2*base):
                base *= 2

              if is_valid_B(base):
                digits = ceil(log_q / log2(base))
                valid_outputs.append({
                    "params": {
                        "N": N,
                        "payload_bits_input": payload_bits_input,
                        "payload_bytes_input": floor(payload_bits_input / 8),
                        "log_d": log_d,
                        "log_q": log_q,
                        "log_q_prime": log_q_prime,
                        "log_p": log_p,
                        "tau": tau,
                        "dim1": dim0,
                        "dim2": dim1,
                        "log_s_working": log_s_working,
                        "s_working": s_working,
                        "num_partitions": num_partitions,
                        "N_working": N_working,
                        "log_max_packed": log_max_packed,
                        "max_packed": max_packed,
                        "log_resize_factor": log_resize_factor,
                        "resize_factor": 2**log_resize_factor,
                        "log_B": ceil(log2(base)),
                        "digits": digits,
                    },
                    "stats": stats(base),
                })
  return valid_outputs

def compute_costs_of_inspire(N_input, payload_bits_input):

  valid_outputs = []

  for param_set in [
      # {
      #     "max_log_q": 29,
      #     "log_d": 10,
      #     "sigma": 3.2 * sqrt(2 * pi),
      # },
      {
          "max_log_q": 56,
          "log_d": 11,
          "sigma": 6.4 * sqrt(2 * pi),
      },
  ]:

    # RLWE parameters (for second layer)
    max_log_q = param_set["max_log_q"]
    log_d = param_set["log_d"]
    sigma = param_set["sigma"]
    d = 2 ** log_d

    max_reduction = 4
    for log_q_0 in range(max_log_q, max_log_q-max_reduction, -1):
      for log_q_1 in range(max_log_q, max_log_q-max_reduction, -1):
          print(f"log_q_0: {log_q_0}, log_q_1: {log_q_1}")
          for log_p_0 in range(4, 24):
            for log_p_1 in range(4, 24):
              # for log_p_2 in range(2, 17):
                for log2_gamma_0 in range(4, 7): #log_d + 1):
                  for log2_gamma_1 in range(4, log_d + 1):
                    for log2_gamma_2 in range(1, 7): #log_d + 1):

                      gamma_0 = 2**log2_gamma_0
                      gamma_1 = 2**log2_gamma_1
                      gamma_2 = 2**log2_gamma_2

                      log_q_prime_1 = 20 
                      log_q_prime_2 = 28 # 2 + log_p_0 + ceil(log2(d + 1))
                      # log_q_prime_2_0 = 2 + log_p_0 + ceil(log2(d + 1))
                      # log_q_prime_2_1 = 2 + log_p_1 + ceil(log2(d + 1))
                      # log_q_prime_2_2 = 2 + log_p_1 + ceil(log2(d + 1))

                      tau = log_q_prime_2 / log_p_1

                      N_large_items = ceil(N_input * payload_bits_input / (log_p_0 * gamma_0))
                      dim0 = ceil(sqrt(N_large_items))
                      dim1 = floor(sqrt(N_large_items))

                      def stats(B_0, B_1, B_2):
                        digits_0 = ceil(log_q_0 / log2(B_0))
                        digits_1 = ceil(log_q_1 / log2(B_1))
                        digits_2 = ceil(log_q_1 / log2(B_2))

                        param_combinations = set()
                        param_combinations.add((gamma_0, log_q_0, B_0, log_p_0))
                        param_combinations.add((gamma_1, log_q_1, B_1, log_p_1))
                        param_combinations.add((gamma_2, log_q_1, B_2, log_p_1))

                        key_size = 0
                        for (gamma, log_q, B, log_p) in param_combinations:
                          digits = ceil(log_q / log2(B))
                          if gamma == d:
                            key_size += 2 * digits * (d * log_q)
                          else:
                            key_size += digits * (d * log_q)

                        first_dim_query = dim0 * log_q_0
                        second_dim_query = dim1 * log_q_1
                        query_size = first_dim_query + second_dim_query

                        ##### InsPIRe without online packing ###
                        if log2_gamma_2 == 0:
                          protocol_type = ProtocolType.DOUBLE_YPIR_MEDIUM_PAYLOAD_WITHOUT_ONLINE_PACKING
                          packed_rlwes = ceil(ceil(d * tau) / gamma_1) * (d * log_q_prime_1 + gamma_1 * log_q_prime_2) 
                          unpacked_lwes = ceil(gamma_0 * tau) * (d * log_q_prime_1 + log_q_prime_2)
                          # packed_rlwes = ceil(ceil(d * tau) / gamma_1) * (d + gamma_1) * log_q_prime_2_1 
                          # unpacked_lwes = ceil(gamma_0 * tau) * (d + 1) * log_q_prime_2_2
                        
                        ##### InsPIRe with online packing #####
                        else:
                          protocol_type = ProtocolType.DOUBLE_YPIR_MEDIUM_PAYLOAD_WITH_ONLINE_PACKING
                          packed_rlwes = ceil(ceil(d * tau) / gamma_1) * (d * log_q_prime_1 + gamma_1 * log_q_prime_2) \
                                          + ceil(ceil(gamma_0 * tau) / gamma_2) * (d * log_q_prime_1 + gamma_2 * log_q_prime_2)
                          # packed_rlwes = ceil(ceil(d * tau) / gamma_1) * (d + gamma_1) * log_q_prime_2_1 \
                          #                 + ceil(ceil(gamma_0 * tau) / gamma_2) * (d + gamma_2) * log_q_prime_2_2
                          unpacked_lwes = 0

                        response_size = packed_rlwes + unpacked_lwes
                        total_size = key_size + query_size + response_size
                        return {
                            "protocol_type": protocol_type.__str__(),
                            "comm": {
                                "key_size": key_size,
                                "query_size": {
                                    "total": query_size,
                                    "first_dim_query": first_dim_query,
                                    "second_dim_query": second_dim_query,
                                },
                                "response_size": {
                                    "total": response_size,
                                    "packed_lwes": packed_rlwes,
                                    "unpacked_lwes": unpacked_lwes,
                                },
                                "total_size": total_size,
                                "key_size_KB": b_to_KB(key_size),
                                "query_size_KB": b_to_KB(query_size),
                                "response_size_KB": b_to_KB(response_size),
                                "total_size_KB": b_to_KB(total_size),
                            },
                        }

                      def is_valid_B_0(B_0):
                        first_layer_rhs = (
                            sqrt(d) * sqrt(dim0) * (2**log_p_0)
                            + sqrt(gamma_0) * ceil(log_q_0 / log2(B_0)) * B_0 * sqrt(d)
                        ) * sigma
                        first_layer_check = log_q_0 - log_p_0 - 1 - 1 > log2(first_layer_rhs)

                        return first_layer_check                        
                      
                      def is_valid_B_1(B_1):
                        second_layer_mask_rhs = (
                            sqrt(d) * sqrt(dim1) * (2**log_p_1)
                            + sqrt(gamma_1) * ceil(log_q_1 / log2(B_1)) * B_1 * sqrt(d)
                        ) * sigma
                        
                        second_layer_mask_check = log_q_1 - log_p_1 - 1 - 1 > log2(second_layer_mask_rhs)

                        return second_layer_mask_check

                      def is_valid_B_2(B_2):
                        second_layer_body_rhs = (
                            sqrt(d) * sqrt(dim1) * (2**log_p_1)
                            + sqrt(gamma_2) * ceil(log_q_1 / log2(B_2)) * B_2 * sqrt(d)
                        ) * sigma
                        
                        second_layer_body_check = log_q_1 - log_p_1 - 1 - 1 > log2(second_layer_body_rhs)

                        return second_layer_body_check               

                      # find largest B that satisfies the condition
                      base_0 = base_1 = base_2 = 2
                      while is_valid_B_0(2*base_0):
                        base_0 *= 2
                      while is_valid_B_1(2*base_1):
                        base_1 *= 2                          
                      while is_valid_B_2(2*base_2):
                        base_2 *= 2

                      digits_0 = ceil(log_q_0 / log2(base_0))
                      digits_1 = ceil(log_q_1 / log2(base_1))
                      digits_2 = ceil(log_q_1 / log2(base_2))

                      while base_0 > 2 and digits_0 == ceil(log_q_0 / log2(base_0/2)):
                        base_0 = base_0 / 2
                      while base_1 > 2 and digits_1 == ceil(log_q_1 / log2(base_1/2)):
                        base_1 = base_1 / 2
                      while base_2 > 2 and digits_2 == ceil(log_q_1 / log2(base_2/2)):
                        base_2 = base_2 / 2

                      if is_valid_B_1(base_0) and is_valid_B_1(base_1) and is_valid_B_2(base_2):
                        valid_outputs.append({
                            "params": {
                                "N_input": N_input,
                                "payload_bits_input": payload_bits_input,
                                "payload_bytes_input": floor(payload_bits_input / 8),
                                "log_d": log_d,
                                "log_q_0": log_q_0,
                                "log_q_1": log_q_1,
                                "log_q_prime_1": log_q_prime_1,
                                "log_q_prime_2": log_q_prime_2,
                                # "log_q_prime_2_0": log_q_prime_2_0,
                                # "log_q_prime_2_1": log_q_prime_2_1,
                                # "log_q_prime_2_2": log_q_prime_2_2,
                                "log_p_0": log_p_0,
                                "log_p_1": log_p_1,
                                "tau": tau,
                                "N_large_items": N_large_items,
                                "dim1": dim0,
                                "dim2": dim1,
                                "payload_size": log_p_0 * gamma_0,
                                "log2_gamma_0": log2_gamma_0,
                                "log2_gamma_1": log2_gamma_1,
                                "log2_gamma_2": log2_gamma_2,
                                "gamma_0": gamma_0,
                                "gamma_1": gamma_1,
                                "gamma_2": gamma_2,
                                "log_B_0": ceil(log2(base_0)),
                                "log_B_1": ceil(log2(base_1)),
                                "log_B_2": ceil(log2(base_2)),
                                "digits_0": digits_0,
                                "digits_1": digits_1,
                                "digits_2": digits_2,
                            },
                            "stats": stats(base_0, base_1, base_2),
                        })
  return valid_outputs



if __name__ == "__main__":
  DB_size_GB = 1
  # for payload_bits_input in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
  for payload_bits_input in [1]:
    N_input = ceil(DB_size_GB * 1024 * 1024 * 1024 * 8 / payload_bits_input)
    # costs_medium_payload = compute_costs_of_double_ypir_medium_payload(
    #     N=N_input,
    #     payload_bits_input=payload_bits_input,
    # )
    # costs_small_payload = compute_costs_of_double_ypir_small_payload(
    #     N=N_input,
    #     payload_bits_input=payload_bits_input,
    # )
    # costs_simple_ypir = compute_costs_of_simple_ypir(
    #     N=N_input,
    #     payload_bits_input=payload_bits_input
    # )
    # costs = costs_medium_payload + costs_small_payload + costs_simple_ypir
    costs = compute_costs_of_inspire(
        N_input=N_input,
        payload_bits_input=payload_bits_input
    )
    costs = sorted(costs, key=lambda x: x["stats"]["comm"]["total_size"])
    print(
        costs[0]["params"]["payload_bits_input"],
        costs[0]["stats"]["protocol_type"],
        costs[0]["stats"]["comm"]["total_size_KB"]
    )
    # print(costs[0])
  
  # dump costs to file, one item per line
  import json

  with open("costs.txt", "w") as f:
    for i, cost in enumerate(costs):
      if i == 1000:
        break
      f.write(json.dumps(cost))
      f.write("\n")

  # print stats with lowest total size
  print(costs[0])

  # print("Generating plots...")

  # # for N=2**20, 2**21, ..., 2^30, compute the best cost and plot it
  # log_N_s = list(range(20, 31))
  # for payload_bytes_input in [32, 64, 128, 256]:

  #   best_cost = []
  #   for log_N in log_N_s:
  #     costs = compute_costs_of_double_ypir_medium_payload(
  #         N=2**log_N, payload_bits_input=payload_bytes_input * 8
  #     )
  #     costs = sorted(costs, key=lambda x: x["stats"]["comm"]["total_size"])
  #     best_cost.append(costs[0]["stats"]["comm"]["total_size_KB"])

  #   plt.plot(log_N_s, best_cost, label=f"size={payload_bytes_input} B")
  # plt.xlabel("log_N")
  # plt.ylabel("total size (KB)")
  # plt.legend()
  # plt.savefig("total_size_vs_N.png")
