from math import ceil
import json

def b_to_B(x):
    return ceil(x/8)

# TODO: Update the offline and online runtimes accordingly
for (db_rows, db_cols, payload_bits, time_ns, offline_time_s) in [
    (2**15, 2**15, 8, 748208904, 213),
    (2**16, 2**17, 8, 2030662346, 1500),
    (2**17, 2**18, 8, 5908738899, 5700),
]:
    k = 2
    l = 6
    n = 4096
    log_q = 90
    log_Q = 32
    log_p = 8

    payload_factor = payload_bits / log_p
    keys_bits = (k + l) * (n * log_q)
    query_bits = db_cols * (log_Q)
    upload_bits = keys_bits + query_bits

    download_bits = payload_factor * 2 * k * (db_rows + n) * log_q + payload_factor * db_rows * (log_Q)
    total_bits = keys_bits + query_bits + download_bits

    database_size_mb = (db_rows * db_cols * payload_bits) / (8 * 1024 * 1024)
    time_ms = time_ns // 1000000

    # round to 2 digits
    throughputGBs = database_size_mb / 1024 / (time_ms / 1000)
    throughputGBs = round(throughputGBs, 2)

    cost = lambda x : round(1000000 * (x[0] * 0.09 / (1024*1024) + (x[1]/1000) * 1.5 / 100000), 1)

    print("--------------------------------------------")
    print("Database Size (MB):", database_size_mb)
    # print("      Payload (KB):", db_rows * log_p / 8192)
    print("  Offline Time (s):", offline_time_s)
    print("         Keys (KB):", keys_bits / 8192)
    print("        Query (KB):", query_bits / 8192)
    print("     Download (KB):", download_bits / 8192)
    print("        Total (KB):", total_bits / 8192)
    print("  Server Time (ms):", time_ms)
    print(" Throughput (GB/s):", throughputGBs)
    print("              Cost:", cost(((total_bits+7)/8192, time_ms)))
    print("--------------------------------------------")


    output = {
        "specs": {
            # "label": "",
            "inputDatabaseSizeMb": database_size_mb,
            # "inputNumItems": 8589934592,
            # "inputItemSizeBits": 1,
            # "polyLen": 2048,
            # "modulusBits": 53,
            # "ptModulus": 65536,
            # "protocolType": "InsPIRe",
            # "secondLevelPackingMask": "InspiRING",
            # "secondLevelPackingBody": "InspiRING",
            # "gamma0": 16,
            # "gamma1": 1024,
            # "gamma2": 16,
            # "interpolateDegree": 0,
            # "resizedDbFirstDim": 8192,
            # "resizedDbSecondDim": 4096,
            # "resizedItemSizeBits": 256,
            # "resizedDatabaseSizeMb": 1024.0,
            # "onlineOnly": false
        },
        "offline": {
            # "uploadBytes": 0,
            # "downloadBytes": 0,
            # "encodeTimeMs": 0,
            "serverTimeMs": offline_time_s * 1000,
            # "clientTimeMs": 0,
            # "simplepirPrepTimeMs": 11007,
            # "simplepirHintBytes": 0,
            # "doublepirHintBytes": 0
        },
        "online": {
            "uploadBytes": b_to_B(upload_bits),
            "uploadKeys": b_to_B(keys_bits),
            "uploadQuery": b_to_B(query_bits),
            "downloadBytes": b_to_B(download_bits),
            "totalBytes": b_to_B(total_bits),
            # "simplepirQueryBytes": 54272,
            # "doublepirQueryBytes": 27136,
            # "simplepirRespBytes": 0,
            # "doublepirRespBytes": 0,
            "serverTimeMs": time_ms,
            # "clientQueryGenTimeMs": 327,
            # "clientDecodeTimeMs": 3,
            # "firstPassTimeUs": 103869,
            # "secondPassTimeUs": 3745,
            # "firstPackTimeUs": 549025,
            # "rgswTimeUs": 0,
            # "packingKeyRotationsTimeUs": 48120,
            # "packingMaskOnlineTimeUs": 28319,
            # "packingBodyOnlineTimeUs": 10178,
            # "totalRingPackingTimeMs": 0,
            # "sqrtNBytes": 0,
            # "allServerTimesMs": [
            # 752,
            # 728,
            # 743,
            # 757,
            # 767
            # ],
            # "stdDevServerTimeMs": 13.215142829345432,
            # "noiseWidth": 0.0
        }
    }
    with open(f'results/hintlesspir/{database_size_mb}.json', 'w') as f:
        json.dump(output, f)