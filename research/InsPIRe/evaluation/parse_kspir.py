from math import ceil, log2
import json

# The key is database size in MB
kspir_stats = {
    256: {
        "offline_ms": 3508,
        "online_us": 230891,
    },
    512: {
        "offline_ms": 6848,
        "online_us": 409777,
    },
    1024: {
        "offline_ms": 14143,
        "online_us": 777786,
    }, 
    2048: {
        "offline_ms": 28932,
        "online_us": 1469723,
    }, 
    4096: {
        "offline_ms": 58915,
        "online_us": 2951009,
    }, 
    8192: {
        "offline_ms": 118288,
        "online_us": 5905015,
    }, 
    16384: {
        "offline_ms": 236361,
        "online_us": 11515959,
    }, 
    32768: {
        "offline_ms": 445654,
        "online_us": 50571281,
    },
}

for db_mb in [
    1024,
    8192,
    32768,
]:
    N = 4096
    log_q = 56
    log_p = 16

    factor = 4

    r = db_mb / 16
    k_auto = 4
    k_ks = 4
    k_rgsw = 2

    N1 = 128
    N2 = N/2 / N1

    num_keys = int(ceil(log2(r/factor)))

    key_packing = num_keys * k_ks * log_q * N
    key_bsgs = N2 * k_ks * log_q * N
    key_rgsw = k_rgsw * 2 * log_q * N
    total_keys = key_packing + key_bsgs + key_rgsw

    sub_query = N/2 * log_q
    total_upload = total_keys + sub_query

    response = factor * (2 * N * log_q)

    total = total_upload + response

    payload = factor * (N * log_p)

    print("--------------------------------------------")
    print("Database Size (MB):", db_mb)
    # print("      Payload (KB):", payload / 8192)
    # print("    Keys BSGS (KB):", key_bsgs / 8192)
    # print(" Keys Packing (KB):", key_packing / 8192)
    # print("    Keys RGSW (KB):", key_rgsw / 8192)
    print("   Total Keys (KB):", total_keys / 8192)

    print("     SubQuery (KB):", sub_query / 8192)
    
    print("  TotalUpload (KB):", total_upload / 8192)
    print("     Download (KB):", response / 8192)
    print("        Total (KB):", total / 8192)
    # print("  Server Time (ms):", time_ms)
    # print("Throughput (GB/s):", throughputGBs)
    print("--------------------------------------------")

    total_upload_bytes = ceil(total_upload / 8)
    total_keys_bytes = ceil(total_keys / 8)
    sub_query_bytes = ceil(sub_query / 8)
    response_bytes = ceil(response / 8)
    total_bytes = ceil(total / 8)

    output = {
        "specs": {
            # "label": "",
            "inputDatabaseSizeMb": db_mb,
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
            "serverTimeMs": kspir_stats[db_mb]['offline_ms'],
            # "clientTimeMs": 0,
            # "simplepirPrepTimeMs": 11007,
            # "simplepirHintBytes": 0,
            # "doublepirHintBytes": 0
        },
        "online": {
            "uploadBytes": total_upload_bytes,
            "uploadKeys": total_keys_bytes,
            "uploadQuery": sub_query_bytes,
            "downloadBytes": response_bytes,
            "totalBytes": total_bytes,
            # "simplepirQueryBytes": 54272,
            # "doublepirQueryBytes": 27136,
            # "simplepirRespBytes": 0,
            # "doublepirRespBytes": 0,
            "serverTimeMs": round(kspir_stats[db_mb]['online_us'] / 1000, -1),
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
    with open(f'results/kspir/{db_mb}.json', 'w') as f:
        json.dump(output, f)