use std::io::{Read, Write};
use std::net::{TcpStream};
use clap::Parser;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use spiral_rs::{
    client::*,
    gadget::*,
};
use spiral_rs::poly::{PolyMatrix, PolyMatrixRaw, PolyMatrixNTT};
use inspire::{client::*, measurement::*, modulus_switch::*, packing::*, params::*, scheme::*};
use std::time::Instant;
// use std::path::Path;

use inspire::commons::*;

/// Run the YPIR scheme with the given parameters
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    // /// Number of items in the database
    // #[clap(long)]
    // num_items: Option<usize>,

    // /// Interpolation degree
    // #[clap(long)]
    // dim0: Option<usize>,

    // /// Interpolation degree
    // #[clap(long)]
    // interpolate_degree: Option<usize>,

    // #[clap(long)]
    // item_size_bits: Option<usize>,

    #[clap(long)]
    which_item: Option<usize>,
}

fn main() {

    let args = Args::parse();
    let Args {
        which_item
    } = args;

    let which_item = which_item.unwrap();

    // --- 1. Initial computation to generate a query ---
    // println!("⚙️  Generating a large query...");
    // let query = "a".repeat(100_000); // Approx 100 KB
        let mut measurement = Measurement::default();

        let num_items = 1 << 15;
        let dim0 = 2048;
        let item_size_bits = 1 * (16 * 2048);

        let (params, interpolate_degree, _) 
            = params_rgswpir_given_input_size_and_dim0(num_items, item_size_bits, dim0);

        let gamma = params.poly_len;
        let packing_type = PackingType::InspiRING;

        let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
        let db_cols = params.instances * params.poly_len;

        // RLWE reduced moduli
        let rlwe_q_prime_1 = params.get_q_prime_1();
        let rlwe_q_prime_2 = params.get_q_prime_2();

        let db_cols_prime = db_cols / gamma;

        // ================================================================
        // QUERY GENERATION PHASE
        // ================================================================
        // let mut queries = Vec::new();

        let mut client = Client::init(&params);

        // let target_idx: usize = rng.gen::<usize>() % (db_rows * db_cols);
        // let target_row = target_idx / db_cols;
        // let target_col = target_idx % db_cols;
        // // let target_sub_col = (target_col % (interpolate_degree * gamma)) / gamma;

        let per = num_items / db_rows;
        let target_row = which_item / per;
        let target_col = (which_item % per) * (item_size_bits / 16);

        let target_sub_col = (target_col % (interpolate_degree * gamma)) / gamma;
        // debug!("Target item: {} ({}, {} ({}))", target_idx, target_row, target_col, target_sub_col);

        let start = Instant::now();
        let sk_reg = client.get_sk_reg();

        let packing_params = PackParams::new_fast(&params, gamma);

        let mut packing_keys = match packing_type {
            PackingType::InspiRING => {
                if gamma <= params.poly_len / 2 {
                    PackingKeys::init(&packing_params, sk_reg, W_SEED)
                } else {
                    PackingKeys::init_full(&packing_params, sk_reg, W_SEED, V_SEED)
                }
            },
            PackingType::CDKS => {
                PackingKeys::init_cdks(&params, sk_reg, STATIC_SEED_2)
            },
            PackingType::NoPacking => {
                panic!("Shouldn't be here");
            }
        };

        let y_client = YClient::new(&mut client, &params);
        let mut ct_gsw_body = PolyMatrixNTT::zero(&params, 1, 2 * params.t_gsw);

        let bits_per = get_bits_per(&params, params.t_gsw);
        for j in 0..params.t_gsw {
            let mut sigma = PolyMatrixRaw::zero(&params, 1, 1);
            let exponent = (2 * params.poly_len * target_sub_col / interpolate_degree) % (2 * params.poly_len);
            sigma.get_poly_mut(0, 0)[exponent % params.poly_len] = if exponent < params.poly_len {
                1u64 << (bits_per * j)
            } else {
                params.modulus - (1u64 << (bits_per * j))
            };
            let sigma_ntt = sigma.ntt();
            let ct = y_client.client().encrypt_matrix_reg(
                &sigma_ntt,
                &mut ChaCha20Rng::from_entropy(),
                &mut ChaCha20Rng::from_seed(RGSW_SEEDS[2*j+1]),
            );
            ct_gsw_body.copy_into(&ct.submatrix(1, 0, 1, 1), 0, 2 * j + 1);
            
            let prod = &y_client.client().get_sk_reg().ntt() * &sigma_ntt;
            let ct = &y_client.client().encrypt_matrix_reg(
                &prod,
                &mut ChaCha20Rng::from_entropy(),
                &mut ChaCha20Rng::from_seed(RGSW_SEEDS[2*j]),
            );
            ct_gsw_body.copy_into(&ct.submatrix(1, 0, 1, 1), 0, 2 * j);
        }        

        let query_row = y_client.generate_query_over_prime(
            SEED_0,
            params.db_dim_1,
            packing_type,
            target_row,
        );

        assert_eq!(query_row.len(), (params.poly_len + 1) * db_rows);
        let query_row_last_row: &[u64] = &query_row[params.poly_len * db_rows..];
        
        assert_eq!(query_row_last_row.len(), db_rows);
        let packed_query_row = pack_query(&params, query_row_last_row);

        measurement.online.client_query_gen_time_ms = start.elapsed().as_millis() as usize;
        let serialized_query = serialize_everything(&params, &mut packing_keys, packed_query_row, ct_gsw_body);
    
    // println!("✅ Query generated ({} bytes).", serialized_query.len());
    // thread::sleep(Duration::from_secs(1));

    // --- 2. Connect to the server and send the query ---
    // println!("\n📞 Connecting to server at 127.0.0.1:8081...");
    match TcpStream::connect("127.0.0.1:8081") {
        Ok(mut stream) => {
            // println!("🤝 Connection successful! Sending data...");
            
            // --- New Protocol: Send length prefix first ---
            // Get the length of the query as a 64-bit unsigned integer.
            let query_len = serialized_query.len() as u64;
            // Convert the u64 length into an array of 8 bytes (big-endian).
            stream.write_all(&query_len.to_be_bytes()).expect("Failed to write length prefix");

            // Now, write the actual query data to the stream.
            // stream.write_all(serialized_query.as_bytes()).expect("Failed to write query data");
            stream.write_all(serialized_query.as_slice()).expect("Failed to write query data");
            // println!("📤 Data sent. Waiting for response...");

            // --- 3. Wait for and process the response ---
            // let mut response = String::new();
            // The server will now send a response and then close the connection,
            // so read_to_string will complete successfully.
            // stream.read_to_string(&mut response).expect("Failed to read response");

            let q_1_bits = (rlwe_q_prime_1 as f64).log2().ceil() as usize;
            let q_2_bits = (rlwe_q_prime_2 as f64).log2().ceil() as usize;
            let total_sz_bits = ((q_2_bits * params.poly_len + q_1_bits * params.poly_len + 63) / 64) * 64; // rounding up to multiple of 64
            let total_sz_bytes = (total_sz_bits + 7) / 8;

            let c = db_cols_prime / interpolate_degree;
            let mut received_data = vec![0u8; c*total_sz_bytes as usize];
            // Read exactly `len` bytes into our vector.
            if stream.read_exact(&mut received_data).is_err() {
                println!("Failed to read the full message from the client.");
                return;
            }
            // println!("\n📥 Response received from server!");

            // The prompt guarantees this division will be clean (no remainder).
            let chunk_size = received_data.len() / c;

            // Use chunks_exact to create an iterator over slices of size `chunk_size`.
            // Then map each slice to an owned Vec<u8> and collect the results.
            let sum_switched: Vec<Vec<u8>> = received_data
                .chunks_exact(chunk_size)
                .map(|chunk| chunk.to_vec())
                .collect();
            
            let mut results = Vec::new();
            for which_poly in 0..c {
                let sum = PolyMatrixRaw::recover_how_many(&params, rlwe_q_prime_1, rlwe_q_prime_2, gamma, sum_switched[which_poly].as_slice());
                results.push(sum);
            }

            let rgsw_ans = results.iter().flat_map(|ct| {
                decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), params.poly_len).as_slice()[..gamma].to_vec()
            }).collect::<Vec<_>>();
            
            let final_result = rgsw_ans.as_slice();

            // /////////////////////////////////////////////
            // // Additional code for testing
            // let filename = " put your file name here ";
            // let read_db = read_file_into_matrix(Path::new(filename), num_items, item_size_bits).unwrap();
            // assert_eq!(read_db[which_item].len(), final_result.len());
            // assert_eq!(read_db[which_item].iter().map(|&x| x as u64).collect::<Vec<_>>().as_slice(), final_result);
            // /////////////////////////////////////////////

            for i in 0..params.poly_len {
                let original_string = format!("{}{}", (final_result[i] >> 8) as u8 as char, final_result[i] as u8 as char);
                print!("{}", original_string);
            } println!();

            // --- 4. Final computation after receiving response ---
            // println!("⚙️  Performing post-response computation...");
            // thread::sleep(Duration::from_secs(1));
            // println!("📄 Server Response: \"{}\"", response.trim());
            // println!("✅ All operations complete.");

        }
        Err(e) => {
            println!("❌ Failed to connect: {}", e);
        }
    }
}
