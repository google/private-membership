use inspire::{packing::ToStr, scheme::run_ypir_batched};

use clap::Parser;
use inspire::packing::PackingType;
use inspire::scheme::ProtocolType;

/// Run the YPIR scheme with the given parameters
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Number of items in the database
    #[clap(long)]
    num_items: Option<usize>,
    /// Size of each item in bits (optional, default 1), values over 8 are
    /// unsupported
    #[clap(long)]
    item_size_bits: Option<usize>,
    /// Number of clients (optional, default 1)
    /// to perform cross-client batching over
    // #[clap(long)]
    // num_clients: Option<usize>,
    /// Number of trials (optional, default 5)
    /// to run the YPIR scheme and average performance measurements over (a
    /// warmup trial is excluded)
    #[clap(long)]
    trials: Option<usize>,
    /// Verbose mode (optional)
    /// if set, run as SimplePIR
    // #[clap(long, short, action)]
    // is_simplepir: bool,
    /// Protocol type (optional)
    /// if set, run as SimplePIR
    #[clap(long, value_parser)]
    protocol_type: Option<ProtocolType>,

    // #[clap(long, short, action)]
    // special_packing: bool,
    /// Second-level packing option (Option1, Option2, Option3)
    #[clap(long, value_parser)]
    second_level_packing_mask: Option<PackingType>,

    #[clap(long, value_parser)]
    second_level_packing_body: Option<PackingType>,

    #[clap(long)]
    performance_factor: Option<usize>,

    /// Small params mode (optional)
    /// if set, small paramseters will be used
    #[clap(long, short, action)]
    small_params: bool,

    /// Output report file (optional)
    /// where results will be written in JSON.
    #[clap(long)]
    out_report_json: Option<String>,
    /// Verbose mode (optional)
    /// if set, the program will print debug logs to stderr.
    #[clap(long, short, action)]
    verbose: bool,

    #[clap(long, short, action)]
    online_only: bool,

    #[clap(
        short,
        long,
        value_name = "gammas",
        help = "A comma-separated list of positive integers"
    )]
    gammas: Option<String>, // Store as String

}

fn main() {
    let args = Args::parse();
    let Args {
        num_items,
        item_size_bits,
        trials,
        out_report_json,
        verbose,
        protocol_type,
        second_level_packing_mask,
        second_level_packing_body,
        performance_factor,
        small_params,
        online_only,
        gammas,
    } = args;

    if verbose {
        println!("Running in verbose mode.");
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Debug)
            .write_style(env_logger::WriteStyle::Always)
            .init();
    } else {
        env_logger::init();
    }

    let num_items = num_items.unwrap_or(67108864);
    let item_size_bits = item_size_bits.unwrap_or(1);
    let gammas = gammas.unwrap_or("16,1024,16".to_string());
    let protocol_type = protocol_type.unwrap_or(ProtocolType::InsPIRe);
    let trials = trials.unwrap_or(3);
    let second_level_packing_mask = second_level_packing_mask.unwrap_or(PackingType::InspiRING);
    let second_level_packing_body = second_level_packing_body.unwrap_or(PackingType::InspiRING);
    let performance_factor = performance_factor.unwrap_or(1);
    let small_params = small_params;


    println!(
        "Protocol={}({}-{}), DB={} KB, trials={}",
        protocol_type,
        second_level_packing_mask.to_str(),
        second_level_packing_body.to_str(),
        (num_items * item_size_bits) / 8192,
        trials
    );

    let gammas_usize_result: Result<Vec<usize>, String> = gammas
    .split(',')
    .map(|n| {
        n.parse::<u64>() // Parse as u64 first
            .map_err(|e| e.to_string())
            .and_then(|num| {
                usize::try_from(num).map_err(|_| format!("{} is too large for usize", num))
            })
    })
    .collect();

    let mut measurement = run_ypir_batched(
        num_items,
        item_size_bits,
        protocol_type,
        second_level_packing_mask,
        second_level_packing_body,
        performance_factor,
        small_params,
        gammas_usize_result.unwrap(),
        trials,
        online_only
    );

    measurement.specs.input_num_items = num_items;
    measurement.specs.input_item_size_bits = item_size_bits;
    measurement.specs.input_database_size_mb = (num_items * item_size_bits) as f64 / (8. * 1024. * 1024.);


    if let Some(out_report_json) = out_report_json {
        println!("Writing report to {}", out_report_json);
        let mut file = std::fs::File::create(out_report_json).unwrap();
        serde_json::to_writer_pretty(&mut file, &measurement).unwrap();
        println!("Report written.");
    } else {
        println!(
            "Measurement completed. See the README for details on what the
            following fields mean."
        );
        println!("Result:");
        println!("{}", serde_json::to_string_pretty(&measurement).unwrap());    
    }
}
