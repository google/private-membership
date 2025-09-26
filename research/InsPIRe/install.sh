sudo apt-get update
sudo apt-get install -y build-essential
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
cargo build --release

cd evaluation
sudo apt install -y python3-virtualenv
virtualenv .venv
. .venv/bin/activate
pip install -r requirements.txt

# building YPIR
git clone https://github.com/menonsamir/ypir.git
cd ypir
cargo run --release -- 1073741824
cd ..

# building hintless_pir
git clone git@github.com:RasoulAM/hintless_pir-clone.git
cd hintless_pir-clone
sudo apt install -y bazel-bootstrap
sudo curl -L https://github.com/bazelbuild/bazelisk/releases/latest/download/bazelisk-linux-amd64 -o /usr/local/bin/bazel
sudo chmod +x /usr/local/bin/bazel
export USE_BAZEL_VERSION=7.1.2
bazel test //...
bazel build //hintless_simplepir:hintless_simplepir_benchmarks
cd ..

# building KSPIR
git clone git@github.com:RasoulAM/kspir-clone.git
cd kspir-clone
git clone https://github.com/intel/hexl.git -b v1.2.5 --single-branch 

sudo apt-get install -y cmake
cd evaluation
cd kspir-clone
cd hexl
mkdir build
cd build
cmake ..
make
sudo make install
cd ../../
mkdir build
cd build
cmake ..
make
