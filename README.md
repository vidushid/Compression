# Compression
Compression of assembled human genomes with gene query


Download data from the Human Pangenome Project (47 samples), GRCh 38 reference, and CHM13 v.1.1 assembly: 
```bash
curl -o HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1 
```
Retrieve data using AGC: 
```bash
git clone https://github.com/refresh-bio/agc 
cd agc && make 
./agc getcol -o ./decompressed HPRC-yr1.agc  #to decompress 95 samples into the folder "decompressed"
```
To compile and run compression code use:
```bash
g++ -o cal calfor50.cpp -lzstd -L./zstd/lib/
./calfor50
```
