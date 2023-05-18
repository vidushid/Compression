# Compression
Compression of assembled human genomes with gene query


Download data from the Human Pangenome Project (47 samples), GRCh 38 reference, and CHM13 v.1.1 assembly: <br />
<pre><code>
curl -o HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1 
</pre> </code>
Retrieve data using AGC: <br />
<pre><code>
git clone https://github.com/refresh-bio/agc <br />
cd agc && make <br />
./agc getcol -o ./decompressed HPRC-yr1.agc <br />
</pre> </code>
