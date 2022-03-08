# README: motif-mark 

# About
Motif-mark searches for and maps the positions of transcription factor binding motifs along gene sequences in input FASTA files. The Python script in this repo, `motif-mark-oop.py`, takes in at least one FASTA file (including those where each sequence is split across multiple lines rather than all being on one line) and at least one text file of motifs to be searched for, where each motif is on a new line. 

For each input FASTA file, the script outputs a one-line sequence version of the FASTA file, along with a PNG image file depicting relative positions of motifs along each gene sequence in the FASTA file. Each output PNG image depicts genes going left to right from 5'->3' and also denotes intron and exon positions. Moreover, each gene in an output PNG image is preceded by a text label that includes the gene name, the chromosome on which the gene is located, and the start and stop positions of the gene along that chromosome.


# How to Use

## 1. Configure Conda environment

For best results, make sure you are using __Python 3.9 or higher__. Create a Conda environment called `motif_mark` and include the following packages and versions:

```
# Name                    Version                   Build  Channel
appnope                   0.1.2            py39h6e9494a_2    conda-forge
asttokens                 2.0.5              pyhd8ed1ab_0    conda-forge
backcall                  0.2.0              pyh9f0ad1d_0    conda-forge
backports                 1.0                        py_2    conda-forge
backports.functools_lru_cache 1.6.4              pyhd8ed1ab_0    conda-forge
black                     22.1.0             pyhd8ed1ab_0    conda-forge
blas                      1.0                         mkl  
bottleneck                1.3.2            py39he3068b8_1  
brotli                    1.0.9                hb1e8313_2  
ca-certificates           2022.2.1             hecd8cb5_0  
cairo                     1.16.0               h8023c5d_1  
certifi                   2021.10.8        py39hecd8cb5_2  
click                     8.0.3            py39h6e9494a_1    conda-forge
cycler                    0.11.0             pyhd3eb1b0_0  
dataclasses               0.8                pyhc8e2a94_3    conda-forge
debugpy                   1.5.1            py39h9fcab8e_0    conda-forge
decorator                 5.1.1              pyhd8ed1ab_0    conda-forge
entrypoints               0.4                pyhd8ed1ab_0    conda-forge
executing                 0.8.2              pyhd8ed1ab_0    conda-forge
fontconfig                2.13.1               ha9ee91d_0  
fonttools                 4.25.0             pyhd3eb1b0_0  
freetype                  2.11.0               hd8bbffd_0  
gettext                   0.21.0               h7535e17_0  
giflib                    5.2.1                haf1e3a3_0  
glib                      2.69.1               h8346a28_1  
icu                       58.2                 h0a44026_3  
intel-openmp              2021.4.0          hecd8cb5_3538  
ipykernel                 6.9.0            py39h71a6800_0    conda-forge
ipython                   8.0.1            py39h6e9494a_0    conda-forge
jedi                      0.18.1           py39h6e9494a_0    conda-forge
jpeg                      9d                   h9ed2024_0  
jupyter_client            7.1.2              pyhd8ed1ab_0    conda-forge
jupyter_core              4.9.1            py39h6e9494a_1    conda-forge
kiwisolver                1.3.2            py39he9d5cce_0  
lcms2                     2.12                 hf1fd2bf_0  
libcxx                    12.0.0               h2f01273_0  
libffi                    3.3                  hb1e8313_2  
libgfortran               3.0.1                h93005f0_2  
libiconv                  1.16                 h1de35cc_0  
libpng                    1.6.37               ha441bb4_0  
libsodium                 1.0.18               hbcb3906_1    conda-forge
libtiff                   4.2.0                h87d7836_0  
libwebp                   1.2.2                h56c3ce4_0  
libwebp-base              1.2.2                hca72f7f_0  
libxml2                   2.9.12               hcdb78fc_0  
llvm-openmp               12.0.0               h0dcd299_1  
lz4-c                     1.9.3                h23ab428_1  
matplotlib                3.5.1            py39hecd8cb5_0  
matplotlib-base           3.5.1            py39hfb0c5b7_0  
matplotlib-inline         0.1.3              pyhd8ed1ab_0    conda-forge
mkl                       2021.4.0           hecd8cb5_637  
mkl-service               2.4.0            py39h9ed2024_0  
mkl_fft                   1.3.1            py39h4ab4a9b_0  
mkl_random                1.2.2            py39hb2f4e1b_0  
munkres                   1.1.4                      py_0  
mypy_extensions           0.4.3            py39h6e9494a_4    conda-forge
ncurses                   6.3                  hca72f7f_2  
nest-asyncio              1.5.4              pyhd8ed1ab_0    conda-forge
numexpr                   2.8.1            py39h2e5f0a9_0  
numpy                     1.21.2           py39h4b4dc7a_0  
numpy-base                1.21.2           py39he0bd621_0  
olefile                   0.46               pyhd3eb1b0_0  
openssl                   1.1.1m               hca72f7f_0  
packaging                 21.3               pyhd3eb1b0_0  
pandas                    1.4.1            py39he9d5cce_0  
parso                     0.8.3              pyhd8ed1ab_0    conda-forge
pathspec                  0.9.0              pyhd8ed1ab_0    conda-forge
pcre                      8.45                 h23ab428_0  
pexpect                   4.8.0              pyh9f0ad1d_2    conda-forge
pickleshare               0.7.5                   py_1003    conda-forge
pillow                    8.4.0            py39h98e4679_0  
pip                       21.2.4           py39hecd8cb5_0  
pixman                    0.40.0               h9ed2024_1  
platformdirs              2.5.0              pyhd8ed1ab_0    conda-forge
prompt-toolkit            3.0.27             pyha770c72_0    conda-forge
ptyprocess                0.7.0              pyhd3deb0d_0    conda-forge
pure_eval                 0.2.2              pyhd8ed1ab_0    conda-forge
pycairo                   1.19.1           py39h06c6e95_0  
pygments                  2.11.2             pyhd8ed1ab_0    conda-forge
pyparsing                 3.0.4              pyhd3eb1b0_0  
python                    3.9.7                h88f2d9e_1  
python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
python_abi                3.9                      2_cp39    conda-forge
pytz                      2021.3             pyhd3eb1b0_0  
pyzmq                     22.3.0           py39h7fec2f1_1    conda-forge
readline                  8.1.2                hca72f7f_1  
scipy                     1.7.3            py39h8c7af03_0  
seaborn                   0.11.2             pyhd3eb1b0_0  
setuptools                58.0.4           py39hecd8cb5_0  
six                       1.16.0             pyh6c4a22f_0    conda-forge
sqlite                    3.37.2               h707629a_0  
stack_data                0.2.0              pyhd8ed1ab_0    conda-forge
tk                        8.6.11               h7bc2e8c_0  
tomli                     2.0.1              pyhd8ed1ab_0    conda-forge
tornado                   6.1              py39h89e85a6_2    conda-forge
traitlets                 5.1.1              pyhd8ed1ab_0    conda-forge
typed-ast                 1.5.2            py39h89e85a6_0    conda-forge
typing_extensions         4.1.1              pyha770c72_0    conda-forge
tzdata                    2021e                hda174b7_0  
wcwidth                   0.2.5              pyh9f0ad1d_2    conda-forge
wheel                     0.37.1             pyhd3eb1b0_0  
xz                        5.2.5                h1de35cc_0  
zeromq                    4.3.4                he49afe7_1    conda-forge
zlib                      1.2.11               h4dc903c_4  
zstd                      1.4.9                h322a384_0
```

## 2. Format input files


## 3. Run the Python script


