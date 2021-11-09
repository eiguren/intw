# intw3.3

Hau intw-ren lehenengo bertsioa da QE 6.4 eta nC-PP gure formatua erabiliz.

Konpilatzeko intw3.3 direktorio nagusian,

```bash
mkdir build
cd build
QE_HOME=/path/to/QE W_HOME=/path/to/W cmake ..
make -j8
```

Intel-eko konpilatzailearekin lan egiteko ondo dago `export FC=ifort CC=icx CXX=icp` aldagaiak eskuragarri izatea.
