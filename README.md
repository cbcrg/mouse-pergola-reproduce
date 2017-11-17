# CB1-mice-Pergola-Reproduce.nf

This repository contains the software, scripts and data to reproduce the results corresponding to the CB1 mice experiment of the Pergola paper.

If you have not install yet [docker](https://www.docker.com/) and [nextflow](https://www.nextflow.io/), follow this [intructions](../README.md)

## Data processing

### Pull docker image
Pull the Docker image use for processing data with Pergola (Pergola and its dependencies installed)

```bash
docker pull pergola/pergola@sha256:8a52116be9bd371ae9bed9c0b36a8fc14634a7e14bbc764cc93905d8566e0939
```

### Clone the repository

```bash
git clone --recursive https://github.com/JoseEspinosa/pergola-paper-reproduce.git
cd pergola-paper-reproduce/cb1_mice
```

**Note**: If you have previously download the repository, then you only need to go to the ``cb1_mice`` folder.

### Data

Data is publicly available in [Zenodo](https://zenodo.org/) as a compressed tarball.

Data can be downloaded and uncompressed using the following command:

```bash
mkdir data
wget -O- https://zenodo.org/record/580312/files/CB1_mice.tar.gz | tar xz -C data
```

### Run nextflow pipeline
Once data is downloaded, it is possible to reproduce all the results using this command:

```bash
NXF_VER=0.24.3 nextflow run CB1_mice-Pergola-Reproduce.nf \ 
    --recordings='data/mice_recordings/' \
    --mappings='data/mappings/b2p.txt' \
    --mappings_bed='small_data/mappings/bed2pergola.txt' \
    --phases='small_data/mice_recordings/exp_phases.csv' \
    --exp_info='exp_info.txt' 
    --image_format='tiff'
    -with-docker
```

## Online visualization

### Downloading shiny-pergola config file

Download the configuration files assigning files to each of the group (wt_food_sc, wt_food_fat, cb1_food_sc, cb1_food_fat)

```bash
wget -O-  https://gist.githubusercontent.com/JoseEspinosa/9e65d54d765d9e5af554d837b3427569/raw/48fb424fb367c570461e7e6c8226abf81ead8ee2/cb1_pergola_conf.txt > exp_info.txt
```

### Downloading and running the shiny-pergola image

Pull the docker image containing the version of shiny-pergola web application used for render the data visualization:

```bash
docker pull pergola/shiny-pergola@sha256:e8791c5f230b612a6f702ac397849163e3a52b923befd1977e4a4c0235e91f72
```

With docker running, launch the image:

```bash
docker run --rm -p 3600:80 -v "$(pwd)":/pergola_data  pergola/shiny-pergola@sha256:e8791c5f230b612a6f702ac397849163e3a52b923befd1977e4a4c0235e91f72 &
```

**Note**: `"$(pwd)"` can be substitute by your absolute path to the folder where the `CB1_mice-Pergola-Reproduce.nf` has been run. 

**Note**: Figure has several snapshots if you want to get exactly the exact the same figure just select it by setting on id.txt file. For instance if you want to reproduce figure **a** you just have to type the following command before running Docker shiny-pergola image.

```bash
echo "cb1_a" > id.txt
```

Go to your web browser and type in your address bar the ip address returned by the following command e.g. http://0.0.0.0:3600

**Note**: In newer Docker versions by default the IP address used is the localhost ``0.0.0.0``. As this IP might be used by other services in the default port (80), we especified the port used in our docker command by ``-p 3600:80``. Old docker versions running in OS might used a different IP address that you can get using the following command (usually this one ``http://192.168.99.100``). You may then enter this IP followed by the same port in the address bar of your browser ``http://192.168.99.100:3600``.

```bash
docker-machine ip default 
```

Et voila, the running container will load the shiny app in your browser.

<img src="/cb1_mice/images/cb1_snapshot_shiny_pergola.png" alt="snapshot shiny-pergola" style="width: 100%;"/>
