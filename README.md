# mouse-pergola-reproduce.nf

![CircleCI status](https://circleci.com/gh/cbcrg/mouse-pergola-reproduce.png?style=shield)
[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.20.0-brightgreen.svg)](http://nextflow.io)

This repository contains the software, scripts and data to reproduce the results corresponding to the CB1 mice experiment of the Pergola paper.

If you have not install yet [docker](https://www.docker.com/) and [nextflow](https://www.nextflow.io/), follow this [intructions](../README.md)

## Clone the repository

```bash
git clone --recursive https://github.com/cbcrg/mouse-pergola-reproduce.git
cd mouse-pergola-reproduce
```

## Data

Data is publicly available in [Zenodo](https://zenodo.org/) as a compressed tarball.

Data can be downloaded and uncompressed using the following command:

```bash
mkdir data
wget -O- https://zenodo.org/record/580312/files/mouse_dataset.tar.gz | tar xz -C data
```

## Pull docker image
Pull the Docker image use for processing data with Pergola (Pergola and its dependencies installed)

```bash
docker pull pergola/pergola@sha256:f7208e45e761dc0cfd3e3915237eb1a96eead6dfa9c8f3a5b2414de9b8df3a3d
```

## Run nextflow pipeline
Once data is downloaded, it is possible to reproduce all the results using this command:

```bash
NXF_VER=0.26.1 nextflow run mouse-pergola-reproduce.nf \
  --recordings='data/mouse_recordings/' \
  --mappings='data/mappings/b2p.txt' \
  --mappings_bed='data/mappings/bed2pergola.txt' \
  --phases='data/phases/exp_phases.csv' \
  --mappings_phase='data/mappings/f2g.txt' \
  --exp_info='data/mappings/exp_info.txt' \
  --image_format='tiff' \
  -with-docker
```

## Online visualization

#### Downloading shiny-pergola config file

Download the configuration files assigning files to each of the group (wt_food_sc, wt_food_fat, cb1_food_sc, cb1_food_fat)

```bash
wget -O-  https://gist.githubusercontent.com/JoseEspinosa/9e65d54d765d9e5af554d837b3427569/raw/48fb424fb367c570461e7e6c8226abf81ead8ee2/cb1_pergola_conf.txt > exp_info.txt
```

#### Downloading and running the shiny-pergola image

Pull the docker image containing the version of shiny-pergola web application used for render the data visualization:

```bash
docker pull pergola/shiny-pergola@sha256:e8791c5f230b612a6f702ac397849163e3a52b923befd1977e4a4c0235e91f72
```

With docker running, launch the image:

```bash
docker run --rm -p 3600:80 -v "$(pwd)":/pergola_data  pergola/shiny-pergola@sha256:e8791c5f230b612a6f702ac397849163e3a52b923befd1977e4a4c0235e91f72 &
```

**Note**: `"$(pwd)"` can be substitute by your absolute path to the folder where the `mouse-pergola-reproduce.nf` has been run. 

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

<img src="/images/cb1_snapshot_shiny_pergola.png" alt="snapshot shiny-pergola" style="width: 100%;"/>
