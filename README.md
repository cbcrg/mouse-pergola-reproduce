# mouse-pergola-reproduce.nf
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1300200.svg)](https://doi.org/10.5281/zenodo.1300200)
![CircleCI status](https://circleci.com/gh/cbcrg/mouse-pergola-reproduce.png?style=shield)
[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.20.0-brightgreen.svg)](http://nextflow.io)

This repository contains the software, scripts and data to reproduce the results corresponding to the high-fat diet mice experiment of the Pergola paper.
The original source of this data set can be found on this [paper](https://onlinelibrary.wiley.com/doi/abs/10.1111/adb.12595)

If you have not install yet [docker](https://www.docker.com/) and [nextflow](https://www.nextflow.io/), follow this [intructions](https://github.com/cbcrg/pergola-reproduce/blob/master/README.md)

## Clone the repository

```bash
git clone --recursive https://github.com/cbcrg/mouse-pergola-reproduce.git
cd mouse-pergola-reproduce
```

## Data

Data is publicly available in [Zenodo](https://zenodo.org/) as a compressed tarball [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1154827.svg)](https://doi.org/10.5281/zenodo.1154827).

Data can be downloaded and uncompressed using the following command:

```bash
mkdir data
wget -O- https://zenodo.org/record/1300200/files/mouse_dataset.tar.gz | tar xz -C data
```

## Pull docker image
Pull the Docker image use for processing data with Pergola (Pergola and its dependencies installed)

```bash
docker pull pergola/pergola-reproduce@sha256:02bf3e701175104a488f40761b856efa1f97e2f2f82af8adae63b24ac2517326
```

## Run nextflow pipeline
Once data is downloaded, it is possible to reproduce all the results using this command:

```bash
NXF_VER=0.30.2 nextflow run mouse-pergola-reproduce.nf \
  --recordings='data/mouse_recordings/' \
  --mappings='data/mappings/b2p.txt' \
  --mappings_bed='data/mappings/bed2pergola.txt' \
  --phases='data/phases/exp_phases.csv' \
  --mappings_phase='data/mappings/f2g.txt' \
  --exp_info='data/mappings/exp_info.txt' \
  --tbl_chromHMM="data/chromHMM_files/cellmarkfiletable" \
  --n_bins_HMM=5 \
  --n_states_HMM=4 \
  --image_format='png' \
  -with-docker
```

##  Results

The previous command generates a results folder that contains the results of the analysis shown in the paper:

* A figure created using [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html) displaying the mouse feeding behavior and mouse accumulated food intakes within 30 minutes time-windows displayed as a heatmap.   
* A figure created using [Sushi](https://bioconductor.org/packages/release/bioc/html/Sushi.html) rendering the same layout as in the previous figure.
* A folder containing all the necessary files to render the raw feeding behavior and the the accumulated food intakes using [IGV](http://software.broadinstitute.org/software/igv/). Data is separated in folders corresponding to each mouse group.
* A heatmap comparing the feeding behavior of high-fat mice group with their controls (fold change). 
* A folder named **feeding_activity_profiles** containing the actograms generated using [deepTools](https://deeptools.readthedocs.io/en/develop/#).
* All files regarding the HMM modeling of the data using [chromHMM](http://compbio.mit.edu/ChromHMM/) can be found on **chromHMM** folder.
* A folder named **files** which contains the complete set of files needed to visualize the data using Shiny-pergola as explained below.

## Online Shiny-Pergola visualization

From the folder where you have run the pipeline, you can visualize the results using Shiny-Pergola introducing the commands described below.
 
#### shiny-pergola config file

To set the membership of each of the files inside the folder with the same name, a config file is needed.
The repository includes a file assigning files to each of the groups (*Control,  HF*)
You can also check the file [here](https://gist.githubusercontent.com/JoseEspinosa/52e7e948fbeb5c9265d567b556628d31/raw/494957605cf99f53484d6c4ed7f6456ec86cb431/hf_mice_conf.txt) 

#### Downloading and running the shiny-pergola image

Pull the docker image containing the version of shiny-pergola web application used for render the data visualization:

```bash
docker pull pergola/shiny-pergola@sha256:e4c470b3916ce6ce0298b231d3a8d18bf72535beb04251c13439351c90797431
```

With docker running, launch the image:

```bash
docker run --rm -p 3600:80 -v "$(pwd)":/pergola_data pergola/shiny-pergola@sha256:e4c470b3916ce6ce0298b231d3a8d18bf72535beb04251c13439351c90797431 &
```

**Note**: `"$(pwd)"` can be substitute by your absolute path to the folder where the `mouse-pergola-reproduce.nf` has been run. 

**Note**: Paper contain several figures related to mouse dataset analysis. If you want to get exactly one of the figure 
just select it by setting the figure on id.txt file. For instance if you want to reproduce **Supplementary figure 4a** 
you just have to type the following command before running Docker shiny-pergola image.

```bash
echo "hf_1" > id.txt
```

The codes for each figure are:

| Figure        | Code   |
| ------------- | ------ |
| Fig. 4a       | hf_1   |
| Fig. 4b       | hf_2   |
| Supp. Fig. S2 | hf_3   |

Go to your web browser and type in your address bar the following ip address: http://0.0.0.0:3600

**Note**: In newer Docker versions by default the IP address used is the localhost ``0.0.0.0``. As this IP might be used 
by other services in the default port (80), we especified the port used in our docker command by ``-p 3600:80``. 
Old docker versions might used a different IP address that you can get using the following command (usually this one ``http://192.168.99.100``). 
You may then enter this IP followed by the same port in the address bar of your browser e.g. ``http://192.168.99.100:3600``.

```bash
docker-machine ip default 
```

Et voila, the running container will load the shiny app in your browser.

<img src="/images/HF_snapshot_shiny_pergola.png" alt="snapshot shiny-pergola" style="width: 100%;"/>
