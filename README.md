# FootTrack

Footprint Analysis for Tracking TF Occupancy and Kinetics

based on TOBIAS ([https://github.com/loosolab/TOBIAS](https://github.com/loosolab/TOBIAS))

Introduction 
------------

Gene regulation is orchestrated by the precise binding of specific transcription factors (TFs) at regulatory elements. However, simultaneously detecting the binding of hundreds of TFs on chromatin remains challenging. We have developed a cytosine deaminase-based footprinting assay for occupancy and organization of transcription factors by sequencing (cFOOT-seq), enabling high-resolution, sensitive, and quantitative genome-wide assessment of TF binding for hundreds of TFs simultaneously.

To comprehensively assess transcription factor occupancy and dynamics under various conditions through TF footprint analysis, we developed a computational framework named FootTrack (Footprint Analysis for Tracking TF Occupancy and Kinetics), building on the TOBIAS. Using known TF binding information and motif information, FootTrack precisely maps TF occupancy and the chromatin landscape around TF motif centers at known binding sites. By integrating cFOOT-seq data and motif information from JASPAR, FootTrack also facilitates de novo prediction of transcription factor binding sites genome-wide.

<img src="/figures/framework.png">

Installation
------------

To install directly from the repository, the python package Cython is needed to build the C-extensions used by FootTrack. You can obtain Cython with `pip install cython`. FootTrack can then be installed using the setup.py script:

```bash
cd FootTrack
conda env create -f foottrack_env.yaml
conda activate FootTrack_ENV
python setup.py install
```

License
------------
This project is licensed under the [MIT license](LICENSE). 
