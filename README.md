# NPBDetect: A Neural Network Model to Detect Bioactivities from AntiSMASH GBK Files

NPBdetect is a neural network model which can be trained and used to predict natural product bioactivity from biosynthetic gene clusters (BGCs).The tool requires output files from antiSMASH (version 7). NPBDetect is able to predict 8 different bioactivities like antibacterial, antifungal and antitumor or cytotoxic, siderophore, antiprotozoal, inhibitor, antiviral, and surfactant.

NPBDetect is available via:

<table>
<tr>
<td><a href="#cli-installation"><img src="https://edent.github.io/SuperTinyIcons/images/svg/powershell.svg" width="70" title="CLI"></a></td>
<td><a href="#google-colab"><img src="https://edent.github.io/SuperTinyIcons/images/svg/colaboratory.svg" width="70" title="Google Colab"></a></td>
<td><a href="#docker"><img src="https://edent.github.io/SuperTinyIcons/images/svg/docker.svg" width="70" title="Docker"></a></td>
</tr>
</table>

----

## Contents

- [Usage](#usage)
- [CLI Installation](#cli-installation)
- [Google Colab](#google-colab)
- [Docker](#docker)

----

## Usage

1. **Detect biosynthetic gene clusters (BGCs) with antiSMASH 7**

   The first step is to generate `.gbk` files using [antiSMASH 7](https://antismash.secondarymetabolites.org). These GBK files are used by NPBDetect to predict bioactivities.

   - **Online Usage**:
     You can use the antiSMASH web service here:  
     `https://antismash.secondarymetabolites.org`

   - **Local Installation**:
     - Download antiSMASH: `https://antismash.secondarymetabolites.org/#!/download`
     - Follow installation instructions: `http://docs.antismash.secondarymetabolites.org/install/`

   Example command to generate a `.gbk` file:

   ```bash
   antismash --output_dir <OUTPUT_FOLDER> --minlength 500 --cb-general --cb-knownclusters --cb-subclusters --fullhmmer --asf --pfam2go --smcog-trees -c 32 --genefinding-tool prodigal <INPUT_file>
   ```

2. **Use NPBDetect**

   NPBDetect can be run in two ways:
   - Online via Google Colab
   - Locally via CLI (Conda environment) or Docker

----

## CLI Installation

1. **Set up a Conda environment**

   ```bash
   conda create -n npbdetect python=3.10
   conda activate npbdetect
   ```

2. **Install required Python packages**

   ```bash
   pip install pandas
   pip install scikit-learn
   pip install biopython
   pip install torch torchaudio torchvision torchtext torchdata
   ```

3. **Clone the NPBDetect repository**

   ```bash
   git clone https://github.com/cbl-nabi/NPBDetect
   cd NPBDetect
   ```

4. **Validate installation and test sample prediction**

   ```bash
   python NPBDetect.py predict \
       -v 1 \
       --gbk test/BGC0000004.gbk \
       --pred HC \
       --out_dir outs/
   ```

----

## Google Colab

NPBDetect is also available online through Google Colaboratory.  
No installation required!

👉 **Launch it here**: [NPBDetect@Google-Colab](https://colab.research.google.com/drive/1cVv9VgYTEHSPsBX-xVdVVNJgyjkXOz9E?usp=sharing#scrollTo=fcDLr4xs3jF2)

----

## Docker

You can also run NPBDetect easily using Docker. The image can be pulled 

```bash
docker pull mantrilabnabi/npbdetect
```

Example command:

```bash
docker run --rm -it \
    --volume <INPUT_DIR>:/data/input \
    --volume <OUTPUT_DIR>:/data/output \
    mantrilabnabi/npbdetect \
    python NPBDetect.py predict \
        -v 1 \
        --gbk /data/input/BGC0000004.gbk \
        --pred HC \
        --out_dir /data/output/
```

Replace `<INPUT_DIR>` and `<OUTPUT_DIR>` with your local paths.

Note: change output_type = from "HC" to "ORG" for 8 bioactivity class predictions.

----

## Output
The prediction/output file will be generated as a csv file containing the probability and prediction values for all bioactivity classes.

If the probability value ≥ 0.5 → prediction = 1

If the probability value < 0.5 → prediction = 0.

----


## License
The software is licensed under MIT License.


----


## Citation

If you use NPBDetect in your research, please cite:

> Hemant Goyat, Dalwinder Singh, Sunaina Paliyal, Shrikant Mantri,  
> [**Predicting biological activity from biosynthetic gene clusters using neural networks**](https://www.biorxiv.org/content/10.1101/2024.06.20.599829v1.full.pdf), 2024.

----
