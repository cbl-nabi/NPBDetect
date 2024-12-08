# NPBDetect: A Neural network model to detect bioactivities from Antismash GBK files

NPBDetect is available in:
<table>
<tr>
<td><img src="https://edent.github.io/SuperTinyIcons/images/svg/powershell.svg" width="100" title="CLI"></td>
<td><img src="https://edent.github.io/SuperTinyIcons/images/svg/colaboratory.svg" width="100" title="Google Colab"></td>
</tr>
<tr>
<td><img src="https://jupyter.org/assets/homepage/main-logo.svg" width="100" title="Jupyter"></td>
<td><img src="https://edent.github.io/SuperTinyIcons/images/svg/docker.svg" width="100" title="Docker"></td>
</tr>
</table>

----

## Contents ##

* [CLI](#cli-installation)
* [Google Colab](#google-colab)
* [Jupyter Notebook](#jupyter-notebook)
* [Docker](#docker)
* [Usage](#usage)

----

## CLI installation ##

1. Setup conda environment

```conda create -n npbdetect python=3.10```

2. Install packages with pip

```
conda activate npbdetect
pip3 install pandas
pip3 install scikit-learn
pip3 install biopython
pip3 install torch torchaudio torchvision torchtext torchdata
```
3. Clone repository
   
```git clone git@```

5. Test samples to validate installation
   
```python NPBDetect.py predict --gbk test/BGC0000004.gbk --pred HC --out_dir outs/```

----

## Google Colab ##
Click on link to use google colab without setting up locally:

*[NPBDetect@google-colab](https://colab.research.google.com/drive/12zXjqk1DFX8Ouv0rYSAWSoaOheluBys9#scrollTo=MiUQ7HHR2rC8)*

----

## Jupyter Notebook ##

Use local jupyter notebook as


----

## Docker ##

For docker users


----

## Usage ##

1. Generate GBK files with antiSMASH 7
   Use following command to generate GBK file:
   
   ```antismash```
   
2. Use any version of NPBDetect to predict bioactivity

----

If you use the code, please cite the paper:

Hemant Goyat, Dalwinder Singh Sunaina Paliyal, Mohit Panghal, Shrikant Mantri , *[Predicting biological activity from biosynthetic gene clusters using neural networks](https://www.biorxiv.org/content/10.1101/2024.06.20.599829v1.full.pdf))*, 2024.
