# KnapsackSearch
 Automated data search in the KNApSAcK database
 
## KNApSAcK Change
February 17, 2020. KNApSAcK changed and KnapsackSearch had to change accordingly.

See family_from_web.py and compounds.py for details.

## Aim
[KNApSAcK](http://www.knapsackfamily.com) is a highly useful source of information about natural products, which is accessible through a web interface.

The goal of the KnapsackSearch python scripts is to create an SDF structure file with molecules related to one or more genera of living organisms, typically those from the same botanical family.

All associations between all organism names and compounds for all given genera names are searched in KNApSAcK through the web and a list of compounds is established.

Data are then extracted from KNApSAcK for each compound, such as elemental formula, atomic, mass, one of more names, SMILES and InChI character strings.

The SMILES character strings are processed by [RDKit](https://www.rdkit.org/) to create a single SDF file that also contains 2D atomic coordinates.

Each compound is associated to a list of carbon-13 NMR chemical shifts predicted
by [nmrshiftdb2](https://nmrshiftdb.nmr.uni-koeln.de/) under the NMRSHIFTDB2_ASSIGNMENT SDF tag.

The whole KnapsackSearch process transforms the file *familyname*_genera.txt into the file *familyname*_knapsack.sdf. 
See file data_flow.txt and comments in the python scripts for the content of the eight intermediate files.

The molecules built from SMILES strings are converted to InChI strings.
The molecules for which this recalculated InChI string is different from
the original one are discarded from the final result file. Discrepancies seem to arise from atom configuration differences.

## Installation

Install rdkit from [Anaconda](https://www.anaconda.com/distribution/) as
recommended in https://www.rdkit.org/docs/Install.html, if not already done.

Install the python *requests* module from the rdkit environment (`conda install requests`)

Install a [Java runtime environment](https://www.java.com/fr/download/) (jre) if not already done.

See below section *MyInstallation* for updated installation instructions.

Linux and Mac: Replace process.py by MacLinux/process.py and predictSdf.bat by MacLinux/predictSdf.
Ensure that the script file predictSdf has execution permission.

[Notepad++](https://notepad-plus-plus.org/downloads/) is the recommended text editor for Windows.

## Usage
 
The user of KnapsackSearch first creates a file named *familyname*_genera.txt, so that *familyname* stands for the nickname of a botanical family such as 'papaver' for the family of Papaveracea.

The file *familyname*_genera.txt may contain lines of three kinds:

1. Blank lines, considered as a comment
2. Lines that start with a # sign, considered as a comment
3. Lines with a single word, standing for a *genus* name, starting with an upper-case letter (A-Z)

Enter command `python -m process familyname` from the rdkit environment.
As an example run `python -m process papaver` to collect data about compounds reported in KNApSAcK from Papaveraceae,
according to the list of genera written in file papaver_genera.txt.
On April 9, 2020, the resulting papaver_knapsack.sdf file contained 458 molecules.
Other examples can be found in the Examples directory.

The list of genera that belong to a given family can be found by means of the [NCBI Taxonomy tool](https://www.ncbi.nlm.nih.gov/taxonomy).

The viewing of molecular structures and their attributes in SDF files is conveniently achieved by means of the [EdiSDF](http://infochim.u-strasbg.fr/spip.php?rubrique41) software.

## Limitations
  
A compound name may be assigned to two different compounds
  
  
## Acknowledgments 

nmrshiftdb2 implementation would not have been possible without the
inspiration and help from Pr. Christoph Steinbeck (University of Jena, Germany)
and Dr. Stefan Kuhn (De Monfort University, Leicester, UK).

# Fake_ACD

## Aim

Fake_ACD creates SDF files with 13C NMR chemical shifts predicted by nmrshiftdb2 and formatted
to be read by ACD software to produce databases with molecular structures and 13C NMR data.

## Test

Any recent python interpreter should work. No need for a particular environment.

The Fake_ACD_Results directory contains the files produced by the following tests.

### sdfrw.py

See [SDFrw](https://github.com/nuzillard/SDFrw) for sdfrw.py. 

`python -m sdfrw quercetin2D.sdf`

creates `copied_quercetin2D.sdf`, a copy of `quercetin2D.sdf`.

This test can be skipped and is only there to ensure that the following ones will have a chance to succeed.

### addnmrsdb.py

`python -m addnmrsdb quercetin2D.sdf`

creates `nmrsdb_quercetin2D.sdf`, a copy of `quercetin2D.sdf` with
added chemical shifts values from nmrshiftdb2.

Input files to addnmrsdb.py are .sdf files with 2D coordinates (z atom coordinates are 0)
with possible information about configuration of stereocenters.

Output files from addnmrsdb.py are .sdf files with an added NMRSHIFTDB2_ASSIGNMENT tag and value
lines, one per carbon atom, formatted like `8, 105.34 \`, stating that carbon 8 (indexing starts at 1)
has an nmrshiftdb2-predicted chemical shift value of 105.34 ppm.

### fakeACD.py

`python -m fakeACD nmrsdb_quercetin2D.sdf`

creates `fake_acd_nmrsdb_quercetin2D.sdf`, a copy of `nmrsdb_quercetin2D.sdf` with
original chemical shifts values from nmrshiftdb2 formatted in the style of `addnmrsdb.py`
and reformatted under the CNMR_SHIFTS tag.

File `fake_acd_nmrsdb_quercetin2D.sdf` can be imported to an ACD DB to produce
database file `fake_acd_nmrsdb_quercetin2D.NMRUDB`.

fakeACD.py can process output files from KnapsackSearch.

### fakefakeACD.py

`python -m fakefakeACD quercetin2D.sdf`

creates `fake_acd_quercetin2D.sdf`, a copy of `quercetin2D.sdf` with
very fake (99.99) chemical shifts values in the style of `addnmrsdb.py`
and reformatted under the CNMR_SHIFTS tag.

File `fake_acd_quercetin2D.sdf` can be imported to an ACD DB to produce
database file `fake_acd_quercetin2D.NMRUDB`.

### rdcharge.py

`python -m rdcharge filename.sdf filename_elec.sdf`

creates `filename_elec.sdf` from `filename.sdf`. This is necessary when .sdf files from RDKit are produced to be 
read by ACD software. `rdcharge.py` corrects the description of electrically charged atom, for which RDKit
issues a non-zero valence information field that is not correctly interpreted by ACD and inhibits
the prediction of chemical shift values. `rdcharge.py` is included in the script `process.py` of KnapsackSearch.

`python -m rdcharge filename.sdf`

applies an in-place correction.

### uniqInChI.py

`python -m uniqInChI filename.sdf filename_uniq.sdf`

creates `filename_uniq.sdf` from `filename.sdf`. This is necessary when .sdf files
contain duplicate compounds, according to the corresponding InChI.

`python -m uniqInChI filename.sdf`

applies an in-place elimination of duplicates.

### tautomer.py

`python -m tautomer filename.sdf filename_tauto.sdf`

creates `filename_tauto.sdf` from `filename.sdf` in which tautomer correction was applied.
Among others aliphatic iminols are converted to amides.

`python -m tautomer filename.sdf`

applies an in-place conversion of tautomers.

# Quick ACD DB with calculated experimental 13C NMR data

The use of the related files, in directory `CNMR_Predict` requires the availability of ACD/Labs CNMR Predictor and DB.

## Aim

Transformation of a .smi file containing SMILES chains and compound names separated by a single space (.smi file with 1 line per compound)
into a database (ACD DB) with "experimental" 13C NMR chemical shifts determined by ACD prediction.

The python script `smi2ACD.py` processes a .smi file to produre a minimal .sdf file
in which the structures can be imported in an ACD Database.
The python script `CNMR_predict.py` transforms a .sdf file with calculated chemical shift values from ACD/Lasbs DB
into another .sdf file in which the calculated values replace the supposedly experimental ones.
See section for Example 1, hereafter.

The python scripts `smi2ACD.py` and `CNMR_predict.py` may be used independently for other purposes.
Note that `smi2ACD.py` assigns 99.99 as a placeholder for the experimental chemical shift value of all carbon atoms.
More realistic values, from nmrshiftdb2, may be obtained by action of 'addnmrsdb.py'
on a .sdf file.

The combination of KnapsackSearch (`process.py`) and `fakeACD.py` produces .sdf files
that can be processed by `CNMR_predict.py`. The resulting files contain
13C NMR chemical shifts calculated by nmrshiftdb2 and by ACD. See section for Example 2, hereafter.

`CNMR_predict.py` produces a .sdf compound library file with tags compatible with its use by the
[MixONat](https://sourceforge.net/projects/mixonat/) software, described in
[this publication](https://dx.doi.org/10.1021/acs.analchem.0c00193).
The action of `ACD_to_DerepCrude.py` on a .sdf file with 13C chemical shifts
formatted for ACD produces a .sdf file suitable with a use by the
[DerepCrude](http://eos.univ-reims.fr/LSD/JmnSoft/DerepCrude/) software, described in
[this publication](https://doi.org/10.1021/acs.jnatprod.6b01063).


## Example 1

The example file `small.smi` contains 2 lines, one for quercetin and the other one for resveratrol.
Their SMILES chains were copied from Wikipedia.

See `tutorial_CNMR_Predict.pdf` for detailed explanations.


Running the example requires an RDKit environment.

1. `python -m smi2ACD small.smi fake_acd_small.sdf`
2. Create DB `fake_acd_small.NMRUDB` and import `fake_acd_small.sdf`
3. Calculate 13C NMR chemical shifts in ACD/Labs DB: Database->Tools->Check Chemical Shifts
4. Export DB as `fake_acd_small_exported.sdf`
5. `python -m CNMR_predict fake_acd_small_exported.sdf calc_acd_small.sdf` copies calculated data as if they were experimental.
6. Create DB `calc_acd_small.NMRUDB` and import `calc_acd_small.sdf`
7. Calculate again 13C NMR chemical shifts: Database->Tools->Check Chemical Shifts

The last step is optional but shows that the "experimental" chemical shift values
are the same as the calculated ones, obviously because they are calculated in the same way.

Files in directory `Small_results` were created from `small.smi` in the following order:

1. `fake_acd_small.sdf` (requires RDKit and `smi2ACD.py`, step 1)
2. `fake_acd_small.NMRUDB` (requires ACD software, steps 2 and 3)
3. `fake_acd_small_exported.sdf` (requires ACD software, step 4)
4. `calc_acd_small.sdf` (requires RDKit and `CNMR_predict.py`, step 5)
5. `calc_acd_small.NMRUDB` (requires ACD software, steps 6 and 7)

## Example 2

Starting for `papaver_knapsack.sdf` as obtained hereabove from `papaver_genera.txt`
an ACD database with 13C chemical shifts from ACD may be produced as follows:

1. `python -m fakeACD.py papaver_knapsack.sdf` creates `fake_acd_papaver_knapsack.sdf`
2. Create DB `fake_acd_papaver_ks.NMRUDB` and import `fake_acd_papaver_knapsack.sdf`
3. Calculate 13C NMR chemical shifts in ACD/Labs DB: Database->Tools->Check Chemical Shifts
4. Export DB as `fake_acd_papaver_ks_exported.sdf`
5. `python -m CNMR_predict.py fake_acd_papaver_ks_exported.sdf calc_acd_papaver_ks.sdf` copies calculated data as if they were experimental.
6. Create DB `calc_acd_papaver_ks.NMRUDB` and import `calc_acd_papaver_ks.sdf`
7. Calculate again 13C NMR chemical shifts: Database->Tools->Check Chemical Shifts

The last created sdf file, `calc_acd_papaver_ks.sdf` is stored in directory `Papaver_result`.

# MyInstallation

This directory reports the details of my computer installation relatively to python, nmrshiftdb and java.

# LOTUS

This directory is about LOTUS, a natural product knowledge base. LOTUS might well supersede KNApSAcK.

# ClassyFire

KNApSAcK provides molecular structures selected according to biological taxonomy.
These structures can be enriched with chemical taxonomy data
from ClassyFire.

[ClassyFire](http://classyfire.wishartlab.com/) provides a chemical classification of compounds from the description of their molecular stucture.
A Python Classyfire API can be downloaded from [here](https://github.com/JamesJeffryes/pyclassyfire).
From a RDKit envivronment run `python setup.py install` to install the *pyclassyfire* package.

The `classyfy.py` module is a very simple wrapper around `pyclassyfire`.
From directory Classyfire:

`python -m classyfy flavonoids.sdf` (call with 1 argument)

creates `flavonoids_classified.sdf` from `flavonoids.sdf` and

`python -m classyfy flavonoids.sdf myfile.sdf` (call with 2 arguments)

creates `myfile.sdf` from `flavonoids.sdf`.

Be patient. The resulting files look like SDF files but are syntactically incorrect and cannot be used directly.

`python -m importClassyf flavonoids.sdf flavonoids_classified.sdf flavonoids_classified_merged.sdf` (call with 3 arguments)

creates `flavonoids_classified_merged.sdf`, a copy of `flavonoids.sdf` to which classification data from `flavonoids_classified.sdf` is appended.

`python -m importClassyf flavonoids.sdf flavonoids_classified.sdf` (call with 2 arguments)

supplements `flavonoids.sdf` with classification data from `flavonoids_classified.sdf`. The original content of `flavonoids.sdf` is overwritten.

`classyfy.py` and `importClassyf.py` rely on SDFrw and not on RDKit for the handling of the SDF files.

See also [pybatchclassyfire](https://pypi.org/project/pybatchclassyfire/).

