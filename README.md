# KnapsackSearch
 Automated data search in the KNApSAcK database

## Aim
 [KNApSAcK](http://www.knapsackfamily.com) is a highly useful source of information about natural products, which is accessible through a web interface.
 
 The goal of the KnapsackSearch python scripts is to create an SDF structure file with molecules related to one or more genera of living organisms, typically those from the same botanical family.
 
 The user of KnapsackSearch first creates a file named *familyname*_genera.txt, so that *familyname* stands for the nickname of a botanical family such as 'papaver' for the family of Papaveracea.
 
 The file *familyname*_genera.txt may contain lines of three kinds:
 
 1. Blank lines, considered as a comment
 2. Lines that start with a # sign, considered as a comment
 3. Lines with a single word, standing for a *genus* name
 
 All associations between all organism names and compounds for all given genera names are searched in KNApSAcK through the web and a list of compounds is established.
 
 Data are then extracted from KNApSAcK for each compound, such as elemental formula, atomic, mass, one of more names, SMILES and InChI character strings.
 
 The SMILES character strings are processed by [RDKit](https://www.rdkit.org/) to create a single SDF file that also contains 2D atomic coordinates.
 
 Each compound is associated to a list of carbon-13 NMR chemical shifts predicted by [nmrshiftdb2](https://nmrshiftdb.nmr.uni-koeln.de/) and presented as NMREDATA_ASSIGNMENT tags according to [NMReDATA](https://nmredata.org/).
 
 The whole KnapsackSearch process transforms the file *familyname*_genera.txt into the file *familyname*.nmredata.sdf. See file data_flow.txt and comments in the python scripts for the content of the eight intermediate files.
 
 The molecules built from SMILES strings are converted to InChI strings. The molecules for which this recalculated InChI string is different from the original one are discarded from the final result file. Discrepancies seem to arise from atom configuration differences.
 
## Installation

 Install rdkit from [Anaconda](https://www.anaconda.com/distribution/) as recommended in https://www.rdkit.org/docs/Install.html, if not already done.
 
 Install the python *requests* module from the rdkit environment (`conda install -c anaconda requests`)
 
 Install a [Java runtime environment](https://www.java.com/fr/download/) (jre) if not already done.
 
 Windows: Change the first part of the command in file predictSdf.bat, so that what is written between double quotes is the path to the java executable. It may be simply "java" if Java was installed by the default way. A non-default way can be to install first Java by the default way, copy the jre folder at the place you like and uninstall Java. This non-default Java installation may be desired to avoid being regularly warned by the messages sent by the Java updating service.
 
 Linux and Mac: Copy the files in the MacLinux folder to the main folder. Ensure that file predictSdf has execution permission.
 
 [Notepad++](https://notepad-plus-plus.org/downloads/) is recommended for text edition works.
 
 The predictSdf.bat file has to be edited in order to work with Linux/Mac. Windows uses ';' as separator between file name where Linux/Mac require ','.
 
 ## Usage
 
 Enter command `python process.py familyname` from the rdkit environment. As an example run `python process.py papaver` to collect data about compounds reported in KNApSAcK from Papaveracea, according to the list of genera written in file papaver_genera.txt. On November 19, 2019, the resulting papaver.nmredata.sdf file contained 476 molecules.

 The viewing of molecular structures and their attributes in SDF files is conveniently achieved by means of the [EdiSDF](http://infochim.u-strasbg.fr/spip.php?rubrique41) software.
 