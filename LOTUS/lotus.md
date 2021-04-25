LOTUS, <https://lotus.naturalproducts.net/>, is a web site about organic natural products
that includes data from the biological taxonomy and bibliographic references.

A simple search for a taxon returns a collection of molecular structures 
that can be downloaded as an SDF file in V3000 format.

This text explains how to build an ACD database with these structures,
supplemented by predicted ^13^C NMR chemical shift values with the help of a template file
that can be used as a reminder.

The reminder files are not intended to be run by
any software but are simply here to indicate what to do and in which order,
using cut-and-paste of character strings from these files to either a python interpreter
or to file dialogs of ACD software.

Template files contain three-star strings that have to be replaced by a
name of your choice, a plant family or genus for example.
We might want to work about the *Oryza sativa* species, namely the rice,
and all three-star strings in templateA.txt can be replaced *oryza*.

On April, 25th 2021, LOTUS collected 235 data records from the search
for *Oryza sativa*.
The resulting downloaded file, named `lotus_simple_search_result.sdf`, was moved
to a directory named `Oryza_sativa`.
A copy of `templateA.sdf` was placed in this directory
and all the occurences of the three-star strings replaced by *oryza*.

The result looks like:

~~~
import fakefakeACD, CNMR_predict

	base	oryza.NMRUDB
	import	lotus_simple_search_result.sdf
	export	oryza.sdf
	close

fakefakeACD.fakefakeACD('oryza.sdf')
		--> fake_acd_oryza.sdf

	base	fake_acd_oryza.NMRUDB
	import	fake_acd_oryza.sdf
		Check Chemical Shifts
	export	fake_acd_oryza_exported.sdf
	close
	
CNMR_predict.CNMR_predict('fake_acd_oryza_exported.sdf', 'true_acd_oryza.sdf')
		--> true_acd_oryza.sdf

	base	lotus_oryza.NMRUDB
	import	true_acd_oryza.sdf
		Check Chemical Shifts
	export	lotus_oryza.sdf
	close
~~~

This list of instructions translates in practice into:

- Open a python interpreter with `Oryza_sativa` as current directory and
with access to rdkit

- Cut-and-paste `import fakefakeACD, CNMR_predict` to the python interpreter and validate this command

- Launch ACD CNMR predictor (and DB).

- Create a new ACD DB named `oryza.NMRUDB`

- Import `lotus_simple_search_result.sdf` from directory `Oryza_sativa` in this DB.

	At this stage some filtering of the DB content may be carried out, depending on the purpose of the search.
Use the reminder file to document the changes you brought to the DB

- Export the DB as `oryza.sdf` and close this DB. Exportation in SDF V2000 format.

- Cut-and-paste `fakefakeACD.fakefakeACD('oryza.sdf')` to the python interpreter to create `fake_acd_oryza.sdf`
with dummy ^13^C NMR chemical shift values inside (99.99 for all carbon atoms).

- Create a new ACD DB named `fake_acd_oryza.NMRUDB`

- Import `fake_acd_oryza.sdf`

- In ACD software, run Database->Tools->Check Chemical Shifts

- Export DB as `fake_acd_oryza_exported.sdf` and close this DB

- Cut-and-paste `CNMR_predict.CNMR_predict('fake_acd_oryza_exported.sdf', 'true_acd_oryza.sdf')` to the python interpreter to create `true_acd_oryza.sdf`

- Create a new ACD database named `lotus_oryza.NMRUDB`

- Import `true_acd_oryza.sdf`

- In ACD software, run Database->Tools->Check Chemical Shifts

- Export DB as `lotus_oryza.sdf` and close DB (or keep it open to start doing something with it)

The DB named `lotus_oryza.NMRUDB` is ready for compound search according to ACD-predicted ^13^C NMR chemical shifts values.

The file `lotus_oryza.sdf.zip`, a zipped archive file of `lotus_oryza.sdf`, is available for reloading in an ACD database

Reminder file `templateB.txt` is almost identical to `templateA.txt` but uses ^13^C NMR chemical shift
prediction by nmrshiftdb2 instead of dummy values.
