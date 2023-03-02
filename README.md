# RetentionData

[![DOI](https://zenodo.org/badge/426138381.svg)](https://zenodo.org/badge/latestdoi/426138381)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JanLeppert.github.io/RetentionData/dev)

This github is a compilation of retention data in gas chromatography (GC).

**Short introduction to retention factor and the 3 models**

## Databases

Different sets of parameters (A, B, C parameters, thermodynamic parameters or K-centric parameters) from different sources are compiled 
in the folder `Databases`. For every source a separate folder is created containing the original data from the sources as comma separated files, reference information as a BibTex file and processed files (e.g. calculation of the other sets of parameters). In the sub-folder `Measurements` in the folder `Databases` measured retention data, which is not publicized yet, is collected. 

[Latest database](https://raw.githubusercontent.com/JanLeppert/RetentionData/main/Databases/database_all.csv)

[Latest database without flagged substances](https://raw.githubusercontent.com/JanLeppert/RetentionData/main/Databases/database_nonflag.csv)

### How to add data
_preliminary_

* contact us
 
* fork the Github

* New data in form of parameter sets (A,B,C) and/or (Tchar, thetachar, DeltaCp) and/or (DeltaS, DeltaH, DeltaCp, Tref) -> data in certain structure? and folders -> notebook 'Convert_Parameters.jl' and save 'AllParam' files

* New data in form of ln(k) over T: ...

* alternative: add the data in folders with structure ... -> run script.jl -> new database constructed (two files, one with and one without flagged substances)

### Structure of filenames and files:

see [Documentation](https://janleppert.github.io/RetentionData/dev/filestructure/)

## Overview of Notebooks of the Project _RetentionData_

Several notebooks, written in the Julia programming language (https://julialang.org/) and using the Pluto.jl package (https://github.com/fonsp/Pluto.jl), are used to
calculate the retention parameters from retention data (`Fit_lnk-T_Outliertest.jl`) or convert available sets of retention parameters into the other sets (`Convert_Parameters.jl`). These notebooks can be found in the folder `notebooks`. Functions, which are used in these notebooks can be found in the file `RetentionData.jl` in the `src` folder. 

- `Convert_Parameters.jl`:
  - Load all files with the keyword `Parameters` from a folder
  - convert the parameter set from the file into the other parameter sets
  - ToDo: if no column for `CAS` number exists in the data, add one with the `CAS`-numbers from ChemicalIdentifiers.jl (also by using the file `shortname.csv`) resp. from the file `missing.csv` (chemical substances, which are not found in the database of ChemicalIdentifiers.jl)
  - optional: save the data of all three parameter sets into a file with the keyword `AllParam`

- `Fit_lnk-T_Outliertest.jl`
  - Load all files with the keyword `lnk-T` resp. `log10k-T` from a folder
  - fit of the ABC-model and the Kcentric model to the data and convert them to the third parameter set
  - export the parameter sets together with some statistics of the fits
  - ToDo: if no column for `CAS` number exists in the data, add one with the `CAS`-numbers from ChemicalIdentifiers.jl (also by using the file `shortname.csv`) resp. from the file `missing.csv` (chemical substances, which are not found in the database of ChemicalIdentifiers.jl)
  - optional: save the data of all three parameter sets into a file with the keyword `AllParam`

- `Notebook_PCA.jl`: 
  - ...

- `Database.jl`: 
  - ...
  




