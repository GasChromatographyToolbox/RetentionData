## Naming convention

### Stationary phases

For the name of a stationary phase no fixed naming structure is given. Usually, the name of a stationary phase consists of an acronym of the manufacturer followed by an acronym for the stationary phase chemistry (often a number). This can be followed by additional features of the phase, like 'HT' for high temperature or 'MS' for mass spectrometer. If a new dataset is added, it should be checked, if the stationary phase is already in the database and in this case the same name should be used.

### Substances

For the assignment of the CAS number (and additional information in the GasChromatographySimulator.jl) the substance name is looked up in a database (ChemicalIdentifiers.jl). Not all synonyms of a substance are included in this database, especially non-english names. Therefore the names of substances should be the IUPAC name or the most common english name. 

Using an IDList (a list of at least two columns, first column has the shortname or number (ID), the second column has a full name or alternative name, which can be found with ChemicalIdentifiers.jl) can help to rename substances without changing the names in the original data (e.g. in the `lnk-T` data). The file should be named `Source_IDList.csv`.

Not all substances are included in ChemicalIdentifiers.jl. For substances for which this is the case the file `missing.csv` in the `data` folder is used, where the substance name together with the needed information (CAS, formula, molar weight, SMILES) is noted.