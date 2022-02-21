# Note

## Structure of filenames:

The general structure of filenames in this Database is the following:

Source\_Datatype\_(Substances)or(Table)\_StationaryPhase\_PhaseRatio\_(AddParameters).csv

with:
- Source: Name of the source, identical with the key of the Bibtex file `references.bib` (name of first author and year). For not yet publicized datasets the name consists of the name of the primary person which measured the data and the year.
- Datatype: Which data is in the file
  - `Parameters`: One of the following sets of thermodynamic parameters:
    - `A, B, C`
    - `DeltaHref`, `DeltaSref`, `DeltaCp`
    - `Tchar`, `thetachar`, `DeltaCp`
  - `lnk_T`: The measured `lnk`-values (logarithm of the retention factor) of the substances for isothermal GC measurements at defined temperatures `T`. The thermodynamic parameters can be estimated from these measurements by fitting of the model.
  - `tR_T`: The measured `tR`-values (retention time) of the substances for isothermal GC measurements at defined temperatures `T`. Also, the hold-up time `tM` must be included (measured or calculated), to calculate `lnk` from this data and estimate the thermodynamic parameters by model fitting. 
- (Substances): Optional. A label for the substances in the file, e.g. `Mix1` or `Alkanes`. 
- (Table): Table + Number of the table in the source.
- StationaryPhase: The name of the stationary phase. Spaces or hyphens are removed from the name. The name is follows in general the name given by the manufacturer.
- PhaseRatio: The nominal phase ratio of the used column according the the data of the manufacturer. Calculated by the approximation $\beta \approx 1/4 d/d_f$
- (AddParameters): Additional parameters (the value follows the parameters without spaces):
  - `Tref`: Reference temperature for the thermodynamic parameters `DeltaHref` and `DeltaSref`
  - `gas`: The type of gas used as stationary phase for the measurement of isothermal chromatograms or thermodynamic parameters, e.g. `He` or `H2`.
  - `d`: The diameter of the column used for the measurements of isothermal chromatograms or thermodynamic parameters, e.g. `0.10` or `0.32`.

