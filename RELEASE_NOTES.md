### Version 1.8.0
__Changes__
- updated dbCAN to v10

### Version 1.7.1
__Changes__
- tidying

### Version 1.7.0
__Changes__
- added phylogenetic marker library App

### Version 1.6.0
__Changes__
- updated HMMER to v3.3.2
- allow multiple target inputs to HMMER_Model_Group_Search()
- fixed bug with EnvBioelement 'backconverted subseq didn't end at expected length' by removing MSA output

### Version 1.5.0
__Changes__
- changed heatmap cells to DIV
- updated dbCAN to v8
- updated dbCAN defaults to match dbCAN2 server E-Value < 1e-15, coverage > 0.35
- added descriptions to dbCAN HMMs in key and link to key from column title
- added option to dbCAN to save hits to explicitly requested or ALL models in a category
- changed gene overlap_perc to model_cov_perc
- made addition of new model groups much easier. Method now general in Utils/HmmerUtils.py
- added EnvBioelement
- added genes hit to rollover behavior in report

### Version 1.4.5
__Changes__
- updated HMMER to v3.3.1

### Version 1.4.4
__Changes__
- added low_val option to heatmap containing methods to set color scale
- fixed bug with display of MSA Groups where columns shifted left

### Version 1.4.3
__Changes__
- add support for AnnotatedMetagenomeAssembly as target
- removed coalesce option for HMMER_Local_MSA_Group()

### Version 1.4.2
__Changes__
- add configuration for genome display name

### Version 1.4.1
__Changes__
- changed heatmap colors to make easier to distinguish and extend to black for max value

### Version 1.4.0
__Changes__
- updated HMMER to v3.3

### Version 1.3.0
__Changes__
- added App HMMER_EnvBioelement_Search()

### Version 1.2.2
__Changes__
- made protein sequence checking stricter
- changed base docker image to sdkbase2

### Version 1.2.1
__Changes__
- changed citation to PLOS format

### Version 1.2.0
__Changes__
- added App HMMER_dbCAN_Search()

### Version 1.1.2
__Changes__
- added HTML profile output for GenomeSet input to HMMER_Local_MSA_Group_Search()
- added overlap threshold to HMMER_MSA_Search() and HMMER_Local_MSA_Group_Search()

### Version 1.1.1
__Changes__
- added selectable MSAs to HMMER_Local_MSA_Group_Search()

### Version 1.1.0
__Changes__
- added HMMER_Local_MSA_Group_Search()

### Version 1.0.3
__Changes__
- changed kbase help contact

### Version 1.0.0
- Initial release
