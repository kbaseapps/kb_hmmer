### Version 1.5.0
__Changes__
- DEBUGGING changed heatmap cells to DIV
- DEBUGGING updated dbCAN to v8
- DEBUGGING updated dbCAN defaults to match dbCAN2 server E-Value < 1e-15, coverage > 0.35
- DEBUGGING added descriptions to dbCAN HMMs in key and link to key from column title
- DEBUGGING added option to dbCAN to save hits to explicitly requested or ALL models in a category
- PENDING added overlap config and rules to allow multiple non-overlapping hits to same gene
- PENDING create DomainAnnotation object
- PENDING added HMM-specific thresholds in config
- PENDING added genes hit to rollover behavior in report

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
