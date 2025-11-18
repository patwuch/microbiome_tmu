# Core Laboratory of Human Microbiome @ Taipei Medical University

Our lab provides end-to-end services from storing biological samples to publication-ready reports. This includes step like DNA extraction, artefact creation, taxa classification, diversity + phylogenetic analysis, and functional metabolomics prediction. Currently we handle:
* Next Generation Sequencing (16S V3-V4 rRNA)
* Third Generation Sequencing (16S Full Length)
* Shotgun Metagenomics
* Anaerobic Bacterium Cultivation


For more details, please visit [our landing page](https://microbiome-in-tmu.mystrikingly.com).

## Table of Contents

This repository consists of three branches.
* __'main'__: This holds the analysis pipeline. Data provenance is tracked via Snakemake.
* __'work'__: This is the staging directory for both your raw data and your expected output plots, artefacts, tables. 
* __'references'__: This holds reference phylogenetic databases used by the analysis pipeline. Version control of this must be maintained by the user themselves.

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
