# Biomarkers and Shiny

This application requires R shiny to run.

Also, a "data" directory must be created in the same directory as this repository
with the following format:

- data/
    - tcga-dataset-name/
        - manifest-file
        - summary/
            - clinical.tsv
            - gdc_sample_sheet.2018-06-02.tsv
        - dir1/
        - dir2/
        - dir3/
        - dir4/
            - FPKM-UQ file from gdc


By running the gdc data transfer tool on the manifest
file under the dataset-name directory the dirX directories
will be created. Add the clinical data and sample sheet
in tsv format to a summary directory under the dataset directory as well.
