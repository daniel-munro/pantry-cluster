"""Add age or sex metadata to table of Pantry covariates"""

import polars as pl

meta_fname = "GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
modalities = [
    "all",
    "alt_polyA",
    "alt_TSS",
    "expression",
    "isoforms",
    "splicing",
    "stability",
]
tissue = "ADPSBQ"

# Load metadata (infer_schema_length=0 loads everything as string)
meta = pl.read_csv(
    meta_fname, separator="\t", columns=["SUBJID", "SEX", "AGE"], infer_schema_length=0
).rename({"SUBJID": "individual", "SEX": "sex", "AGE": "age"})

# Load covariates
for modality in modalities:
    covar_fname = f"{tissue}/intermediate/covar/{modality}.nometa.covar.tsv"
    covar = pl.read_csv(covar_fname, separator="\t", infer_schema_length=0)
    # Columns are covar ID and then individual IDs. Add metadata covar as a new row
    for name in ['age']:#['sex', 'age']:
        to_add = (
            meta.filter(pl.col("individual").is_in(covar.columns[1:]))
            .with_columns(pl.lit(name).alias("ID"))
            .pivot(on="individual", values=name, index="ID")
        )
        covar2 = covar.vstack(to_add)
        covar2.write_csv(
            f"{tissue}/intermediate/covar/{modality}.{name}.covar.tsv", separator="\t"
        )
