# Batch processing

Using the `pdemtools` suite of functions, it is possible to efficiently download and preprocess batches of ArcticDEM and REMA imagery.

An example of a batch download/coregistration script for ArcticDEM strips is [made available at the `pdemtools` GitHub](https://github.com/trchudley/pdemtools/blob/main/batch/batch_download_and_coregister.py). Coregistration is performed against the ArcticDEM mosiac, which is in turn coregistered against ICESat-2.