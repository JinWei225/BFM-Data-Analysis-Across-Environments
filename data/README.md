The dataset (bfm_data.csv) used in this analysis can be get from https://dx.doi.org/10.21227/bqsk-d072

`bfm_data.csv` contains one row per BFM report, with:

- `session_id` — unique identifier of a recording session
- `environment` — `open`, `foil`, or `nofoil` (shielding condition)
- `subject` — participant identifier (`collin`, `kenny`, `matthew`, `abel`, `ivan`)
- `activity` — `standing` or `walking`
- `timestamp` — capture time (UTC in the raw file)
- `receiver_address`, `destination_address`, `transmitter_address`, `source_address` — MAC addresses
- `SCIDX_<k>_Ratio_Mag` — per-subcarrier magnitude of the feedback matrix
- `SCIDX_<k>_..._Phase` — per-subcarrier phase of the feedback matrix