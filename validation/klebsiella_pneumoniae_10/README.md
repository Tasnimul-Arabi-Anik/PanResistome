# Klebsiella pneumoniae 10-Genome v0.5.0 Docker/GHCR Validation

This validation subset contains the first 10 records from the BioProject-diverse
100-record `Klebsiella pneumoniae` validation input. It is intended as a compact
biological validation for v0.5.0 interpretation outputs, not as a replacement
for the 100-record or 300-record scale validations.

The run uses the public GHCR Docker image and keeps heavyweight stages bounded:
CheckM2, GTDB-Tk, QUAST, and ANI are disabled; Mash, AMRFinderPlus, native PanR2
feature runners, PanR2 comprehensive reporting, and PanR2 handoff export are
enabled.

See `V0_5_0_DOCKER_VALIDATION_RESULTS.md` for the completed validation result.
