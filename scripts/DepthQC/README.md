# Optional base-level depth-ratio utilities

`RunSampleBaseDepthRatio.pl` and `CalcDepthRatio.pl` are retained as optional
diagnostic utilities. They are not called by `HCMExonDel.pl` and do not create
standard S2 outputs.

The utilities can be run manually against the per-gene files in
`S2_Depth/gene_depth/`. For example:

```bash
perl scripts/DepthQC/RunSampleBaseDepthRatio.pl \
  --InDir results/SAMPLE/S2_Depth/gene_depth \
  --OutDir results/SAMPLE/DepthQC/base_depth_ratio \
  --Perl /usr/bin/perl \
  --Script scripts/DepthQC/CalcDepthRatio.pl
```

These optional outputs are intended for detailed inspection or custom
visualization. They are not used by the standard candidate-evidence merge.
