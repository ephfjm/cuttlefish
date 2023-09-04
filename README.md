# Cuttlefish.r: An Rscript for Structural Variant Tool Comparison

## Quickstart

```
Rscript cuttlefish.r \
tool_1.vcf \
tool_2.vcf \
tool_1_name \
tool_2_name \
start_coord_bool \
end_coord_bool \
coord_buffer \
output_dir
```

Currently only support VCFs produced by Sniffles2, CuteSV and NanoVar, with accepted "tool names" as "snif", "cutesv" and "nv", respectively.

- tool_1.vcf - VCF file produced by tool 1
- tool_2.vcf - VCF file produced by tool 2
- tool_1_name - snif/cutesv/nv
- tool_2_name - snif/cutesv/nv
- start_coord_bool - TRUE/FALSE (Enable or diasble SV start coordinate comparison)
- end_coord_bool - TRUE/FALSE (Enable or diasble SV end coordinate comparison)
- coord_buffer - INT (Number of bases to create a buffer window for each SV start/end coordinate for comparison. e.g. if coord_buffer=50, then start1-50 <= start2 <= start1+50 needs to be satisfied for a hit)
- output_dir - PATH (Directory to save output files)

Example run:
```
Rscript cuttlefish.r
sniffles.vcf \
nanovar.vcf \
snif \
nv \
TRUE \
FALSE \
50 \
./output
```

## Output
- sim_dat1.csv (Similar SV hits from VCF file 1)
- unique_dat1.csv (Unique SVs from VCF file 1)
- sim_dat2.csv (Similar SV hits from VCF file 2)
- unique_dat2.csv (Unique SVs from VCF file 2)
