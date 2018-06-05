# pp224_phasing
### Programming Project for Haplotype phasing

## Command line usage:
```
 python segment_phasing.py <input_genotypes> <step> <overlap>
```
### Example:
```
 python src/segment_phasing.py data/ex1_3k.txt 5 5
```
Note: Most likely haplotypes are calculated for each segment. 
	segment_length = step + overlap

## To test switching accuracy:
```
 Rscript calculate_switch_accuracy.R <expected_haplotypes> <real_haplotypes>
```
### Example:
```
 /usr/local/bin/Rscript src/calculate_switch_accuracy.R data/ex2_3k_seg10_ovl5.txt data/ex2_3k_s.txt 
```