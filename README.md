# pp224_phasing
### Programming Project for Haplotype phasing

## Command line usage:
```
 python3 segment_phasing.py <input_genotypes> <step> <overlap>
```
### Example 1, 2, or 3:
```
 python3 src/segment_phasing.py data/example_data_1.txt 5 5
```
### Test Data Set 1:
```
 python3 src/segment_phasing.py data/test_data_1.txt 5 5
```
### Test Data Set 2:
```
 python3 src/segment_phasing.py data/test_data_2.txt 5 5
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
