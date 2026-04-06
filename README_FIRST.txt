FINAL CLEAN MATLAB PACKAGE FOR THE ECG ASSIGNMENT
=================================================

This folder is the clean final version.

WHAT TO KEEP
------------
Keep ONLY these files together in one MATLAB Online folder:
- all .m files from this package
- ecg_mother.mat
- Rva1.mat, Rva2.mat, Rva3.mat
- Ra1.mat, Ra2.mat, Ra3.mat
- indnonrecur.mat, indrecur.mat

DELETE OR IGNORE older mixed files from previous attempts.

ORDER TO RUN
------------
Best option: run one section at a time.
1) run_part11_pregnancy_ecg
2) run_part12_af_spectral
3) run_part21_synthetic_mranc
4) run_part22_fecg_mranc
5) run_part23_af_single_mranc
6) run_part24_af_full_mranc
7) run_part31_synthetic_bss
8) run_part32_fecg_bss
9) run_part33_af_single_bss
10) run_part34_af_full_bss

OR run everything with:
- run_all_assignment

FILES BY QUESTION
-----------------
Section 1.1  -> run_part11_pregnancy_ecg.m
Section 1.2  -> run_part12_af_spectral.m
Section 2.1  -> run_part21_synthetic_mranc.m
Section 2.2  -> run_part22_fecg_mranc.m
Section 2.3  -> run_part23_af_single_mranc.m
Section 2.4  -> run_part24_af_full_mranc.m
Section 3.1  -> run_part31_synthetic_bss.m
Section 3.2  -> run_part32_fecg_bss.m
Section 3.3  -> run_part33_af_single_bss.m
Section 3.4  -> run_part34_af_full_bss.m

IMPORTANT NOTES
---------------
1) AF SAMPLING FREQUENCY
The AF .mat files do not explicitly store fs.
The package currently uses:
    cfg.fsAF = 250;
inside assignment_config.m

Why 250 Hz?
Because each AF record has 15000 samples, so 250 Hz corresponds to 60 seconds.
If your professor gave another fs, change ONLY that single line in:
    assignment_config.m

2) SYNTHETIC SECTIONS
The exact "Exercises 5-7" synthetic definitions from the course were not pasted here.
So the synthetic scripts are fully runnable and answer the requested analyses,
but they are generic recreations, not guaranteed to match the professor's exact
course-generated synthetic signals sample-by-sample.

3) ICA
This package includes a self-contained ICA implementation:
    simple_fastica.m
So it works even if RobustICA is not installed.
If your professor strictly wants RobustICA, you can still use this package for
everything else and later replace the ICA line only.

OUTPUT FILES
------------
Each main script saves a result file:
- results_part11_pregnancy_ecg.mat
- results_part12_af_spectral.mat
- results_part21_synthetic_mranc.mat
- results_part22_fecg_mranc.mat
- results_part23_af_single_mranc.mat
- results_part24_af_full_mranc.mat
- results_part31_synthetic_bss.mat
- results_part32_fecg_bss.mat
- results_part33_af_single_bss.mat
- results_part34_af_full_bss.mat

MINIMUM SET IF YOU ARE SHORT ON TIME
------------------------------------
If you want the essential real-data sections first, run:
- run_part11_pregnancy_ecg
- run_part12_af_spectral
- run_part22_fecg_mranc
- run_part23_af_single_mranc
- run_part24_af_full_mranc
- run_part32_fecg_bss
- run_part33_af_single_bss
- run_part34_af_full_bss
