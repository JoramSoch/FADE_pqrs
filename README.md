# FADE_pqrs

**Computational modelling of subsequent memory reports**

This code belongs to the preprint "A novel approach for modelling subsequent memory reports by separating decidedness, recognition and confidence" by Soch et al. (2022), submitted to *Journal of Experimental Psychology: Learning, Memory, and Cognition* and publicly available from *bioRxiv*. For instructions how to process these data, see below.

- Preprint: https://psyarxiv.com/u5t82/
- Data: https://github.com/JoramSoch/FADE_pqrs/blob/main/data/logfiles_FADE.mat
- Code: https://github.com/JoramSoch/FADE_pqrs
- Toolbox: https://github.com/JoramSoch/pqrs


### Requirements

This code was developed and run using the following software:
- Windows 10 Professional 64-bit
- [MATLAB R2018a/R2021a](https://de.mathworks.com/help/matlab/release-notes.html) (Version 9.4/9.10)
- [MATLAB Statistics and ML Toolbox](https://de.mathworks.com/products/statistics.html) (Version 11.7)
- [pqrs](https://github.com/JoramSoch/pqrs) package (as on GitHub)
- [bonf_holm](https://www.mathworks.com/matlabcentral/fileexchange/28303-bonferroni-holm-correction-for-multiple-comparisons) function (Version 1.1.0.0)


### Instructions

For re-running analyses reported in the paper, you need to perform the following steps:
1. Create a folder on your computer that hereafter is referred to as the "study directory".
2. Make sure that you have downloaded the pqrs package to some folder on your computer (see above).
3. Download the [analysis scripts](https://github.com/JoramSoch/FADE_pqrs/archive/main.zip) and place them into a sub-folder of the study directory called "FADE_pqrs".
4. Open MATLAB, set your current directory to this sub-folder and (if necessary) edit the pqrs directory [in line 13](https://github.com/JoramSoch/FADE_pqrs/blob/main/Figures_res_all.m#L13) of `Figures_res_all.m`.
5. Finally, open `analyses_FADE_pqrs.m` and run this script.

* When running this, the Figures 1B, 5, 6, 7, 8 and 9, as they appear in the paper, should be displayed.
* Moreover, the summary statistics shown on Figure 10 should also be printed to the command window.
* For creating Figures 7 and 8, you need to download the MATLAB function `bonf_holm.m` from [MATLAB Central](https://www.mathworks.com/matlabcentral/fileexchange/28303-bonferroni-holm-correction-for-multiple-comparisons).