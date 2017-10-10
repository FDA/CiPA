# Milnes protocol data for 12 CiPA training drugs
Drug effects on dynamic hERG current block were measured using a modified Milnes protocol (Li *et al.* 2017, Milnes *et al.* 2010). These data were obtained at the U.S. Food and Drug Administration in White Oak, MD by Jiansong Sheng, Phu Tran, and Wendy Wu. The data have been postprocessed and formatted for model fitting by Sara Dutta and Kelly Chang.

## File format
Files are in comma-separated value (CSV) format with the following headers:

* **time**: time in msec during the sweep
* **frac**: fractional current
* **conc**: drug concentration in nM
* **exp**: experiment (cell) number
* **sweep**: the sweep number

Files should be saved in this directory as DRUG.csv, where "DRUG" is the drug name to be specified during fitting (case-sensitive).

## References
* Li, Z., Dutta, S., Sheng, J., Tran, P.N., Wu, W., Chang, K., et al. (2017). Improving the In Silico Assessment of Proarrhythmia Risk by Combining hERG (Human Ether-à-go-go-Related Gene) ChannelDrug Binding Kinetics and Multichannel Pharmacology. Circulation: Arrhythmia and Electrophysiology 10(2), e004628. doi: 10.1161/circep.116.004628.
* Milnes, J.T., Witchel, H.J., Leaney, J.L., Leishman, D.J., and Hancox, J.C. (2010). Investigating dynamic protocol-dependence of hERG potassium channel inhibition at 37 degrees C: Cisapride versus dofetilide. J Pharmacol Toxicol Methods 61(2), 178-191. doi: 10.1016/j.vascn.2010.02.007.
