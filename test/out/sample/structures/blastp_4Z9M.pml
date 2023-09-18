load blastp_4Z9M.pdb
hide all
bg_color white
set transparency, 0.8
zoom center, 50;
hide everything, show cartoon, chain A
show surface, chain A and blastp_4Z9M
color 3, chain A
show cartoon, chain A
hide everything, show cartoon, chain B
show surface, chain B and blastp_4Z9M
color 5, chain B
show cartoon, chain B
hide everything, show cartoon, chain C
show surface, chain C and blastp_4Z9M
color 6, chain C
show cartoon, chain C
hide everything, show cartoon, chain D
show surface, chain D and blastp_4Z9M
color 7, chain D
show cartoon, chain D
hide everything, show cartoon, chain E
show surface, chain E and blastp_4Z9M
color 8, chain E
show cartoon, chain E
hide everything, show cartoon, chain F
show surface, chain F and blastp_4Z9M
color 9, chain F
show cartoon, chain F
hide everything, show cartoon, chain G
show surface, chain G and blastp_4Z9M
color 10, chain G
show cartoon, chain G
hide everything, show cartoon, chain H
show surface, chain H and blastp_4Z9M
color 11, chain H
show cartoon, chain H
dist intra_21_valid , resid 292 and blastp_4Z9M and chain A and name cb, resid 140 and blastp_4Z9M and chain A and name cb
dist intra_21_10_valid , resid 292 and blastp_4Z9M and chain B and name cb, resid 140 and blastp_4Z9M and chain B and name cb
dist intra_21_19_valid , resid 292 and blastp_4Z9M and chain C and name cb, resid 140 and blastp_4Z9M and chain C and name cb
dist intra_21_28_valid , resid 292 and blastp_4Z9M and chain D and name cb, resid 140 and blastp_4Z9M and chain D and name cb
dist intra_21_37_valid , resid 292 and blastp_4Z9M and chain E and name cb, resid 140 and blastp_4Z9M and chain E and name cb
dist intra_21_46_valid , resid 292 and blastp_4Z9M and chain F and name cb, resid 140 and blastp_4Z9M and chain F and name cb
dist intra_21_55_valid , resid 292 and blastp_4Z9M and chain G and name cb, resid 140 and blastp_4Z9M and chain G and name cb
dist intra_21_64_valid , resid 292 and blastp_4Z9M and chain H and name cb, resid 140 and blastp_4Z9M and chain H and name cb
dist intra_22_same , resid 353 and blastp_4Z9M and chain A and name cb, resid 353 and blastp_4Z9M and chain A and name cb
dist intra_22_19_same , resid 353 and blastp_4Z9M and chain C and name cb, resid 353 and blastp_4Z9M and chain C and name cb
dist intra_22_28_same , resid 353 and blastp_4Z9M and chain D and name cb, resid 353 and blastp_4Z9M and chain D and name cb
dist intra_22_37_same , resid 353 and blastp_4Z9M and chain E and name cb, resid 353 and blastp_4Z9M and chain E and name cb
dist intra_22_46_same , resid 353 and blastp_4Z9M and chain F and name cb, resid 353 and blastp_4Z9M and chain F and name cb
dist intra_22_64_same , resid 353 and blastp_4Z9M and chain H and name cb, resid 353 and blastp_4Z9M and chain H and name cb
show dashes
set dash_gap, 0.1
color blue, intra*_valid
color red, intra*_out_range
color red, intra*_overlaps
color red, intra*_same
color blue, inter*_valid
color red, inter*_out_range
color red, inter*_overlaps
color red, inter*_same
set dash_width, 9
center
save blastp_4Z9M.pse
png blastp_4Z9M.png
