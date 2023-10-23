load blastp_4Z9M.pdb
hide all
bg_color white
set transparency, 0.8
zoom center, 50;
hide everything, show cartoon, chain B
show surface, chain B and blastp_4Z9M
color 3, chain B
show cartoon, chain B
hide everything, show cartoon, chain A
show surface, chain A and blastp_4Z9M
color 5, chain A
show cartoon, chain A
dist intra_21_10_valid , resid 292 and blastp_4Z9M and chain B and name cb, resid 140 and blastp_4Z9M and chain B and name cb
dist intra_22_same , resid 353 and blastp_4Z9M and chain A and name cb, resid 353 and blastp_4Z9M and chain A and name cb
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
