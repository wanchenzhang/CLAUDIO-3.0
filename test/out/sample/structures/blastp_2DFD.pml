load blastp_2DFD.pdb
hide all
bg_color white
set transparency, 0.8
zoom center, 50;
hide everything, show cartoon, chain A
show surface, chain A and blastp_2DFD
color 1, chain A
show cartoon, chain A
hide everything, show cartoon, chain B
show surface, chain B and blastp_2DFD
color 3, chain B
show cartoon, chain B
hide everything, show cartoon, chain C
show surface, chain C and blastp_2DFD
color 5, chain C
show cartoon, chain C
hide everything, show cartoon, chain D
show surface, chain D and blastp_2DFD
color 6, chain D
show cartoon, chain D
dist intra_1_same , resid 87 and blastp_2DFD and chain A and name cb, resid 87 and blastp_2DFD and chain A and name cb
dist inter_1_valid , resid 87 and blastp_2DFD and chain A and name cb, resid 87 and blastp_2DFD and chain B and name cb
dist inter_2_valid , resid 87 and blastp_2DFD and chain A and name cb, resid 87 and blastp_2DFD and chain C and name cb
dist inter_3_valid , resid 87 and blastp_2DFD and chain A and name cb, resid 87 and blastp_2DFD and chain D and name cb
dist inter_4_valid , resid 87 and blastp_2DFD and chain B and name cb, resid 87 and blastp_2DFD and chain C and name cb
dist inter_5_valid , resid 87 and blastp_2DFD and chain B and name cb, resid 87 and blastp_2DFD and chain D and name cb
dist inter_6_valid , resid 87 and blastp_2DFD and chain C and name cb, resid 87 and blastp_2DFD and chain D and name cb
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
save blastp_2DFD.pse
png blastp_2DFD.png
