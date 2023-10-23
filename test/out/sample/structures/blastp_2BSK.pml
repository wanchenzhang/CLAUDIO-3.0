load blastp_2BSK.pdb
hide all
bg_color white
set transparency, 0.8
zoom center, 50;
hide everything, show cartoon, chain C
show surface, chain C and blastp_2BSK
color 3, chain C
show cartoon, chain C
hide everything, show cartoon, chain D
show surface, chain D and blastp_2BSK
color 5, chain D
show cartoon, chain D
dist inter_18_5_valid , resid 58 and blastp_2BSK and chain C and name cb, resid 45 and blastp_2BSK and chain D and name cb
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
save blastp_2BSK.pse
png blastp_2BSK.png
