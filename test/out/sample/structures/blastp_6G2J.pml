load blastp_6G2J.pdb
hide all
bg_color white
set transparency, 0.8
zoom center, 50;
hide everything, show cartoon, chain O
show surface, chain O and blastp_6G2J
color 3, chain O
show cartoon, chain O
hide everything, show cartoon, chain P
show surface, chain P and blastp_6G2J
color 5, chain P
show cartoon, chain P
hide everything, show cartoon, chain W
show surface, chain W and blastp_6G2J
color 6, chain W
show cartoon, chain W
dist intra_24_valid , resid 144 and blastp_6G2J and chain P and name cb, resid 154 and blastp_6G2J and chain P and name cb
dist inter_15_valid , resid 208 and blastp_6G2J and chain O and name cb, resid 89 and blastp_6G2J and chain W and name cb
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
save blastp_6G2J.pse
png blastp_6G2J.png
