load blastp_7UOL.pdb
hide all
bg_color white
set transparency, 0.8
zoom center, 50;
hide everything, show cartoon, chain E
show surface, chain E and blastp_7UOL
color 3, chain E
show cartoon, chain E
dist intra_23_26_valid , resid 274 and blastp_7UOL and chain E and name cb, resid 279 and blastp_7UOL and chain E and name cb
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
save blastp_7UOL.pse
png blastp_7UOL.png
