load blastp_2W6E.pdb
hide all
bg_color white
set transparency, 0.8
zoom center, 50;
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
save blastp_2W6E.pse
png blastp_2W6E.png
