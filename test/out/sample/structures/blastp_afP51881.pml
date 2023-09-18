load blastp_afP51881.pdb
hide all
bg_color white
set transparency, 0.8
zoom center, 50;
hide everything, show cartoon, chain A
show surface, chain A and blastp_afP51881
color 3, chain A
show cartoon, chain A
dist intra_25_valid , resid 147 and blastp_afP51881 and chain A and name cb, resid 33 and blastp_afP51881 and chain A and name cb
dist intra_26_valid , resid 166 and blastp_afP51881 and chain A and name cb, resid 155 and blastp_afP51881 and chain A and name cb
dist intra_27_valid , resid 92 and blastp_afP51881 and chain A and name cb, resid 23 and blastp_afP51881 and chain A and name cb
dist intra_28_valid , resid 147 and blastp_afP51881 and chain A and name cb, resid 63 and blastp_afP51881 and chain A and name cb
dist intra_29_valid , resid 147 and blastp_afP51881 and chain A and name cb, resid 245 and blastp_afP51881 and chain A and name cb
dist intra_30_valid , resid 92 and blastp_afP51881 and chain A and name cb, resid 105 and blastp_afP51881 and chain A and name cb
dist intra_31_valid , resid 49 and blastp_afP51881 and chain A and name cb, resid 147 and blastp_afP51881 and chain A and name cb
dist intra_32_valid , resid 33 and blastp_afP51881 and chain A and name cb, resid 49 and blastp_afP51881 and chain A and name cb
dist intra_33_valid , resid 147 and blastp_afP51881 and chain A and name cb, resid 43 and blastp_afP51881 and chain A and name cb
dist intra_34_valid , resid 163 and blastp_afP51881 and chain A and name cb, resid 245 and blastp_afP51881 and chain A and name cb
dist intra_35_valid , resid 166 and blastp_afP51881 and chain A and name cb, resid 245 and blastp_afP51881 and chain A and name cb
dist intra_36_valid , resid 43 and blastp_afP51881 and chain A and name cb, resid 245 and blastp_afP51881 and chain A and name cb
dist intra_37_valid , resid 147 and blastp_afP51881 and chain A and name cb, resid 163 and blastp_afP51881 and chain A and name cb
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
save blastp_afP51881.pse
png blastp_afP51881.png
