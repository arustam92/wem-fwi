Input from a tty

Spike:   arustam@oas.stanford.edu   Mon Jul 30 13:27:49 2018
		sets next: in="/homes/sep/arustam/scratch/Spike@nRp6G0"
		n1=2000
		n2=2
		n3=1
		n4=1
		n5=1
		n6=1
		n7=1
		o1=0
		o2=0
		o3=0
		o4=0
		o5=0
		o6=0
		o7=0
		d1=0.004
		d2=0.1
		d3=0.1
		d4=0.1
		d5=0.1
		d6=0.1
		d7=0.1
		label1="sec"
		label2="km"
		label3="km"
		label4="km"
		label5="km"
		title="Impulsive_plane_wave"
		esize=4
		data_format="xdr_float"



Filter:   arustam@oas.stanford.edu   Mon Jul 30 13:28:51 2018
		sets next: in="/homes/sep/arustam/scratch/Filter@hZJy3G"
	REMARK: data and filter are different lengths
	REMARK: filter length will be padded for filtering
	REMARK: data length will be padded for filtering
		nt_after_padding_for_filtering=8192
		data_format="xdr_float"



Pad:   arustam@oas.stanford.edu   Sun Aug  5 20:58:28 2018
		sets next: in="/homes/sep/arustam/scratch/Pad@GWFbFV"
	 	n1=2000 o1=0 d1=0.004 label1="" 

	 	n2=2 o2=0 d2=0.1 label2="" 

		data_format="xdr_float"



Filter:   arustam@oas.stanford.edu   Sun Aug  5 20:58:28 2018
		sets next: in="/homes/sep/arustam/scratch/Filter@PYz3iY"
	WARNING: data and filter have different d1's, using data's
	REMARK: filter length will be padded for filtering
	REMARK: data length will be padded for filtering
		nt_after_padding_for_filtering=8192
		data_format="xdr_float"

