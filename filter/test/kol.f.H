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



Filter:   arustam@oas.stanford.edu   Mon Jul 30 13:28:45 2018
		sets next: in="/homes/sep/arustam/scratch/Filter@5YGE2P"
	REMARK: data and filter are different lengths
	REMARK: filter length will be padded for filtering
	REMARK: data length will be padded for filtering
		nt_after_padding_for_filtering=8192
		data_format="xdr_float"



bin/Kolmogorov.x:   arustam@oas.stanford.edu   Sun Aug  5 20:58:28 2018
		sets next: in="/homes/sep/arustam/scratch/Kolmogorov.x@2tVHsX"
	 	n1=2000 o1=0 d1=0.004 label1="sec" 

	 	n2=2 o2=0 d2=0.1 label2="km" 

	 	n3=1 o3=0 d3=1 label3="none" 

	 	n4=1 o4=0 d4=1 label4="none" 

	 	n5=1 o5=0 d5=1 label5="none" 

	 	n6=1 o6=0 d6=1 label6="none" 

	 	n7=1 o7=0 d7=1 label7="none" 

	 	n8=1 o8=0 d8=1 label8="none" 

		data_format="xdr_float"

