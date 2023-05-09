
export void lanczos_fwd(uniform int ib1,uniform int ie1,
												uniform float filter3,uniform float filter2, uniform float filter1[],
                        uniform float mod_re[], uniform float dat_re[]) {

    float re = 0;
    float im = 0;
    int index;

    foreach (j1=ib1...ie1) {
      if (filter1[j1] != 0) {
        re += mod_re[2*j1] * filter3 * filter2 * filter1[j1];
        im += mod_re[2*j1+1] * filter3 * filter2 * filter1[j1];
      }
    }
    dat_re[0] += reduce_add(re);
    dat_re[1] += reduce_add(im);
  }

export void lanczos_adj(uniform int ib1,uniform int ie1,
												uniform float filter3,uniform float filter2, uniform float filter1[],
                        uniform float mod_re[], uniform float dat_re[]) {

    float re = 0;
    float im = 0;

    foreach (j1=ib1...ie1) {
      if (filter1[j1] != 0) {
        mod_re[2*j1] += dat_re[0] * filter3 * filter2 * filter1[j1];
        mod_re[2*j1+1] += dat_re[1] * filter3 * filter2 * filter1[j1];
      }
    }
  }



	export void lanczos_fwd_2d(uniform int ib1,uniform int ie1,
													uniform float filter2, uniform float filter1[],
	                        uniform float mod_re[], uniform float dat_re[]) {

	    float re = 0;
	    float im = 0;
	    int index;

	    foreach (j1=ib1...ie1) {
	      if (filter1[j1] != 0) {
	        re += mod_re[2*j1] * filter2 * filter1[j1];
	        im += mod_re[2*j1+1] * filter2 * filter1[j1];
	      }
	    }
	    dat_re[0] += reduce_add(re);
	    dat_re[1] += reduce_add(im);
	  }

	export void lanczos_adj_2d(uniform int ib1,uniform int ie1,
													uniform float filter2, uniform float filter1[],
	                        uniform float mod_re[], uniform float dat_re[]) {

	    float re = 0;
	    float im = 0;

	    foreach (j1=ib1...ie1) {
	      if (filter1[j1] != 0) {
	        mod_re[2*j1] += dat_re[0] * filter2 * filter1[j1];
	        mod_re[2*j1+1] += dat_re[1] * filter2 * filter1[j1];
	      }
	    }
	  }
