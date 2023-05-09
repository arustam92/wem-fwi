
export void ps_fwd(uniform int nk, uniform float k[], uniform float w, uniform float dz,
uniform float s_re, uniform float s_im, uniform float eps,
uniform float mod_re[], uniform float dat_re[]) {

    	float a, b, re, im, c;
    foreach (ik=0...nk) {
			a = w*w*s_re- k[ik]*k[ik];
			b = w*w*(s_im-eps*s_re);
			c = sqrt(a*a + b*b);
			re = sqrt((c+a)/2);
			im = -sqrt((c-a)/2);
			dat_re[2*ik] += exp(im*dz) * (mod_re[2*ik]*cos(re*dz) + mod_re[2*ik+1]*sin(re*dz));
			dat_re[2*ik+1] += exp(im*dz) * (-mod_re[2*ik]*sin(re*dz) + mod_re[2*ik+1]*cos(re*dz));
		    }
  }

export void ps_adj(uniform int nk, uniform float k[], uniform float w, uniform float dz,
uniform float s_re, uniform float s_im, uniform float eps,
uniform float mod_re[], uniform float dat_re[]) {

    	float a, b, re, im, c;
    foreach (ik=0...nk) {
			a = w*w*s_re - k[ik]*k[ik];
			b = w*w*(s_im-eps*s_re);
			c = sqrt(a*a + b*b);
			re = sqrt((c+a)/2);
			im = -sqrt((c-a)/2);
			mod_re[2*ik] += exp(im*dz) * (dat_re[2*ik]*cos(re*dz) - dat_re[2*ik+1]*sin(re*dz));
			mod_re[2*ik+1] += exp(im*dz) * (dat_re[2*ik]*sin(re*dz) + dat_re[2*ik+1]*cos(re*dz));
		    }
  }
