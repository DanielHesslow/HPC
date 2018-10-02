#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"
#include "math.h"



typedef struct Cpx
{
	double re;
	double im;
} Cpx;


Cpx cpx_mul(Cpx a, Cpx b)
{
	return (Cpx) {
		.re = a.re * b.re - a.im * b.im,
		.im = a.re * b.im + b.re * a.im
	};
}

Cpx cpx_add(Cpx a, Cpx b)
{
	return (Cpx){
		.re = a.re + b.re,
		.im = a.im + b.im,
	};
}

Cpx cpx_sub(Cpx a, Cpx b)
{
	return (Cpx){
		.re = a.re - b.re,
		.im = a.im - b.im,
	};
}


Cpx cpx_pow(Cpx x, int p)
{
	Cpx ret = {.re = 1,.im = 0};
	for(int i = 0; i<p;i++) 
	{
		ret = cpx_mul(ret, x);
	}
	return ret;
}

Cpx cpx_pow_fast(Cpx x, int p)
{
	Cpx ret = {.re = 1,.im = 0};
	Cpx x_pow_2 = x;
	// calculate the powers of two's separately then mul together. 
	while(p)
	{
		if(p&1) ret = cpx_mul(ret,x_pow_2);
		x_pow_2 = cpx_mul(x_pow_2, x_pow_2);	
		p>>=1;
	}
	return ret;
}

Cpx cpx_pow_fast_2(Cpx x, int p)
{
	//avoids the mul with 1+0i
	if(!p) return (Cpx){.re = 1,.im = 0};
	Cpx x_pow = x;
	while(!(p&1))
	{
		x_pow = cpx_mul(x_pow, x_pow);	
		p>>=1;
	}
	
	Cpx ret = x_pow;
	x_pow = cpx_mul(x_pow, x_pow);	
	p>>=1;
	while(p)
	{
		if(p&1) ret = cpx_mul(ret,x_pow);
		x_pow = cpx_mul(x_pow, x_pow);	
		p>>=1;
	}
	return ret;
}

Cpx cpx_pow_unrolled(Cpx x, int p)
{
	Cpx xx;
	Cpx xxxx;
	switch(p){
		case 0: return (Cpx){.re=1,.im=0};
		case 1: return x;
		case 2: return cpx_mul(x,x);
		case 3: return cpx_mul(cpx_mul(x,x),x);
		case 4: xx = cpx_mul(x,x); return cpx_mul(xx,xx);
		case 5: xx = cpx_mul(x,x); return cpx_mul(cpx_mul(xx,xx),x);
		case 6: xx = cpx_mul(x,x); return cpx_mul(cpx_mul(xx,xx),xx);
		case 7: xx = cpx_mul(x,x); return cpx_mul(cpx_mul(xx,xx),cpx_mul(xx,x));
		case 8: xx = cpx_mul(x,x); xxxx = cpx_mul(xx,xx); return cpx_mul(xxxx,xxxx);
		default: return cpx_pow_fast_2(x,p);
	}
}

double cpx_magnitude_sq(Cpx x)
{
	return x.re *x.re + x.im*x.im;
}
	
Cpx cpx_mul_real(Cpx x, double r)
{
	x.re *= r;
	x.im *= r;
	return x;
}

Cpx cpx_div_real(Cpx x, double r)
{
	return cpx_mul_real(x, 1.0/r);
}

Cpx newton_iteration(int degree, Cpx x)
{
	Cpx x_pow = cpx_pow_unrolled(x,degree-1);
	Cpx fx, min_inv_dfx;
	{
		fx = cpx_mul(x_pow, x);
		fx.re -= 1;
	}

	{
		Cpx d_x_pow = cpx_mul_real(x_pow, degree);
		min_inv_dfx = cpx_div_real(d_x_pow, cpx_magnitude_sq(d_x_pow));
		min_inv_dfx.re = -min_inv_dfx.re;
	}

	return cpx_add(x, cpx_mul(fx, min_inv_dfx));
}

void print_cpx(Cpx c){
	printf("%f%+fi\n",c.re,c.im);
}


double cpx_distance_sq(Cpx a, Cpx b){
	return cpx_magnitude_sq(cpx_sub(a,b));
}



int main(void)
{
	int degree = 7;
	int l = 4000;

	const double pi = M_PI;

	int *num_its= malloc(sizeof(int)*l*l); 
	int *root   = malloc(sizeof(int)*l*l);
	Cpx *xs      = malloc(sizeof(Cpx)*l*l);
	Cpx *correct_roots = malloc(sizeof(Cpx)*degree);
	for(int i = 0; i<degree;i++) {
		correct_roots[i] = (Cpx){cos(2.0*pi*i/(double)degree),sin(2.0*pi*i/(double)degree)};
	}


	for(int x = 0; x < l;x++){
		for(int y = 0; y < l;y++){
			xs[x+y*l] = (Cpx){.re = 4*(x / (double)l-0.5), 4*(y/(double)l - 0.5)};
		}
	}

	
	for(int px = 0; px< l*l;px++){
		for(int its = 0;; its++){
			int done = 0;
			for(int j=0;j<degree;j++) {
				if(cpx_distance_sq(xs[px], correct_roots[j]) < 1e-6){
					root[px] = j+1;
					num_its[px] = its+1;
					done = 1;
					break;
				}
			}
			if(done) break;
			if(fabs(xs[px].re) > 1e10 || fabs(xs[px].im) > 1e10 || cpx_magnitude_sq(xs[px]) < 1e-6 )
			{
				root[px] = degree+1;
				num_its[px] = its+1;
				break;
			}

			xs[px] = newton_iteration(degree,xs[px]);
		}
	}



	char *cols[] = {"0 0 0", "255 0 0", "0 255 0", "0 0 255","255 0 255","0 255 255","255 255 0",  "255 255 255", "200 100 150", "200 150 100", "150 200 100"};

	FILE *root_file= fopen("newton_attractors_xd.ppm", "w");
	fprintf(root_file, "P3\n%d %d\n%d\n", l, l, 255);
	for(int i = 0; i< l; i++){
		for(int j = 0; j< l; j++){
			fprintf(root_file, "%s ", cols[root[i+l*j]]);
		}
		fprintf(root_file, "\n");
	}
	fclose(root_file);

	return 0;
}











