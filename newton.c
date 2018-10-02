#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"
#include "math.h"
#include "unistd.h"
#include "pthread.h"


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

Cpx cpx_inv(Cpx x)
{
	x = cpx_div_real(x, cpx_magnitude_sq(x));
	x.im = -x.im;
	return x;
}

void print_cpx(Cpx c){
	printf("%f%+fi\n",c.re,c.im);
}


double inv_degree;
Cpx newton_iteration(int degree, Cpx x)
{
	// calculates x - x^d/(dx^(d-1))
	//            = x - 1/d * (x-1/x^(d-1))
	//            = x * (1-1/d) +  1/((d-1)x^(d-1))

	Cpx inv_dfx = cpx_mul_real(cpx_inv(cpx_pow_unrolled(x, degree-1)),inv_degree);
	return cpx_add(cpx_mul_real(x,1-inv_degree), inv_dfx);
}

double cpx_distance_sq(Cpx a, Cpx b){
	return cpx_magnitude_sq(cpx_sub(a,b));
}

int degree = 7;
int l = 4000;
int num_threads = 1;
int *num_its;
int *root;
Cpx *correct_roots;
int num_threads;



void *process(void *data)
{
	int thread_index = *(int *)data;
	int64_t px_start = (l*(int64_t)l*thread_index) / (int64_t)num_threads;
	
	int64_t px_end =  ((int64_t)l*(int64_t)l*((int64_t)thread_index+1)) / (int64_t)num_threads;

	int64_t px_x = px_start % l;
	int64_t px_y = px_start / l;
	for(int64_t px = px_start; px< px_end;px++){
		Cpx x = (Cpx){.re = 4*(px_x / (double)l-0.5), 4*(px_y/(double)l - 0.5)};

		for(int its = 0;; its++){
			int done = 0;
			for(int j=0;j<degree;j++) {

				if(cpx_distance_sq(x, correct_roots[j]) < 1e-6){
					root[px] = j+1;
					num_its[px] = its+1;
					done = 1;
					break;
				}
			}
			if(done) break;
			if(fabs(x.re) > 1e10 || fabs(x.im) > 1e10 || cpx_magnitude_sq(x) < 1e-6 ){
				root[px] = degree+1;
				num_its[px] = its+1;
				break;
			}
			x = newton_iteration(degree,x);
		}

		if(++px_x == l){
			px_x = 0;
			++px_y; 
		}
	}
	return NULL;
}


int main(int argc, char *argv[])
{

	int opt;
	while((opt = getopt(argc, argv, "t:l:")) != -1)
	{
		switch(opt)
		{
		 case 'l':
		 	l = atoi(optarg);
		 	break;
		 case 't':
		 	num_threads = atoi(optarg);
		 	break;
		 break;
		}
	}

 	degree = atoi(argv[argc-1]);
	inv_degree = 1.0/degree;

	const double pi = M_PI;

	num_its 		= malloc(sizeof(int)*l*l); 
	root    		= malloc(sizeof(int)*l*l);
	correct_roots 	= malloc(sizeof(Cpx)*degree);
	
	for(int i = 0; i<degree;i++) {
		correct_roots[i] = (Cpx){cos(2.0*pi*i/(double)degree),sin(2.0*pi*i/(double)degree)};
	}



	pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
	int *thread_args = malloc(num_threads * sizeof(int));
	
  	for(int i = 0;i<num_threads;i++){
  		thread_args[i] = i;
		pthread_create(&threads[i],NULL,process,&thread_args[i]);
  	}

	for(int i=0;i<num_threads;i++){
    	int result_code=pthread_join(threads[i],NULL);
    	if(result_code) printf("error joining");
  	}

	char *cols[] = {"0 0 0 ", "255 0 0 ", "0 255 0 ", "0 0 255 ","255 0 255 ","0 255 255 ","255 255 0 ",  "255 255 255 ", "200 100 150 ", "200 150 100 ", "150 200 100 "};

	char buffer[100];
	sprintf(buffer, "newton_attractors_x%d.ppm",degree);
	FILE *root_file= fopen(buffer, "wb");
	sprintf(buffer, "newton_convergence_x%d.ppm",degree);
	FILE *conv_file= fopen(buffer, "wb");
	fprintf(root_file, "P3\n%d %d\n%d\n", l, l, 255);
	fprintf(conv_file, "P3\n%d %d\n%d\n", l, l, 255);

	printf("done calculating\n");

	char gray_scale[255][20];
	for(int i = 0; i<256;i++) sprintf(gray_scale[i], "%3d %3d %3d ", i,i,i);

	for(int64_t i = 0; i< l*l; i++){
		fwrite(cols[root[i]], 1, strlen(cols[root[i]]), root_file);
		int x = num_its[i]*2;
		x = x > 255? 255:x ;
		fwrite(gray_scale[x], 1, 12, conv_file);
	} 
	fclose(root_file);
	fclose(conv_file);



	return 0;
}











