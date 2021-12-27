#define NSPECIES 6
#define NPARAM 12
#define NREACT 14



__constant__ int smatrix[]={
    1.0,    -1.0,    -1.0,    1.0,    0.0,    0.0,
    0.0,    0.0,    1.0,    -1.0,    0.0,    0.0,
    -1.0,    1.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    1.0,    -1.0,    0.0,    0.0,
    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,
    -1.0,    0.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    -1.0,    0.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    -1.0,    0.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    -1.0,    0.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    -1.0,    0.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,
    0.0,    0.0,    0.0,    0.0,    0.0,    -1.0,
};


__device__ void hazards(int *y, float *h, float t, int tid){

    h[0] = tex2D(param_tex,4,tid)*y[1]*y[2];
    h[1] = tex2D(param_tex,5,tid)*y[0]*y[3];
    h[2] = tex2D(param_tex,11,tid)*y[0];
    h[3] = tex2D(param_tex,6,tid)*y[4]*y[3];
    h[4] = tex2D(param_tex,8,tid);
    h[5] = tex2D(param_tex,7,tid);
    h[6] = tex2D(param_tex,9,tid)*y[0];
    h[7] = tex2D(param_tex,9,tid)*y[1];
    h[8] = tex2D(param_tex,9,tid)*y[2];
    h[9] = tex2D(param_tex,9,tid)*y[3];
    h[10] = tex2D(param_tex,10,tid)*__powf(y[3],tex2D(param_tex,2,tid))/(__powf(tex2D(param_tex,3,tid),tex2D(param_tex,2,tid))+__powf(y[3],tex2D(param_tex,2,tid)));
    h[11] = tex2D(param_tex,9,tid)*y[4];
    h[12] = tex2D(param_tex,1,tid)*__powf(y[3],tex2D(param_tex,2,tid))/(__powf(tex2D(param_tex,3,tid),tex2D(param_tex,2,tid))+__powf(y[3],tex2D(param_tex,2,tid)));
    h[13] = tex2D(param_tex,9,tid)*y[5];

}

