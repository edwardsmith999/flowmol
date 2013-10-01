//----------------------------------------------------------------------------
.intel_syntax noprefix
.text

//-----------------------------------------------------------------------------
// this heaviside function generates its own register constants
// double  heaviside_a1 (double arg);
.globl heaviside_a1

heaviside_a1:
   mov     rax,0x3ff0000000000000
   xorpd   xmm1,xmm1                # xmm1: constant 0.0
   cmplesd xmm1,xmm0                # xmm1: mask (all Fs or all 0s)
   movq    xmm0,rax                 # xmm0: constant 1.0
   andpd   xmm0,xmm1
   retq

//-----------------------------------------------------------------------------
// this heaviside function uses register constants passed from caller
// double  heaviside_a2 (double arg, double const0, double const1);
.globl heaviside_a2

heaviside_a2:
   cmplesd xmm1,xmm0                # xmm1: mask (all Fs or all 0s)
   movsd   xmm0,xmm2                # xmm0: constant 1.0
   andpd   xmm0,xmm1
   retq

//-----------------------------------------------------------------------------
