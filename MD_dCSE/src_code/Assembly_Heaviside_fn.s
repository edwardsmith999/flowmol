//----------------------------------------------------------------------------
.intel_syntax noprefix
.text

//-----------------------------------------------------------------------------
// this heaviside function generates its own register constants
// double  heaviside_a1 (double arg);
.globl heaviside_a1

heaviside_a1:
   xorpd   xmm1,xmm1                # Generate constant with XOR so xmm1=0.0
   cmplesd xmm1,xmm0                # Compare xmm0 (arg) to see if it is greater than xmm1 (zero). Store in xmm1
   mov     rax,0x3ff0000000000000   # Load general purpose register rax with value 1.0
   movq    xmm0,rax                 # then move quadword so xmm0=1.0 (can't set xmm0=1.0 directly)
   andpd   xmm0,xmm1                # As compare returns all zeros/ones, AND with xmm0=1 to set return xmm0=0.0/1.0
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
