/*
!C ********************************************************************
!Description:

  Float function conf_test_dble.c

  Computes clear sky confidence value for an individual spectral test,
  but where low confidence values fall in the middle of the confidence
  interval and high confidence values are found at either end.
  Input two threshold values ('midpt') with associated confidence
  limits ('locut', 'hicut') for each.  "Low" and "high" cutoffs refer to
  low or high confidence ends of an interval and not necessarily to 
  absolute values.  Routine calculates the confidence based on an "S"
  function.  One may change the shape of the function by changing
  'power' and/or 'midpt'.

!Input arguments:
  float val        parameter to be tested (BT, BTD, etc.)
  float locut      low confidence cutoffs (< locut = confidence of 0.0)
  float hicut      high confidence cutoffs (> hicut = confidence of 1.0)
  float power      S function power (1.0 = linear funtion)
  float midpt      midpoints of S function (= confidence of 0.5)

!Output arguments:
  none

  float confidence   output clear sky confidence (0.0 - 1.0)
                     (returned through function call)

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <math.h>


float conf_test_dble(float val, float locut[2], float hicut[2],
                     float power, float midpt[2])


{

/*  Declarations */

    float confidence;
    float coeff;
    float alpha1;
    float alpha2;
    float gamma1;
    float gamma2;
    float beta1;
    float beta2;
    float range;
    float s1;
    float c;

/*  Initializations */

    coeff = powf(2.0, (power - 1.0));

    gamma1 = hicut[0];
    gamma2 = hicut[1];
    alpha1 = locut[0];
    alpha2 = locut[1];
    beta1 = midpt[0];
    beta2 = midpt[1];

/*  Find if interval between inner cutoffs passes or fails test. */
    if(alpha1 - gamma1 > 0.0) {

/*    Inner region fails test. */
/*    Check for value beyond function range. */
      if(val > alpha1 && val < alpha2) {

        c = 0.0; }

      else if(val < gamma1 || val > gamma2) {

        c = 1.0; }

      else if(val <= alpha1) {

/*      Value is within range of lower set of limits. */
        if(val >= beta1) {
          range = 2.0 * (beta1 - alpha1);
          s1 = (val - alpha1) / range;
          c = coeff * powf(s1, power); }
        else {
          range = 2.0 * (beta1 - gamma1);
          s1 = (val - gamma1) / range;
          c = 1.0 - (coeff * (powf(s1, power)));
        }
      }

      else {

/*      Value is within range of upper set of limits. */
        if(val <= beta2) {
          range = 2.0 * (beta2 - alpha2);
          s1 = (val - alpha2) / range;
          c = coeff * (powf(s1, power)); }
        else {
          range = 2.0 * (beta2 - gamma2);
          s1 = (val - gamma2) / range;
          c = 1.0 - (coeff * (powf(s1, power)));
        }
      }

    }

    else {

/*    Inner region passes test. */
/*    Check for value beyond function range. */

      if(val > gamma1 && val < gamma2) {

        c = 1.0; }

      else if(val < alpha1 || val > alpha2) {

        c = 0.0; }

      else if(val <= gamma1) {

/*      Value is within range of lower set of limits. */
        if(val <= beta1) {
          range = 2.0 * (beta1 - alpha1);
          s1 = (val - alpha1) / range;
          c = coeff * powf(s1, power); }
        else {
          range = 2.0 * (beta1 - gamma1);
          s1 = (val - gamma1) / range;
          c = 1.0 - (coeff * (powf(s1, power)));
        }
      }

      else {

/*      Value is within range of upper set of limits. */
        if(val >= beta2) {
          range = 2.0 * (beta2 - alpha2);
          s1 = (val - alpha2) / range;
          c = coeff * powf(s1, power); }
        else {
          range = 2.0 * (beta2 - gamma2);
          s1 = (val - gamma2) / range;
          c = 1.0 - (coeff * (powf(s1, power)));
        }
      }

    }
        
/*  Force confidence values to be between 0 and 1. */
    if(c > 1.0) c = 1.0;
    if(c < 0.0) c = 0.0;

    confidence = c;

    return confidence;

}
