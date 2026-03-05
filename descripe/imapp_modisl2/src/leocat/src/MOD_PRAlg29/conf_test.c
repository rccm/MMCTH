/*
!C ********************************************************************
!Description:

  Float function conf_test.c

  Computes clear sky confidence value for an individual spectral test.
  Input single threshold value ('midpt') with associated confidence
  limits ('locut', 'hicut').  "Low" and "high" cutoffs refer to low
  or high confidence ends of an interval and not necessarily to 
  absolute values.  Routine calculates the confidence based on an "S"
  function.  One may change the shape of the function by changing
  'power' and/or 'midpt'.

!Input arguments:
  float val         parameter to be tested (BT, BTD, etc.)
  float locut       low confidence cutoff (< locut = confidence of 0.0)
  float hicut       high confidence cutoff (> hicut = confidence of 1.0)
  float power       S function power (1.0 = linear funtion)
  float midpt       midpoint of S function (= confidence of 0.5)

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


float conf_test(float val, float locut, float hicut, float power,
                float midpt)


{

/*  Declarations */

    int flipped;

    float confidence;
    float coeff;
    float alpha;
    float gamma;
    float beta;
    float range;
    float s1;
    float c;

/*  Initializations */

    coeff = powf(2.0, (power - 1.0));

/*  Find if high confidence cutoff is the numerically lowest or highest
    value. */
    if(hicut > locut) {
      gamma = hicut;
      alpha = locut;
      flipped = 0;
    } else {
      gamma = locut;
      alpha = hicut;
      flipped = 1;
    }

    beta = midpt;

/*  If value ('val') lies outside function range, confidence is either 
    0.0 or 1.0. */

    if(flipped == 0 && val > gamma) {
      c = 1.0; }

    else if(flipped == 0 && val < alpha) {
      c = 0.0; }

    else if(flipped == 1 && val > gamma) {
      c = 0.0; }

    else if(flipped == 1 && val < alpha) {
      c = 1.0; }

    else {

/*    Otherwise, value ('val') is within the function range and the
      confidence is between 0.0 and 1.0. */

      if(val <= beta) {

        range = 2.0 * (beta - alpha);
        s1 = (val - alpha) / range;
        if(flipped == 0) c = coeff * powf(s1, power);
        if(flipped == 1) c = 1.0 - (coeff * powf(s1, power)); }

      else {

        range = 2.0 * (beta - gamma);
        s1 = (val - gamma) / range;
        if(flipped == 0) c = 1.0 - (coeff * powf(s1, power));
        if(flipped) c = coeff * powf(s1, power);
 
      }

    }

/*  Force confidence values to be between 0 and 1. */
    if(c > 1.0) c = 1.0;
    if(c < 0.0) c = 0.0;

    confidence = c;

    return confidence;

}
