#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "transit.h"
 
void dbl_free(double *ptr){
  /* 
      Called by python to free a double pointer
  */ 
  free(ptr);
} 
 
double modulus(double x, double y) {
  /*
      The arithmetic modulus, x mod y
  */
  return x - y * floor(x / y);
} 
 
double ellec(double k) {
  /*
      Computes polynomial approximation for the complete elliptic
      integral of the second kind (Hasting's approximation).
  */
  double m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2;
  m1 = 1.-k*k;
  a1 = 0.44325141463;
  a2 = 0.06260601220;
  a3 = 0.04757383546;
  a4 = 0.01736506451;
  b1 = 0.24998368310;
  b2 = 0.09200180037;
  b3 = 0.04069697526;
  b4 = 0.00526449639;
  ee1 = 1.+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
  ee2 = m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1./m1);
  return ee1 + ee2;
}

double ellk(double k) {
  /*
      Computes polynomial approximation for the complete elliptic
      integral of the first kind (Hasting's approximation).
  */
  double a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ek1,ek2,m1;
  m1 = 1.-k*k;
  a0 = 1.38629436112;
  a1 = 0.09666344259;
  a2 = 0.03590092383;
  a3 = 0.03742563713;
  a4 = 0.01451196212;
  b0 = 0.5;
  b1 = 0.12498593597;
  b2 = 0.06880248576;
  b3 = 0.03328355346;
  b4 = 0.00441787012;
  ek1 = a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
  ek2 = (b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1);
  return ek1 - ek2;
}
 
double rc(double x, double y, int *err) { 
  /* 
        Carlson's degenerate elliptic integral
        (C) Copr. 1986-92 Numerical Recipes Software "U,6VV'.
  */ 
  double alamb,ave,s,w,xt,yt;   
  *err = ERR_NONE;
  if (x < 0.0 || y == 0.0 || (x+fabs(y)) < RC_TINY || (x+fabs(y)) > RC_BIG ||   
                             (y<-RC_COMP1 && x > 0.0 && x < RC_COMP2)) { 
    *err = ERR_RC;
    return 0.;   
  }
  if (y > 0.0) {   
    xt=x;   
    yt=y;   
    w=1.0;   
  } else {   
    xt=x-y;   
    yt = -y;   
    w=sqrt(x)/sqrt(xt);   
  }   
  do {   
    alamb=2.0*sqrt(xt)*sqrt(yt)+yt;   
    xt=0.25*(xt+alamb);   
    yt=0.25*(yt+alamb);   
    ave=RC_THIRD*(xt+yt+yt);   
    s=(yt-ave)/ave;   
  } while (fabs(s) > RC_ERRTOL);   
  return w*(1.0+s*s*(RC_C1+s*(RC_C2+s*(RC_C3+s*RC_C4))))/sqrt(ave);   
}  

double rj(double x, double y, double z, double p, int *err) {    
  /* 
      Carlson's elliptic integral of the third kind
      (C) Copr. 1986-92 Numerical Recipes Software "U,6VV'. 
  */
  double a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,   
         ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;   
  *err = ERR_NONE;
  if (DMIN(DMIN(x,y),z) < 0.0 || DMIN(DMIN(x+y,x+z),DMIN(y+z,fabs(p))) < RJ_TINY   
                              || DMAX(DMAX(x,y),DMAX(z,fabs(p))) > RJ_BIG) {
    *err = ERR_RJ;
    return 0.;    
  }  
  sum=0.0;   
  fac=1.0;   
  if (p > 0.0) {   
    xt=x;   
    yt=y;   
    zt=z;   
    pt=p;   
  } else {   
    xt=DMIN(DMIN(x,y),z);   
    zt=DMAX(DMAX(x,y),z);   
    yt=x+y+z-xt-zt;   
    a=1.0/(yt-p);   
    b=a*(zt-yt)*(yt-xt);   
    pt=yt+b;   
    rho=xt*zt/yt;   
    tau=p*pt/yt;   
    rcx=rc(rho,tau,err);  
    if (*err != ERR_NONE) return 0.;
  }   
  do {   
    sqrtx=sqrt(xt);   
    sqrty=sqrt(yt);   
    sqrtz=sqrt(zt);   
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;   
    alpha=SQR(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);   
    beta=pt*SQR(pt+alamb);   
    sum += fac*rc(alpha,beta,err);  
    if (*err !=0) return 0.; 
    fac=0.25*fac;   
    xt=0.25*(xt+alamb);   
    yt=0.25*(yt+alamb);   
    zt=0.25*(zt+alamb);   
    pt=0.25*(pt+alamb);   
    ave=0.2*(xt+yt+zt+pt+pt);   
    delx=(ave-xt)/ave;   
    dely=(ave-yt)/ave;   
    delz=(ave-zt)/ave;   
    delp=(ave-pt)/ave;   
  } while (DMAX(DMAX(fabs(delx),fabs(dely)),   
    DMAX(fabs(delz),fabs(delp))) > RJ_ERRTOL);   
  ea=delx*(dely+delz)+dely*delz;   
  eb=delx*dely*delz;   
  ec=delp*delp;   
  ed=ea-3.0*ec;   
  ee=eb+2.0*delp*(ea-ec);   
  ans=3.0*sum+fac*(1.0+ed*(-RJ_C1+RJ_C5*ed-RJ_C6*ee)+eb*
      (RJ_C7+delp*(-RJ_C8+delp*RJ_C4))+delp*ea*(RJ_C2-delp*RJ_C3)-RJ_C2*delp*ec)/
      (ave*sqrt(ave));   
  if (p <= 0.0) ans=a*(b*ans+3.0*(rcx-rf(xt,yt,zt,err)));
  if (*err != ERR_NONE) return 0.;   
  return ans;   
}   
  
double rf(double x, double y, double z, int *err) {  
  /* 
      Carlson's elliptic integral of the first kind
      (C) Copr. 1986-92 Numerical Recipes Software "U,6VV'.
  */ 
  double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;   
  *err = ERR_NONE;
  if (DMIN(DMIN(x,y),z) < 0.0 || DMIN(DMIN(x+y,x+z),y+z) < RF_TINY ||   
      DMAX(DMAX(x,y),z) > RF_BIG) {  
    *err = ERR_RF;
    return 0.;  
  }
  xt=x;   
  yt=y;   
  zt=z;   
  do {   
    sqrtx=sqrt(xt);   
    sqrty=sqrt(yt);   
    sqrtz=sqrt(zt);   
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;   
    xt=0.25*(xt+alamb);   
    yt=0.25*(yt+alamb);   
    zt=0.25*(zt+alamb);   
    ave=RF_THIRD*(xt+yt+zt);   
    delx=(ave-xt)/ave;   
    dely=(ave-yt)/ave;   
    delz=(ave-zt)/ave;   
  } while (DMAX(DMAX(fabs(delx),fabs(dely)),fabs(delz)) > RF_ERRTOL);   
  e2=delx*dely-delz*delz;   
  e3=delx*dely*delz;   
  return (1.0+(RF_C1*e2-RF_C2-RF_C3*e3)*e2+RF_C4*e3)/sqrt(ave);   
}   

double sgn(double x) {
  /* 
      Returns the sign of x
  */ 
  return (x > 0) - (x < 0);
}
 
double TrueAnomaly(double E, double ecc) {
  /* 
      The true anomaly as a function of the eccentric anomaly E
      and the eccentricity ecc
  */ 
  if (ecc == 0.) return E;
  else return 2. * atan2(pow(1. + ecc, 0.5) * sin(E / 2.), 
                         pow(1. - ecc, 0.5) * cos(E / 2.));
}

double EccentricAnomalyFast(double dMeanA, double dEcc, double tol, int maxiter) {
  /* 
      Adapted from Rory Barnes, based on Murray & Dermott 
  */
  
  double dEccA;
  double di_1, di_2, di_3 = 1.0, fi, fi_1, fi_2, fi_3;
  double lo = -2 * PI;
  double up = 2 * PI;
  double next;
  int iter;
  
  if (dEcc == 0.) return dMeanA;                                                      // The trivial circular case
  dEccA = dMeanA + sgn(sin(dMeanA))*0.85*dEcc;

  for (iter = 1; iter <= maxiter; iter++) {
    fi = dEccA - dEcc*sin(dEccA) - dMeanA;
    if (fi > 0)
      up = dEccA;
    else
      lo = dEccA;
    fi_1 = 1.0 - dEcc*cos(dEccA);
    fi_2 = dEcc*sin(dEccA);
    fi_3 = dEcc*cos(dEccA);
    di_1 = -fi / fi_1;
    di_2 = -fi / (fi_1 + 0.5*di_1*fi_2);
    di_3 = -fi / (fi_1 + 0.5*di_2*fi_2 + 1./6.*di_2*di_2*fi_3);
    next = dEccA + di_3;
    
    if (fabs(dEccA - next) < tol) 
      break;
      
    if ((next > lo) && (next < up)) 
      dEccA = next;
    else 
      dEccA = (lo + up) / 2.;
      
    if ((fabs(dEccA - lo) < tol) || (fabs(dEccA - up) < tol))
      break;
  }
  
  if (iter >= maxiter) 
    return -1.;                                                                       // Solver didn't converge
  else
    return dEccA;
}

double EccentricAnomaly(double M, double e, double tol, int maxiter) {
  /*  
      A simpler version of the Kepler solver, borrowed from
      https://github.com/lkreidberg/batman/blob/master/c_src/_rsky.c
  */
  
  double E = M, eps = tol;                                                            // Kreidberg: eps = 1.0e-7;
  
  if (e == 0.) return M;                                                              // The trivial circular case
  
	while(fabs(E - e*sin(E) - M) > eps) E = E - (E - e*sin(E) - M)/(1.0 - e*cos(E));
	return E;
	
}

int Compute(TRANSIT *transit, LIMBDARK *limbdark, SETTINGS *settings, ARRAYS *arr){
  /*
      Compute the transit model
  */    
  double au, bu, u1, u2, c1, c2, c3, c4;
  double omega, per, RpRs, aRs, inc, w, ecc, fi, tperi0, t;
  double dt, tmp;
  double x1, x2, x3, x4, kap1, kap0, lambdae, lambdad, lam, q, Kk, Ek, n, Pk, etad;
  int i, s;
  int np = 0, nm = 0, npctr = 0, nmctr = 0;
  int iErr = ERR_NONE;

  arr->time = malloc(settings->maxpts*sizeof(double)); 
  arr->flux = malloc(settings->maxpts*sizeof(double)); 
  arr->M = malloc(settings->maxpts*sizeof(double)); 
  arr->E = malloc(settings->maxpts*sizeof(double)); 
  arr->f = malloc(settings->maxpts*sizeof(double)); 
  arr->r = malloc(settings->maxpts*sizeof(double)); 
  arr->x = malloc(settings->maxpts*sizeof(double)); 
  arr->y = malloc(settings->maxpts*sizeof(double)); 
  arr->z = malloc(settings->maxpts*sizeof(double)); 
  arr->b = malloc(settings->maxpts*sizeof(double)); 
  arr->calloc = 1;

  if (settings->exppts % 2) return ERR_EXP_PTS;                                       // Verify user input: Must be even!

  if (limbdark->ldmodel == QUADRATIC) {                                               // Verify user input: Limb darkening model
    u1 = limbdark->u1;
    u2 = limbdark->u2;
    if (isnan(u1) || isnan(u2)) return ERR_LD;
  } else if (limbdark->ldmodel == KIPPING) {
    au = sqrt(limbdark->q1);
    bu = 2*limbdark->q2;
    u1 = au*bu;
    u2 = au*(1 - bu);    
    if (isnan(u1) || isnan(u2)) return ERR_LD;
  } else if (limbdark->ldmodel == NONLINEAR) {
    // TODO: Implement this!
    return ERR_NOT_IMPLEMENTED;
  } else {
    return ERR_NOT_IMPLEMENTED;
  }
  
  per = transit->per;                                                                 // Orbital period
  if (!(per > 0.)) return ERR_PER;
  
  RpRs = transit->RpRs;                                                               // Planet radius in units of stellar radius
  if (!((RpRs > 0.) && (RpRs < 1.))) return ERR_RADIUS;
  
  if (isnan(transit->MpMs)) transit->MpMs = 0.;                                       // We'll assume the secondary is massless
  
  if (isnan(transit->rhos)) {                                                         // Stellar density
    if (isnan(transit->aRs)) return ERR_RHOS_ARS;
    else aRs = transit->aRs;
  } else {
    if (transit->rhos <= 0.) return ERR_RHOS;
    aRs = pow(((G * transit->rhos * (1. + transit->MpMs) * 
          pow(per * DAYSEC, 2)) / (3. * PI)), 1./3.);                                 // Semi-major axis in units of stellar radius
    transit->aRs = aRs;
  }
  
  inc = acos(transit->bcirc / aRs);                                                   // Orbital inclination
  
  if (isnan(transit->esw) || isnan(transit->ecw)) {                                   // Eccentricity and longitude of pericenter
    if (isnan(transit->ecc)) return ERR_ECC_W;
    if ((transit->ecc != 0) && isnan(transit->w)) 
      return ERR_ECC_W;
    else if (transit->ecc == 0)
      transit->w = 0;
    else if (isnan(transit->w))
      return ERR_ECC_W;                           
    if ((transit->ecc < 0) || (transit->ecc >= 1)) return ERR_ECC_W;
    if ((transit->w < 0) || (transit->w >= 2 * PI)) return ERR_ECC_W;
    w = transit->w;
    ecc = transit->ecc;
  } else {
    w = atan2(transit->esw, transit->ecw);
    ecc = sqrt(transit->esw * transit->esw + transit->ecw * transit->ecw);
    if ((ecc < 0.) || (ecc >= 1.)) return ERR_BAD_ECC;
    transit->ecc = ecc;
    transit->w = w;
  }
  
  // HACK: My definition of omega in the equations below is apparently
  // off by 180 degrees from Laura Kreidberg's in BATMAN. This isn't elegant,
  // but the two models agree now that I added the following line:
  w = w - PI;
  
  fi = (3. * PI / 2.) - w;                                                            // True anomaly at transit center (Shields et al. 2015)
  tperi0 = per * sqrt(1. - ecc * ecc) / (2. * PI) * (ecc * sin(fi) / 
           (1. + ecc * cos(fi)) - 2. / sqrt(1. - ecc * ecc) * 
           atan2(sqrt(1. - ecc * ecc) * tan(fi/2.), 1. + ecc));                       // Time of pericenter passage (Shields et al. 2015)

  omega = 1. - u1/3. - u2/6.;                                                         // See Mandel and Agol (2002)
  dt = settings->exptime / settings->exppts;                                          // The time step
  
  for (s = -1; s <= 1; s+=2) {                                                        // Sign: -1 or +1
    t = 0.;
    for (i = settings->maxpts/2; ((i < settings->maxpts) && (i >= 0)) ; i+=s) {                           // Loop over all points. Start from transit center and go left, then right
         
      /*
      --- ORBITAL SOLUTION ---
      */
      
      arr->time[i] = t;
      arr->M[i] = 2. * PI / per * (arr->time[i] - tperi0);                            // Mean anomaly
      if (settings->kepsolver == MDFAST)
        arr->E[i] = EccentricAnomalyFast(arr->M[i], ecc, settings->keptol, 
                                         settings->maxkepiter);                       // Eccentric anomaly
      else
        arr->E[i] = EccentricAnomaly(arr->M[i], ecc, settings->keptol, 
                                     settings->maxkepiter);
      if (arr->E[i] == -1) return ERR_KEPLER;
      arr->f[i] = TrueAnomaly(arr->E[i], ecc);                                        // True anomaly
      arr->r[i] = aRs * (1. - ecc * ecc)/(1. + ecc * cos(arr->f[i]));                 // Star-planet separation in units of stellar radius
      if (arr->r[i] - RpRs < 1.) return ERR_STAR_CROSS;                               // Star-crossing orbit!
      arr->b[i] = arr->r[i] * sqrt(1. - pow(sin(w + arr->f[i]) * sin(inc), 2.));      // Instantaneous impact parameter                                   
      arr->x[i] = arr->r[i] * cos(w + arr->f[i]);                                     // Cartesian sky-projected coordinates
      arr->z[i] = arr->r[i] * sin(w + arr->f[i]);
      if (arr->b[i] * arr->b[i] - arr->x[i] * arr->x[i] < 1.e-10) 
        arr->y[i] = 0.;                                                               // Prevent numerical errors
      else {
        tmp = modulus(arr->f[i] + w, 2 * PI);                                         // TODO: Verify this modulus
        arr->y[i] = sqrt(arr->b[i] * arr->b[i] - arr->x[i] * arr->x[i]);
        if (!((0 < tmp) && (tmp < PI))) arr->y[i] *= -1;
      }
      t += s*dt;                                                                      // Increment the time
      
      if (!settings->fullorbit) {                                                     // We're only calculating stuff during transit
        if ((arr->b[i] > 1. + RpRs) || (arr->z[i] > 0)) {                             // Check if we're done transiting, or if it's a secondary eclipse (which we ignore)
          arr->flux[i] = 1.;                                                          // That's easy!
          if (s == -1) {
            nm = i;                                                                   // We're going to truncate the array at this index on the left
            nmctr++;                                                                  // We want to add exppts/2 points on each side of the transit
            if (nmctr == settings->exppts/2) break;                                   // since we'll eventually need those points for binning
          } 
          else if (s == 1) {
            np = i;                                                                   // We're going to truncate the array at this index on the right
            npctr++;
            if (npctr == settings->exppts/2 + 1) break;                               // Note the + 1 on this line to ensure the same number of points on the left and on the right
          }
          continue;
        }
      } else {
        if (fabs(t) >= per/2.) {                                                      // We're going to calculate the full orbit, but we know the flux is 1.
          if (s == -1) nm = i + 1;
          else if (s == 1) np = i - 1;
          break;
        } else {
          if ((arr->b[i] > 1. + RpRs) || (arr->z[i] > 0)) {                           // Ignoring secondary eclipse
            arr->flux[i] = 1.;
            continue;
          }
        }
      
      }

      /*
      --- TRANSIT FLUX ---
      */
      
      x1 = pow(RpRs - arr->b[i], 2.);                                                 // Set up some quantities to compute the transit flux
      x2 = pow(RpRs + arr->b[i], 2.);                                                 // The following is adapted from Eric Agol's fortran routines
      x3 = RpRs * RpRs - arr->b[i]*arr->b[i];
      x4 = RpRs * RpRs + arr->b[i]*arr->b[i];
      
      // 1. Compute lambdae
      if (RpRs >= 1. && arr->b[i] <= RpRs - 1.) {                                     // [ONE] Occulting object completely occults source
        lambdae=1.;
      } else if (arr->b[i] > 1. - RpRs) {                                             // [TWO] Occultor is crossing the limb. Equation (26)
        kap1 = acos(fmin((1. - x3) / 2. / arr->b[i], 1.));
        kap0 = acos(fmin((x4 - 1.) / 2 / RpRs / arr->b[i], 1.));
        lambdae = RpRs * RpRs * kap0 + kap1;
        lambdae -= 0.5*sqrt(fmax(4. * arr->b[i] * arr->b[i] - pow(1. - x3, 2.), 0.));
        lambdae /= PI;
      } else if (arr->b[i] <= 1. - RpRs) {                                            // [THREE] Occultor is crossing the star
        lambdae = RpRs * RpRs;
      }
      
      // 2. Compute lambdad and etad
      if (RpRs >= 1. && arr->b[i] <= RpRs - 1.) {                                     // [ONE] Occulting object completely occults source
        lambdad=1.;
        etad=1.;
      } else if ((arr->b[i] > 0.5 + fabs(RpRs - 0.5) && arr->b[i] < 1. + RpRs) || 
                 (RpRs > 0.5 && arr->b[i] > fabs(1. - RpRs) * 
                 1.0001 && arr->b[i] < RpRs)) {                                       // [TWO] The occultor partly occults the star and crosses the limb
        n = 1./x1 - 1.;
        
        if (1. + n > RJ_BIG){
          // When the impact parameter approaches RpRs, x1 tends to zero and
          // n tends to infinity. The old approach was to set n = RJ_BIG - 1,
          // but this introduces its own set of issues. Here instead we use the
          // equations in Table 3, Case V.
          // NOTE: TODO: Verify this section. Not yet tested.
          if (RpRs == 0.5) {
            lambdad = 1. / 3. - 4. / PI / 9.;
            etad = 3. / 32.;
          } else {
            lam = 0.5 * PI;
            q = 0.5 / RpRs;
            Kk = ellk(q);
            Ek = ellec(q);
            lambdad = 1. / 3. + 16. * RpRs / 9. / PI * (2. * RpRs * RpRs - 1.) * Ek - 
                      (32. * pow(RpRs, 4) - 20. * RpRs * RpRs + 3.) / 9. / PI / 
                      RpRs * Kk;
            etad = 1. / 2. / PI * (kap1 + RpRs * RpRs * (RpRs * RpRs + 2. * 
                   arr->b[i] * arr->b[i]) * kap0 - (1. + 5. * RpRs * RpRs + 
                   arr->b[i] * arr->b[i]) / 4. * sqrt((1. - x1) * (x2 - 1.)));
          }
        } else {
          // Business as usual.
          lam = 0.5 * PI;
          q = sqrt((1. - x1)/ 4. / arr->b[i] / RpRs);
          Kk = ellk(q);
          Ek = ellec(q);
          Pk = Kk - n / 3. * rj(0., 1. - q * q, 1., 1. + n, &iErr);
          if (iErr != ERR_NONE) return iErr;
          lambdad = 1. / 9. / PI / sqrt(RpRs * arr->b[i]) * (((1. - x2) * 
                    (2. * x2 + x1 - 3.) - 3. * x3 * (x2 - 2.)) * Kk + 4. * 
                    RpRs * arr->b[i] * ( arr->b[i] * arr->b[i] + 7. * RpRs * 
                    RpRs - 4.) * Ek - 3. * x3 / x1 * Pk);                             // Equation (34), lambda_1
          if (arr->b[i] < RpRs) lambdad += 2./3.;
          etad = 1. / 2. / PI * (kap1 + RpRs * RpRs * 
                (RpRs * RpRs + 2. * arr->b[i] * arr->b[i]) * kap0 - 
                (1. + 5. * RpRs * RpRs + arr->b[i] * arr->b[i]) / 4. * 
                sqrt((1. - x1) * (x2 - 1.)));                                         // Equation (34), eta_1
        }
      } else if (RpRs <= 1. && arr->b[i] <= (1. - RpRs) * 1.0001) {                   // [THREE] Occultor is crossing the star
          n = x2 / x1 - 1.;
          
          if (1. + n > RJ_BIG) {
            // When the impact parameter approaches RpRs, x1 tends to zero and
            // n tends to infinity. The old approach was to set n = RJ_BIG - 1,
            // but this introduces its own set of issues. Here instead we use the
            // equations in Table 3, Case VI.
            lam = 0.5 * PI;
            q = 2. * RpRs;
            Kk = ellk(q);
            Ek = ellec(q);
            lambdad = 1. / 3. + 2. / 9. / PI * (4. * (2. * RpRs * RpRs - 1.) * Ek + 
                     (1. - 4. * RpRs * RpRs) * Kk);
            etad = RpRs * RpRs / 2. * (RpRs * RpRs + 2. * arr->b[i] * arr->b[i]);
          } else {
            // Business as usual.
            lam = 0.5 * PI;
            q = sqrt((x2 - x1) / (1. - x1));
            Kk = ellk(q);
            Ek = ellec(q);
            Pk = Kk - n / 3. * rj(0., 1. - q * q, 1., 1. + n, &iErr);
            if (iErr != ERR_NONE) return iErr;
            lambdad = 2. / 9. / PI / sqrt(1. - x1) * ((1. - 5. * arr->b[i] * 
                      arr->b[i] + RpRs * RpRs + x3 * x3) * Kk + (1. - x1) * (arr->b[i] 
                      * arr->b[i] + 7. * RpRs * RpRs - 4.) * Ek - 3. * x3 / x1 * Pk); // Equation (34), lambda_2   
            if (arr->b[i] < RpRs) lambdad += 2./3.;
            if (fabs(RpRs + arr->b[i] - 1.) <= 1.e-4)
              lambdad = 2. / 3. / PI * acos(1. - 2. * RpRs) - 4. / 9. / PI * 
                      sqrt(RpRs * (1. - RpRs)) * (3. + 2. * RpRs - 8. * RpRs * RpRs);
            etad = RpRs * RpRs / 2. * (RpRs * RpRs + 2. * arr->b[i] * arr->b[i]);     // Equation (34), eta_2
          }
      }
      
      arr->flux[i] = 1. - ((1. - u1 - 2. * u2) * lambdae + (u1 + 2. * u2) * 
                      lambdad + u2 * etad) / omega;                                   // Finally, the transit flux (baseline = 1.)

    }
  }
  
  if ((nm == 0) || (np == 0)) return ERR_MAX_PTS;                                     // We didn't reach the edge of the transit within settings->maxpts
  if ((nm >= settings->maxpts/2 - settings->exppts/2 - 1) && 
      (np <= settings->maxpts/2 + settings->exppts/2 +  1)) 
    return ERR_NO_TRANSIT;                                                            // There's no transit!
  arr->nstart = nm;                                                                   // first index
  arr->nend = np + 1;                                                                 // one plus last index
  settings->computed = 1;                                                             // Set the flag
	return iErr;
}

int Bin(TRANSIT *transit, LIMBDARK *limbdark, SETTINGS *settings, ARRAYS *arr) {
  int iErr = ERR_NONE;
  int i, j, ep, nb, hx; 
  double sum;

  arr->bflx = malloc(settings->maxpts*sizeof(double)); 
  arr->balloc = 1;
  
  if (!settings->computed) return ERR_NOT_COMPUTED;                                   // Must compute first!
  ep = settings->exppts;                                                              // Shortcut for exppts
  hx = ep/2;                                                                          // The number of extra points on each side of the transit
  nb = ep + 1;                                                                        // Actual number of points in bin must be odd, but user doesn't need to know this!
  
  if (settings->binmethod == RIEMANN) {
    arr->bflx[arr->nstart] = (arr->flux[arr->nstart + hx] + ep) / nb;                 // Set the leftmost bin
  
    for (i = arr->nstart + 1; i < arr->nstart + hx + 1; i++)                          // For these guys, the left edge of the exposure window starts prior to where we've
      arr->bflx[i] = arr->bflx[i - 1] + (arr->flux[i + hx] - 1.) / nb;                // calculated flux values, but we know that the flux is all 1.0 out here
  
    for (i = arr->nstart + hx + 1; i < arr->nend - hx; i++)                           // Intelligent summation to compute bins
      arr->bflx[i] = arr->bflx[i - 1] + 
                    (arr->flux[i + hx] - arr->flux[i - 1 - hx]) / nb;
  
    for (i = arr->nend - hx; i < arr->nend; i++)                                      // Again, deal with edge effects
      arr->bflx[i] = arr->bflx[i - 1] + (1. - arr->flux[i - 1 - hx]) / nb;
  
  } else if (settings->binmethod == TRAPEZOID) {
    arr->bflx[arr->nstart] = 1. + 0.5 / ep * (arr->flux[arr->nstart + hx] - 1.);      // Set the leftmost bin

    for (i = arr->nstart + 1; i < arr->nstart + hx + 1; i++)
      arr->bflx[i] = arr->bflx[i - 1] + 1. / (2 * ep) * (arr->flux[i + hx] + 
                     arr->flux[i + hx - 1] - 2.);                                     
  
    for (i = arr->nstart + hx + 1; i < arr->nend - hx; i++)
      arr->bflx[i] = arr->bflx[i - 1] + 1. / (2 * ep) * (arr->flux[i + hx] + 
                     arr->flux[i + hx - 1] - arr->flux[i - hx] - 
                     arr->flux[i - hx -1]);                                           // We're essentially doing the same intelligent summation as above
  
    for (i = arr->nend - hx; i < arr->nend; i++)
      arr->bflx[i] = arr->bflx[i - 1] + 1. / (2 * ep) * (2. - 
                     arr->flux[i - hx] - arr->flux[i - hx -1]);
    
  } else if (settings->binmethod == -1) {                                             // DEBUG: This is the old trapezoid routine. Use for testing only
    arr->bflx[arr->nstart] = 1. + 0.5 / ep * (arr->flux[arr->nstart + hx] - 1.);      // Set the leftmost bin

    for (i = arr->nstart + 1; i < arr->nstart + hx + 1; i++) {
      sum = 0;
      for (j = i - hx + 1; j < i + hx; j++) {
        if (j >= 0) sum += arr->flux[j];
        else sum += 1.;
      }
      arr->bflx[i] = 1. / (nb - 1) * (0.5 * (1. + arr->flux[i + hx]) + sum);
    }
    
    for (i = arr->nstart + hx + 1; i < arr->nend - hx; i++) {
      sum = 0;
      for (j = i - hx + 1; j < i + hx; j++) sum += arr->flux[j];
      arr->bflx[i] = 1. / (nb - 1) * (0.5 * (arr->flux[i - hx] + 
                                      arr->flux[i + hx]) + sum);
    } 
  
    for (i = arr->nend - hx; i < arr->nend; i++) {
      sum = 0;
      for (j = i - hx + 1; j < i + hx; j++) {
        if (j < arr->nend - arr->nstart) sum += arr->flux[j];                         // TODO: Verify this line (are we off by 1?)
        else sum += 1.;
      }
      arr->bflx[i] = 1. / (nb - 1) * (0.5 * (arr->flux[i - hx] + 1.) + sum);
    }
  
  } else {
	  return ERR_NOT_IMPLEMENTED;
	}
  
  settings->binned = 1;                                                               // Set the flag
  return iErr;


}

int Interpolate(double *t, int ipts, int array, TRANSIT *transit, LIMBDARK *limbdark, SETTINGS *settings, ARRAYS *arr) {
  
  double f1, f0, t1, t0, ti;
  int i, j, nt;
  int iErr = ERR_NONE;
  double *f;
  double fill_value;

  if (!(transit->ntrans))
    if (isnan(transit->t0)) return ERR_T0;                                            // User didn't specify t0!
  
  arr->iarr = malloc(ipts*sizeof(double));                                            // The interpolated array 
  arr->ialloc = 1;
  
  if (!settings->computed) {
    iErr = Compute(transit, limbdark, settings, arr);                                 // Compute the raw transit model if necessary
    if (iErr != ERR_NONE) return iErr;
  } 
  if ((array == ARR_BFLX) && (!settings->binned)) {
    iErr = Bin(transit, limbdark, settings, arr);                                     // Bin the transit if necessary
    if (iErr != ERR_NONE) return iErr;
  }
  
  // Select which array to interpolate
  
  if (array == ARR_FLUX) {
    f = arr->flux;
    fill_value = 1.;
  } else if (array == ARR_BFLX) {
    f = arr->bflx;
    fill_value = 1.;
  } else if (array == ARR_M) {
    f = arr->M;
    fill_value = NAN;
  } else if (array == ARR_E) {
    f = arr->E;
    fill_value = NAN;
  } else if (array == ARR_F) {
    f = arr->f;
    fill_value = NAN;
  } else if (array == ARR_R) {
    f = arr->r;
    fill_value = NAN;
  } else if (array == ARR_X) {
    f = arr->x;
    fill_value = NAN;
  } else if (array == ARR_Y) {
    f = arr->y;
    fill_value = NAN;
  } else if (array == ARR_Z) {
    f = arr->z;
    fill_value = NAN;
  } else if (array == ARR_B) {
    f = arr->b;
    fill_value = NAN;
  } else
    return ERR_NOT_IMPLEMENTED;
  
  j = 0;                                                                              // The interpolation index
  nt = 0;                                                                             // The transit number
  transit->tN[transit->ntrans] = HUGE;                                                // A bit of a hack, but essential to get the last transit right below
    
  for (i = 0; i < ipts; i++) {
    
    if (!(transit->ntrans))
      ti = modulus(t[i]-transit->t0-transit->per/2., transit->per) - transit->per/2.; // Find the folded time, assuming strict periodicity
    else {
      for (; nt < transit->ntrans; nt++) {                                            // Find the folded time given all of the transit times
        if (fabs(t[i] - transit->tN[nt]) < fabs(t[i] - transit->tN[nt + 1])) {
          ti = t[i] - transit->tN[nt];
          break;
        }
      }
    }
    
    if ((ti < arr->time[arr->nstart]) || (ti >= arr->time[arr->nend-1])) {            // The case ti == arr->time[arr->nend-1] is pathological,
      arr->iarr[i] = fill_value;                                                      // but we're technically overestimating the flux slightly
      continue;                                                                       // in the zero-probability event that this does occur
    }
                                                                                              
    // Now we find [j, j + 1], the indices bounding the data point
    
    if (settings->intmethod == SMARTINT) {                                            // Increment j intelligently. NOTE: time array must be sorted!
      if (j > 0) j += settings->exppts * (t[i] - t[i - 1])/settings->exptime;         
      j = j % (arr->nend - arr->nstart);
      
      if (arr->time[arr->nstart + j + 1] <= ti) {                                     // We undershot; let's loop until we get the right index
        for (; j < arr->nend - arr->nstart - 1; j++) {
          if (arr->time[arr->nstart + j + 1] > ti) break;
        }
      } else {                                                                        // We either overshot or got it right
        for (; j >= 0; j--) {
          if (arr->time[arr->nstart + j] < ti) break;
        }
      }
    } else if (settings->intmethod == SLOWINT) {                                      // Brain-dead slow interpolation, useful if time array isn't sorted
      for (j = 0; j < arr->nend - arr->nstart - 1; j++) {
        if (arr->time[arr->nstart + j + 1] > ti) break;
      }
    } else {
      return ERR_NOT_IMPLEMENTED;
    }
    
    t0 = arr->time[arr->nstart + j];                                                  // Interpolation bounds
    t1 = arr->time[arr->nstart + j + 1];
    f0 = f[arr->nstart + j];
    f1 = f[arr->nstart + j + 1];
  
    arr->iarr[i] = f0 + (f1 - f0) * (ti - t0) / (t1 - t0);                            // A simple linear interpolation
    
  }
  
  arr->ipts = ipts;
  
  return iErr;

}