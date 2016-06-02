/* Tempering */



"use strict";



/* create a SimTemp object */
function SimTemp(beta, lnz)
{
  this.n = beta.length;
  this.beta = beta;
  this.hist = newarr(this.n);
  if ( lnz ) {
    this.lnz = lnz;
  } else {
    this.lnz = newarr(this.n);
  }
  this.usum = newarr(this.n);
  this.u2sum = newarr(this.n);
  this.lnzref = newarr(this.n);
  this.uref = newarr(this.n);
  this.cvref = newarr(this.n);
}



SimTemp.prototype.jump = function(itp, ep, type, ave)
{
  var jtp, n = this.n, dlnz;

  if ( type == 1 ) {
    // choose any temperature other than `itp`
    jtp = ( itp + 1 + Math.floor(rand01() * (n - 1)) ) % n;
  } else {
    // choose one of the neighbors
    if ( rand01() > 0.5 ) {
      jtp = itp + 1;
      if ( jtp >= n ) return itp;
    } else {
      jtp = itp - 1;
      if ( jtp < 0 ) return itp;
    }
  }

  var dbeta = this.beta[itp] - this.beta[jtp];
  if ( !ave ) {
    dlnz = this.lnz[itp] - this.lnz[jtp];
  } else {
    if ( ave === 1 ) {
      var eav1 = this.usum[itp] / this.hist[itp];
      var eav2 = (this.hist[jtp] > 0) ? this.usum[jtp] / this.hist[jtp] : eav1;
      dlnz = -dbeta * 0.5 * (eav1 + eav2);
    } else {
      var eav = (this.usum[itp] + this.usum[jtp])
              / (this.hist[itp] + this.hist[jtp]);
      dlnz = -dbeta * eav;
    }
  }

  // compute the acceptance probability
  var x = dbeta * ep + dlnz;
  if ( x > 0 ) {
    return jtp;
  } else {
    var r = rand01();
    return r < Math.exp(x) ? jtp : itp;
  }
}



SimTemp.prototype.add = function(itp, ep)
{
  this.hist[itp] += 1;
  this.usum[itp] += ep;
  this.u2sum[itp] += ep * ep;
}
