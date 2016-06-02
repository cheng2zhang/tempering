


"use strict";



function transform(x, l)
{
  var n = x.length;
  var xyz = newarr2d(n, 3), xc = [l * 0.5, l * 0.5, l * 0.5], xi = [0, 0, 0];

  for ( var i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    vmxv(xyz[i], viewmat, xi);
    vinc(xyz[i], xc);
  }
  return xyz;
}



function sortbyz(x)
{
  var i, j, k, l, n = x.length;
  var xyz = newarr2d(n, 3), rt = newarr(D);
  var idmap = newarr(n);
  var invmap = newarr(n);

  for ( i = 0; i < n; i++ ) {
    idmap[i] = i;
    // i:         index of the output array `xyz`
    // idmap[i]:  index of the input array `x`
    // so xyz[i] --> x[ idmap[i] ];
    invmap[i] = i;
  }

  // use bubble sort
  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[i]);
  }

  for ( i = 0; i < n; i++ ) {
    // find the ith smallest z
    k = i;
    var zmin = x[ idmap[i] ][2];
    for ( j = i + 1; j < n; j++ ) {
      if ( x[ idmap[j] ][2] < zmin ) {
        k = j;
        zmin = x[ idmap[j] ][2];
      }
    }
    if ( k != i ) {
      // before
      //  xyz[i] --> x[ idmap[i] ]
      //  xyz[k] --> x[ idmap[k] ]
      // after
      //  xyz[i] --> x[ idmap[k] ]
      //  xyz[k] --> x[ idmap[i] ]
      l = idmap[i];
      idmap[i] = idmap[k];
      idmap[k] = l;
    }
  }

  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[ idmap[i] ]);
  }
  // compute the inverse map
  for ( i = 0; i < n; i++ ) {
    invmap[ idmap[i] ] = i;
  }
  return [xyz, idmap, invmap];
}



// return the z-scaling factor
function getzscale(r, zmin, zmax, ortho)
{
  if ( ortho ) {
    return 0.8;
  } else {
    var zf = (r[2] - zmin) / (zmax - zmin);
    return 0.8 + 0.2 * zf;
  }
}



/* draw all atoms in the box */
function ljdraw3d(lj, target, xin, userscale, color)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // the system dimension is L + two radii
  var scale = userscale * Math.min(width, height) / (lj.l + 1.0);
  var ortho = grab("orthographic").checked;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  var ret, xyz, idmap, invmap, zmin, zmax;
  var i, j, ic, id, jd;
  var xi, yi, scli, xj, yj, sclj;

  // draw the box
  var l = lj.l;
  var xbox = [[0,0,0], [0,0,l], [0,l,0], [0,l,l], [l,0,0], [l,0,l], [l,l,0], [l,l,l]];
  var xboxt = transform(xbox, l); // apply the rotation matrix
  ret = sortbyz(xboxt); // sort the corners by the z index
  xyz = ret[0];
  idmap = ret[1];
  invmap = ret[2];
  zmax = xyz[7][2];
  zmin = xyz[0][2];
  //console.log(invmap, xbox, xboxt);
  ctx.lineWidth = 1;
  ctx.strokeStyle = '#c0c0c0';
  var boxedges = [[0,1], [2,3], [4,5], [6,7],
                  [0,2], [1,3], [4,6], [5,7],
                  [0,4], [1,5], [2,6], [3,7]];

  for ( ic = 0; ic < boxedges.length; ic++ ) {
    i = invmap[ boxedges[ic][0] ];
    scli = scale * getzscale(xyz[i], zmin, zmax, ortho);
    xi = Math.floor(  (xyz[i][0] - lj.l * 0.5) * scli + width  * 0.5 );
    yi = Math.floor( -(xyz[i][1] - lj.l * 0.5) * scli + height * 0.5 );

    j = invmap[ boxedges[ic][1] ];
    sclj = scale * getzscale(xyz[j], zmin, zmax, ortho);
    xj = Math.floor(  (xyz[j][0] - lj.l * 0.5) * sclj + width  * 0.5 );
    yj = Math.floor( -(xyz[j][1] - lj.l * 0.5) * sclj + height * 0.5 );

    drawLine(ctx, xi, yi, xj, yj);
    //console.log(xi, yi, xj, yj, i, j, boxedges[ic], xyz[i], xyz[j], scli, sclj);
  }

  // start to draw particles
  var xt = transform(xin, lj.l); // apply the rotation matrix
  ret = sortbyz(xt); // sort particles by the z order
  xyz = ret[0];
  idmap = ret[1];
  invmap = ret[2];
  // xyz[i]           --> xt[ idmap[i] ]
  // xyz[ invmap[i] ] --> xt[ i ]
  zmax = xyz[lj.n - 1][2];
  zmin = xyz[0][2];

  // draw each particle
  for ( i = 0; i < lj.n; i++ ) {
    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    // make closer particles larger
    var zscl = getzscale(xyz[i], zmin, zmax, ortho);
    var scl = scale * zscl;
    var x = Math.floor(  (xyz[i][0] - lj.l * 0.5) * scl + width  * 0.5 );
    var y = Math.floor( -(xyz[i][1] - lj.l * 0.5) * scl + height * 0.5 );
    var rad = 0.5;
    var rz = Math.floor( rad * scl );
    zscl = getzscale(xyz[i], zmin, zmax, false);
    var color2 = darkenColor(color, zscl);
    var spotcolor = lightenColor(color2, 0.3);
    paintBall(ctx, x, y, rz, color2, spotcolor);
  }
}



