/* Handle web interface */



"use strict";



// parameters
var n = 55;
var rho = 0.5;
var tpmin = 2.5;
var tpmax = 3.7;
var tpcnt = 12;
var rcdef = 1000.0;
var mddt = 0.002;
var thdt = 0.02;
var nstepspsmd = 100; // number of steps per second for MD
var nstepspfmd = 10;  // number of steps per frame for MD
var nstepspsmc = 10000; // number of steps per second for MC
var nstepspfmc = 1000;  // number of steps per frame for MC
var simulmethod = "MC";
var mcamp = 0.2;

var wl_lnf0 = 0.01;
var wl_flatness = 0.3;
var wl_frac = 0.5;

var timer_interval = 100; // in milliseconds

var lj = null;
var simtemp = null;
var wl = null;
var ljtimer = null;

var beta = null; // temperature array
var itp = 0; // current temperature index

var nstepsmd = 0;
var hmctot = 0.0;
var hmcacc = 0.0;

var nstepsmc = 0;
var mctot = 0.0;
var mcacc = 0.0;

var histplot = null;
var vplot = null;

var color = "#407f10"; // color of ball



function getparams()
{
  n = get_int("n", 55);
  rho = get_float("density", 0.7);
  tpmin = get_float("tpmin", 2.5);
  tpmax = get_float("tpmax", 3.7);
  tpcnt = get_int("tpcnt", 12);
  rcdef = get_float("rcutoff", 1000.0);

  simulmethod = grab("simulmethod").value;
  mddt = get_float("mddt", 0.002);
  thdt = get_float("thermostatdt", 0.01);
  nstepspsmd = get_int("nstepspersecmd", 100);
  nstepspfmd = nstepspsmd * timer_interval / 1000;
  nstepsmd = 0;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;
  nstepsmc = 0;

  wl_lnf0 = get_float("wl_lnfinit", 0.01);
  wl_flatness = get_float("wl_flatness", 0.3);
  wl_frac = get_float("wl_frac", 0.5);

  mousescale = get_float("animationboxscale");
}



function changescale()
{
  mousescale = get_float("animationboxscale");
  paint();
}



/* return a string of the current simulation details */
function getsinfo()
{
  var s = "";
  var flatness = wl.getflatness();
  s += 'WL stage ' + wl.stage + '.<br>';
  s += '<span class="math">ln <i>f</i> </span>: ' + wl.lnf.toExponential(3) + '.<br>';
  s += 'flatness: ' + roundto(flatness * 100, 2) + '%.<br>';
  return s;
}



function domd()
{
  var istep, sinfo = "";

  var tp = 2.0;
  for ( istep = 1; istep <= nstepspfmd; istep++ ) {
    lj.vv(mddt);
    lj.vrescale(tp, thdt);

    //wl.add( lj.csize );
    wl.updatelnf();
  }
  wl.trimv();
  nstepsmd += nstepspfmd;
  sinfo += "step " + nstepsmd + ".<br>";
  //sinfo += "hmcacc: " + roundto(100.0 * hmcacc / hmctot, 2) + "%.<br>";
  sinfo += getsinfo();
  return sinfo;
}



function domc()
{
  var istep, sinfo = "";

  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    // do a step of Metropolis algorithm
    mctot += 1.0;
    mcacc += lj.metro(mcamp, beta[itp]);
    // temperature transition
    itp = simtemp.jump(itp, lj.epot, 1);
    simtemp.add(itp, lj.epot);
    wl.add( itp );
    if ( wl.updatelnf() ) {
      console.log("acc", mcacc/mctot, ", amp", mcamp);
      mctot = 1e-30;
      mcacc = 0;
    }
  }
  wl.trimv();
  nstepsmc += nstepspfmc;
  sinfo += "step: " + nstepsmc + ".<br>";
  sinfo += "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%.<br>";
  sinfo += getsinfo();
  return sinfo;
}



// normalize the histogram such that the maximum is 1.0
function normalize_hist(arr, n)
{
  var hs = newarr(n);
  var i, m = 0;
  for ( i = 0; i < n; i++ ) {
    m = Math.max(m, 1.0 * arr[i]);
  }
  for ( i = 0; i < n; i++ ) {
    hs[i] = 1.0 * arr[i] / m;
  }
  return hs;
}



/* update the histogram plot */
function updatehistplot(wl)
{
  var i;
  var dat = "Temperature,Histogram (all time),Histogram (this stage)\n";

  var chhs = normalize_hist(wl.hh, wl.n);
  var chs  = normalize_hist(wl.h, wl.n);
  for ( i = 0; i < wl.n; i++ )
    dat += "" + (beta[i]) + "," + chhs[i] + "," + chs[i] + "\n";
  if ( histplot === null ) {
    var h = grab("animationbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      xlabel: '<small>Temperature</small>',
      ylabel: '<small>Histogram</small>',
      includeZero: true,
      drawPoints: true,
      axisLabelFontSize: 10,
      pointSize: 2,
      xRangePad: 2,
      plotter: barChartPlotter,
      width: w,
      height: h
    };
    histplot = new Dygraph(grab("histplot"), dat, options);
  } else {
    histplot.updateOptions({ file: dat });
  }
}



/* update the cluster potential plot */
function updatevplot(wl)
{
  var i;
  var dat = "Temperature,potential energy,reference\n";
  for ( i = 0; i < wl.n; i++ ) {
    var uave = simtemp.usum[i] / (simtemp.hist[i] + 1e-30) / n;
    dat += "" + (beta[i]) + "," + uave + ","
        + (simtemp.uref[i]/n) + "\n";
  }
  if ( vplot === null ) {
    var h = grab("animationbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      //title: 'Adaptive potential',
      xlabel: '<small><i>&beta;<sub>i</sub></i></small>',
      ylabel: '<small><i>U<sub>i</sub></i></small>',
      includeZero: false,
      drawPoints: true,
      axisLabelFontSize: 10,
      pointSize: 2,
      xRangePad: 2,
      width: w,
      height: h
    };
    vplot = new Dygraph(grab("vplot"), dat, options);
  } else {
    vplot.updateOptions({ file: dat });
  }
}



function paint()
{
  if ( !lj ) return;
  var r = Math.floor( 255 * itp / tpcnt );
  var color = "rgb(" + r + ", 0, " + (255 - r) + ")";
  ljdraw3d(lj, "animationbox", lj.x, mousescale, color);
}



function pulse()
{
  var sinfo;

  // the following parameter might have been changed
  if ( wl ) {
    wl.flatness = get_float("wl_flatness", 0.3);
    wl.frac = get_float("wl_frac", 0.5);
  }

  if ( simulmethod === "MD" ) {
    sinfo = domd();
  } else if ( simulmethod === "MC" ) {
    sinfo = domc();
  }
  grab("sinfo").innerHTML = sinfo;

  paint();
  updatehistplot(wl);
  updatevplot(wl);
}



function stopsimul()
{
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
  }
  grab("pause").innerHTML = "&#9724;";
  hmctot = 1e-30;
  hmcacc = 0.0;
  mctot = 1e-30;
  mcacc = 0.0;
  munit(viewmat);
}



function pausesimul()
{
  if ( !lj ) return;
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
    grab("pause").className = "glyphicon glyphicon-play";
  } else {
    ljtimer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").className = "glyphicon glyphicon-pause";
  }
}


function init_simtemp()
{
  var i, bmin = 1.0/tpmax, bmax = 1.0/tpmin;
  var eos = new Array(tpcnt);

  // initialize the temperature array
  beta = new Array(tpcnt);
  for ( i = 0; i < tpcnt; i++ ) {
    beta[i] = bmin + (bmax - bmin) * i / (tpcnt - 1);
    eos[i] = lj_eos3dPVEhBH(rho, 1/beta[i]);
  }

  // create a simulated tempering object
  simtemp = new SimTemp(beta, wl.v);

  // attach the standard data
  for ( i = 0; i < tpcnt; i++ ) {
    simtemp.lnzref[i] = -eos[i][2] * beta[i] * n;
    simtemp.uref[i] = eos[i][0] * n;
  }
  for ( i = tpcnt - 1; i >= 0; i-- ) {
    simtemp.lnzref[i] -= simtemp.lnzref[0];
  }
  itp = 0;
}



function startsimul()
{
  stopsimul();
  getparams();
  wl = new WL(0, tpcnt - 1, 1, false, wl_lnf0, wl_flatness, wl_frac, 1.0, 0);
  lj = new LJ(n, D, rho, rcdef);
  lj.force();
  init_simtemp();
  installmouse("animationbox", "animationboxscale");
  ljtimer = setInterval(
    function(){ pulse(); },
    timer_interval);
}



// on a mouse click
function pausesimul2()
{
  // skip a mouse-move
  if ( mousemoved > 0 ) {
    return;
  }
  if ( !lj ) {
    startsimul();
  } else if ( mousemoved === 0 ) {
    pausesimul();
  }
}



function changeparams()
{
  if ( ljtimer !== null ) {
    startsimul();
  }
}



function showtab(who)
{
  who = grab(who);
  var par = who.parentNode;
  var c = par.childNodes;
  var i, iwho, k = 0;

  // arrange the tabs
  for ( i = 0; i < c.length; i++ ) {
    if ( c[i].className === "params-panel" ) {
      if ( c[i] !== who ) {
        c[i].style.zIndex = k;
        k += 1;
      } else {
        iwho = k;
      }
    }
  }
  who.style.zIndex = k;

  // arrange the clickable tab titles
  k += 1;
  var pt = grab("tabheadings");
  pt.style.zIndex = k;
  var ct = pt.childNodes, ik = 0;
  for ( i = 0; i < ct.length; i++ ) {
    if ( ct[i].tagName ) {
      if ( ik === iwho ) {
        ct[i].style.fontWeight = "bold";
        ct[i].style.borderTop = "2px solid #c0c0d0";
      } else {
        ct[i].style.fontWeight = "normal";
        ct[i].style.borderTop = "0px solid #e0e0f0";
      }
      ik += 1;
    }
  }
}



function resizecontainer(a)
{
  var canvas = grab("animationbox");
  var ctx = canvas.getContext("2d");
  var w, h;
  if ( a === null || a === undefined ) {
    w = canvas.width;
    h = canvas.height;
  } else {
    a = parseInt( grab(a).value );
    w = h = a;
    canvas.width = w;
    canvas.height = h;
  }
  ctx.font = "24px Verdana";
  ctx.fillText("Click to start", w/2-40, h/2-10);

  var hsbar = 30; // height of the global scaling bar
  var hcbar = 40; // height of the control bar
  var htbar = 30; // height of the tabs bar
  var wr = h*3/4; // width of the plots
  var wtab = 700; // width of the tabs
  var htab = 280;

  grab("simulbox").style.width = "" + w + "px";
  grab("simulbox").style.height = "" + h + "px";
  grab("simulbox").style.top = "" + hsbar + "px";
  grab("controlbox").style.top = "" + (h + hsbar) + "px";
  grab("animationboxscale").style.width = "" + (w - 100) + "px";
  histplot = null;
  grab("histplot").style.left = "" + w + "px";
  grab("histplot").style.width = "" + wr + "px";
  grab("vplot").style.top = "" + hcbar + "px";
  grab("histplot").style.height = "" + h/2 + "px";
  vplot = null;
  grab("vplot").style.left = "" + w + "px";
  grab("vplot").style.width = "" + wr + "px";
  grab("vplot").style.top = "" + (h/2 + hcbar) + "px";
  grab("vplot").style.height = "" + h/2 + "px";
  grab("tabsrow").style.top = "" + (h + hsbar + hcbar) + "px";
  grab("tabsrow").style.width = "" + wtab + "px";
  var c = grab("container").childNodes;
  var i;
  /* tabs */
  for ( i = 0; i < c.length; i++ ) {
    if ( c[i].className === "params-panel" ) {
      c[i].style.top = "" + (h + hsbar + hcbar + htbar) + "px";
      c[i].style.width = "" + (wtab - 20) + "px";
      c[i].style.height = "" + htab + "px";
    }
  }
  grab("sinfo").style.top = "" + (h + hsbar + hcbar + htbar) + "px";
  grab("sinfo").style.left = "" + (wtab + 10) + "px";
  grab("sinfo").style.width = "" + (w + wr - wtab - 20) + "px";
  grab("container").style.height = "" + (h + hsbar + hcbar + htbar + htab) + "px";
  grab("container").style.width = "" + (w + wr) + "px";
}



function init()
{
  resizecontainer();
  showtab("system-params");
}

