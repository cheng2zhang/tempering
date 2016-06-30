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
var nstepspsmd = 100; // number of steps per second for MD
var nstepspfmd = 10;  // number of steps per frame for MD
var nstepspsmc = 10000; // number of steps per second for MC
var nstepspfmc = 1000;  // number of steps per frame for MC
var simulmethod = "MC";
var weightmethod = "WL";
var weightmethod_id = 0; // WL
var mcamp = 0.2;

var thtype = "v-rescale"; // thermostat type
var vresdamp = 20.0;  // velocity-rescaling damping factor
var langdamp = 1.0;  // Langevin-dynamics damping factor
var zeta = null;
var zmass = null;

var forcescaling = false; // use force-scaling

var wl_lnf0 = 1.0;
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
var ketot = 0.0;
var kesum = 0.0;

var nstepsmc = 0;
var mctot = 0.0;
var mcacc = 0.0;

var histplot = null;
var vplot = null;

var usercode_init = null;
var usercode_tpmove = null;
var usercode_update = null;
var usercode_obj = null;



function getparams()
{
  n = get_int("n", 55);
  rho = get_float("density", 0.7);
  tpmin = get_float("tpmin", 2.5);
  tpmax = get_float("tpmax", 3.7);
  tpcnt = get_int("tpcnt", 12);
  rcdef = get_float("rcutoff", 1000.0);

  simulmethod = grab("simulmethod").value;
  weightmethod = grab("weightmethod").value;
  if ( weightmethod === "WL" ) {
    weightmethod_id = 0;
  } else if ( weightmethod === "AVE" ) {
    weightmethod_id = 1;
  } else {
    weightmethod_id = 2;
  }
  mddt = get_float("mddt", 0.002);
  nstepspsmd = get_int("nstepspersecmd", 100);
  nstepspfmd = nstepspsmd * timer_interval / 1000;
  nstepsmd = 0;

  thtype = grab("thermostattype").value;
  vresdamp = get_float("vresdamp", 20.0);
  langdamp = get_float("langdamp", 1.0);
  var nhclen = get_int("nhclen", 5);
  var nhcmass1 = get_float("nhcmass1", 1.0);
  var nhcmass2 = get_float("nhcmass2", 1.0);
  var i;
  zeta = newarr(nhclen);
  zmass = newarr(nhclen);
  for ( i = 0; i < nhclen; i++ ) {
    zeta[i] = 0.0;
    zmass[i] = ( i === 0 ) ? nhcmass1 : nhcmass2;
  }

  forcescaling = grab("forcescaling").checked;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;
  nstepsmc = 0;

  wl_lnf0 = get_float("wl_lnfinit", 1.0);
  wl_flatness = get_float("wl_flatness", 0.3);
  wl_frac = get_float("wl_frac", 0.5);

  mousescale = get_float("animationboxscale");
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

  var vpot = null;
  if ( weightmethod === "WL" ) {
    wl = new WL(0, tpcnt, 1, false, wl_lnf0, wl_flatness, wl_frac, 1.0, 0);
    vpot = wl.v;
  } else {
    wl = null;
  }

  // create a simulated tempering object
  simtemp = new SimTemp(beta, vpot);

  // attach the standard data
  for ( i = 0; i < tpcnt; i++ ) {
    simtemp.lnzref[i] = -eos[i][2] * beta[i] * n;
    simtemp.uref[i] = eos[i][0] * n;
  }
  for ( i = tpcnt - 1; i >= 0; i-- ) {
    simtemp.lnzref[i] -= simtemp.lnzref[0];
  }
  itp = 0;

  if ( usercode_init ) {
    usercode_obj = usercode_init(beta);
  }
}



/* return a string of the current simulation details */
function getsinfo()
{
  var s = "";
  if ( weightmethod === "WL" ) {
    var flatness = wl.getflatness();
    s += 'WL stage ' + wl.stage + '.<br>';
    s += '<span class="math">ln <i>f</i> </span>: ' + wl.lnf.toExponential(3) + '.<br>';
    s += 'flatness: ' + roundto(flatness * 100, 2) + '%.<br>';
  }
  return s;
}



function thermostat(tp, dt)
{
  if ( thtype === "v-rescale" ) {
    return lj.vrescale(tp, dt * vresdamp);
  } else if ( thtype === "Nose-Hoover" ) {
    return lj.nhchain(tp, dt, zeta, zmass);
  } else if (thtype === "Langevin" ) {
    return lj.langevin(tp, dt * langdamp);
  }
}



// temperature transition
function tpmove(itp, ep)
{
  if ( usercode_tpmove ) {
    if ( !usercode_obj )
      usercode_obj = usercode_init(beta);
    return usercode_tpmove(itp, ep, usercode_obj);
  } else {
    return simtemp.jump(itp, ep, 1, weightmethod_id);
  }
}



// update data
function update(itp, ep)
{
  simtemp.add(itp, ep);
  if ( usercode_update ) {
    if ( !usercode_obj )
      usercode_obj = usercode_init(beta);
    usercode_update(itp, ep, usercode_obj);
  } else {
    if ( weightmethod === "WL" ) {
      wl.add( itp );
      var ret = wl.updatelnf();
      if ( ret && simulmethod === "MC" ) {
        console.log("acc", mcacc/mctot, ", amp", mcamp);
        mctot = 1e-30;
        mcacc = 0;
      }
    }
  }
}



function domd()
{
  var istep, sinfo = "";
  var tp_thstat = (tpmax + tpmin) / 2, ke;

  for ( istep = 1; istep <= nstepspfmd; istep++ ) {
    var tp = 1 / beta[itp], fs = 1;
    if ( forcescaling ) {
      tp = tp_thstat;
      fs = beta[itp] * tp_thstat;
    }
    thermostat(tp, 0.5 * mddt);
    lj.vv_fs(mddt, fs);
    ke = thermostat(tp, 0.5 * mddt);

    // temperature transition
    var jtp = tpmove(itp, lj.epot);
    if ( jtp != itp ) {
      var s = Math.sqrt(beta[itp]/beta[jtp]), i;
      if ( !forcescaling ) {
        // scale velocities if not in force scaling
        for ( i = 0; i < n; i++ ) {
          lj.v[i][0] *= s;
          lj.v[i][1] *= s;
          lj.v[i][2] *= s;
        }
        ke *= s * s;
      }
      itp = jtp;
    }

    kesum += ke / tp;
    ketot += 1;
    update(itp, lj.epot);
  }
  if ( weightmethod === "WL" ) {
    wl.trimv();
  }
  nstepsmd += nstepspfmd;
  sinfo += "step " + nstepsmd + ".<br>";
  sinfo += "<i>E<sub>K</sub></i>/<i>E<sub>K</sub></i><sup>ref</sup>: "
         + roundto(2.0 * kesum / ketot / lj.dof, 4) + ".<br>";
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
    itp = tpmove(itp, lj.epot);
    update(itp, lj.epot);
  }
  if ( weightmethod === "WL" ) {
    wl.trimv();
  }
  nstepsmc += nstepspfmc;
  sinfo += "step: " + nstepsmc + ".<br>";
  sinfo += "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%.<br>";
  sinfo += getsinfo();
  return sinfo;
}



// normalize the histogram such that the maximum is 1.0
function normalize_hist()
{
  var n = simtemp.n, hs = newarr(n), i, hm = 0;
  for ( i = 0; i < n; i++ ) {
    hm = Math.max(hm, 1.0 * simtemp.hist[i]);
  }
  for ( i = 0; i < n; i++ ) {
    hs[i] = 1.0 * simtemp.hist[i] / hm;
  }
  return hs;
}



/* update the histogram plot */
function updatehistplot()
{
  var i;
  var dat = "Temperature,Histogram\n";

  var chs = normalize_hist();
  for ( i = 0; i < simtemp.n; i++ )
    dat += "" + beta[i] + "," + chs[i] + "\n";
  if ( histplot === null ) {
    var h = grab("animationbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      xlabel: '<small>Temperature</small>',
      ylabel: '<small>Histogram</small>',
      includeZero: true,
      drawPoints: true,
      axisLabelFontSize: 10,
      plotter: barChartPlotter,
      width: w,
      height: h
    };
    histplot = new Dygraph(grab("histplot"), dat, options);
  } else {
    histplot.updateOptions({ file: dat });
  }
}



/* update the potential energy plot */
function updatevplot()
{
  var i;
  var dat = "Temperature,potential energy,reference\n";
  for ( i = 0; i < simtemp.n; i++ ) {
    var uave = simtemp.usum[i] / (simtemp.hist[i] + 1e-30) / n;
    dat += "" + beta[i] + "," + uave + ","
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



function changescale()
{
  mousescale = get_float("animationboxscale");
  paint();
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
  updatehistplot();
  updatevplot();
}



// stop simulation and reset some data
function stopsimul()
{
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
  }
  grab("pause").innerHTML = "&#9724;";
  ketot = 1e-30;
  kesum = 0.0;
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



function startsimul()
{
  stopsimul();
  getparams();
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



var usercode_ref = new Array();

usercode_ref["stlj_WL_init"] = "" +
"// This function accepts an array of inverse temperatures\n" +
"// and return an object for future use\n" +
"function init(beta) {\n" +
"  var n = beta.length, i;\n" +
"  var varr = new Array(n);\n" +
"  var harr = new Array(n);\n" +
"  for ( i = 0; i < n; i++ ) {\n" +
"    varr[i] = 0;\n" +
"    harr[i] = 0;\n" +
"  }\n" +
"  // return the WL object\n" +
"  return { beta  : beta,\n" +
"           v     : varr,\n" +
"           hist  : harr,\n" +
"           tot   : 0,\n" +
"           lnf   : 1.0};\n" +
"}\n";

usercode_ref["stlj_WL_tpmove"] = "" +
"// This function accepts the current temprature index,\n" +
"// current potential energy, and the user object,\n" +
"// returns the new temperature index\n" +
"function tpmove(itp, ep, obj) {\n" +
"  var ntp = obj.beta.length;\n" +
"  var jtp = Math.floor(itp + 1 + rand01() * (ntp - 1)) % ntp;\n" +
"  var de = -ep * obj.beta[jtp] - obj.v[jtp]\n" +
"           +ep * obj.beta[itp] + obj.v[itp];\n" +
"  if ( de < 0 && rand01() > Math.exp(de) ) {\n" +
"    jtp = itp; // abandon the move\n" +
"  }\n" +
"  return jtp;\n" +
"}\n";

usercode_ref["stlj_WL_update"] = "" +
"// This function accepts the current temprature index,\n" +
"// current potential energy, and the user object,\n" +
"// and updates the user object\n" +
"function update(itp, ep, obj) {\n" +
"  var n = obj.beta.length;\n" +
"  obj.hist[itp] += 1;\n" +
"  obj.tot += 1;\n" +
"  obj.lnf = n / (obj.tot + 10000);\n" +
"  obj.v[itp] += obj.lnf;\n" +
"}\n";

// ace editor
var aceeditor = null;

function geteditorvalue()
{
  if ( aceeditor ) {
    return aceeditor.getValue();
  } else {
    return document.getElementById("usercode").value;
  }
}

function seteditorvalue(s)
{
  if ( aceeditor ) {
    aceeditor.setValue(s);
  } else {
    document.getElementById("usercode").value = s;
  }
}



function usercode_funcname()
{
  var set = grab("userset").value;
  if ( set.substring(0, 8) === "userset_" ) {
    set = set.substring(8);
  }
  var func = grab("userfunc").value;
  if ( func.substring(0, 9) === "userfunc_" ) {
    func = func.substring(9);
  }
  return "stlj_" + set + "_" + func;
}

// compile the user-defined function
function usercode_compile()
{
  stopsimul();
  var s = geteditorvalue();
  var p0 = s.indexOf("function(");
  var p1 = s.indexOf("(", p0) + 1;
  var p9 = s.indexOf(")", p1);
  var id0 = s.indexOf("{", p9);
  var id1 = s.lastIndexOf("}");
  // extracting the function body including the braces
  var funcbody = s.substring(id0, id1 + 1);
  var funcid = grab("userfunc").value;
  var npar = (funcid === "userfunc_init") ? 1 : 3;
  var thefunc, var1, var2, var3;
  if ( npar === 1 ) {
    var1 = s.substring(p1, p9).trim();
    thefunc = new Function(var1, funcbody);
  } else if ( npar === 3 ) {
    var p2 = s.indexOf(",", p1);
    var p3 = s.indexOf(",", p2 + 1);
    var1 = s.substring(p1, p2).trim();
    var2 = s.substring(p2 + 1, p3).trim();
    var3 = s.substring(p3 + 1, p9).trim();
    thefunc = new Function(var1, var2, var3, funcbody);
  }
  if ( funcid === "userfunc_init" ) {
    usercode_init = thefunc;
  } else if ( funcid === "userfunc_tpmove" ) {
    usercode_tpmove = thefunc;
  } else if ( funcid === "userfunc_update" ) {
    usercode_update = thefunc;
  }
}


function usercode_changeset()
{
  //usercode_init = ;
  //usercode_tpmove = ;
  //usercode_update = ;
}

// reset the user's code
function usercode_reset(what)
{
  localStorage.removeItem(what);
  //seteditorvalue( );
}

// save the user's code
function usercode_save()
{
  var funcname = usercode_funcname();
  var usercode = geteditorvalue();
  localStorage.setItem(funcname, usercode);
}

// load the user's code
function usercode_load(what)
{
  // remind the user to save the code
  var funcname = usercode_funcname();
  //usercode_load(funcname);
  var code = usercode_ref[funcname];
  if ( code ) {
    seteditorvalue( code );
  }
  //var usercode = localStorage.getItem(what);
  //if ( usercode ) {
  //  seteditorvalue( usercode );
  //  //alert("loading " + what + " " + usercode);
  //}
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
  var htab = 450;

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

