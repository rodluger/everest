.. raw:: html
  
  <head>
  <script src="https://code.jquery.com/jquery-1.10.2.js"></script>
  </head>
    
  <!-- Title --!>
  <h2 style="display: inline;"><span id="target-name">Loading...</span></h2>
  <br><br>
  
  <!-- Info --!>
  <div class="wy-table-responsive" id = "target-table">
  <table border="1" class="docutils">
  <colgroup>
  <col width="19%">
  <col width="81%">
  </colgroup>
  <tbody valign="top">
  <tr class="row-odd"><td style="font-weight:bold;">Mission</td>
  <td id="target-mission"></td>
  </tr>
  <tr class="row-even"><td style="font-weight:bold;" id="target-season-label">Season</td>
  <td id="target-season"></td>
  </tr>
  <tr class="row-odd"><td style="font-weight:bold;">Cadence</td>
  <td id="target-cadence">Long</td>
  </tr>
  <tr class="row-even"><td style="font-weight:bold;">Magnitude</td>
  <td id="target-mag"></td>
  </tr>
  <tr class="row-odd"><td style="font-weight:bold;">Saturated</td>
  <td id="target-sat"></td>
  </tr>
  <tr class="row-even"><td style="font-weight:bold;">Raw 6hr CDPP</td>
  <td id="target-cdppr"></td>
  </tr>
  <tr class="row-odd"><td style="font-weight:bold;">EVEREST 6hr CDPP</td>
  <td id="target-cdpp"></td>
  </tr>
  <tr class="row-even"><td style="font-weight:bold;">Improvement factor</td>
  <td id="target-impr"></td>
  <tr class="row-odd"><td style="font-weight:bold;">Links</td>
  <td id="target-links">
    <span id="target-links-lc">
      [<a id="target-fits-link-lc" href = "#">LC FITS</a>]
      &nbsp;
      [<a id="target-dvs-link-lc" href = "#">LC DVS</a>]
    </span>
    <span id="target-links-sc" style="display:none;">
      &nbsp;
      [<a id="target-fits-link-sc" href = "#">SC FITS</a>]
      &nbsp;
      [<a id="target-dvs-link-sc" href = "#">SC DVS</a>]
    </span>
  </td>
  </tr>
  </tbody>
  </table>
  </div>
  
  <!-- The magical stuff --!>
  <script>
    
    function isInt(value) {
      return !isNaN(value) && (function(x) { return (x | 0) === x; })(parseFloat(value))
    }
    
    // Reads query strings
    function getParameterByName(name, url) {
      if (!url) {
        url = window.location.href;
      }
      name = name.replace(/[\[\]]/g, "\\$&");
      var regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)"),
        results = regex.exec(url);
      if (!results) return null;
      if (!results[2]) return '';
      return decodeURIComponent(results[2].replace(/\+/g, " "));
    }
    
    // Read query strings
    var id = getParameterByName("id");
    var mission = getParameterByName("mission");
    
    // Constants
    var mast_root = "http://archive.stsci.edu/missions/hlsp/everest/v2/";
    var version = "2.0";
    
    // Initialize
    document.title = "Loading...";
    document.getElementById("target-table").style.display = 'none';
    
    // Mission-specific
    if (id != null) {
      if (mission == "k2") {
        if ((id.length == 9) && isInt(id)) {
          
          // Check if short cadence exists
          $.get("sc.csv")
            .done(function( data ) {
              if (data.indexOf(id) != -1) {
                // Short cadence is available
                document.getElementById("target-cadence").innerHTML = ("Long | Short");
                document.getElementById("target-links-sc").style.display = 'inline';
              }
            });
          
          // Figure out the campaign number
          var campaigns = ["c00", "c01", "c02", "c03", "c04", "c05", "c06", "c07", "c08"];
          var done = 0;
          for (var i = 0; i < campaigns.length; i++) {
            $.get(campaigns[i] + ".csv")
              .done(function( data ) {
                var start = data.indexOf(id);
                if (start != -1) {
                
                  // We found the target; get the data
                  var campaign = data.substr(0,3);
                  if (campaign.substr(1,2) < 10)
                    var campaign_int = campaign.substr(2,1)
                  else
                    var campaign_int = campaign.substr(1,2)
                  var stop = data.indexOf("\n", start);
                  var info = data.slice(start, stop).split(',');
                  var mag = info[1];
                  var cdppr = info[2];
                  var cdpp = info[3];
                  if (info[4] == "1")
                    var saturated = "Yes";
                  else
                    var saturated = "No";

                  // Set the name
                  document.getElementById("target-name").innerHTML = "EPIC " + id;
                  document.title = "EPIC " + id;
          
                  // Path to MAST
                  var path = mast_root + campaign + "/" + id.substr(0,4) + "00000" + "/" + id.substr(4,5) + "/";
          
                  // Path to files
                  var fits = path + "hlsp_everest_k2_llc_" + id + "-" + campaign + "_kepler_v" + version + "_lc.fits";
                  var fits_sc = path + "hlsp_everest_k2_llc_" + id + "-" + campaign + "_kepler_v" + version + "_sc.fits";
                  var dvs = path + "hlsp_everest_k2_llc_" + id + "-" + campaign + "_kepler_v" + version + "_dvs.pdf";
                  var dvs_sc = path + "hlsp_everest_k2_llc_" + id + "-" + campaign + "_kepler_v" + version + "_dvs_sc.fits";
          
                  // Set the info
                  document.getElementById("target-mission").innerHTML = ("K2");
                  document.getElementById("target-season").innerHTML = (campaign_int);
                  document.getElementById("target-season-label").innerHTML = ("Campaign");
                  document.getElementById("target-mag").innerHTML = (mag);
                  document.getElementById("target-sat").innerHTML = (saturated);
                  document.getElementById("target-cdppr").innerHTML = (cdppr + " ppm");
                  document.getElementById("target-cdpp").innerHTML = (cdpp + " ppm");
                  document.getElementById("target-impr").innerHTML = ( (cdppr / cdpp).toFixed(2) );
                  document.getElementById("target-fits-link-lc").setAttribute('href', fits);
                  document.getElementById("target-fits-link-sc").setAttribute('href', fits_sc);
                  document.getElementById("target-dvs-link-lc").setAttribute('href', dvs);
                  document.getElementById("target-dvs-link-sc").setAttribute('href', dvs_sc);
                  
                  // Make visible
                  document.getElementById("target-table").style.display = 'inline';
                }
                
                // Check if we've gone through all campaigns. If still no match, throw error
                done++;
                if (done == campaigns.length - 1) {
                  if (document.title == "Loading...") {
                    document.getElementById("target-name").innerHTML = "Target not found.";
                    document.title = "Target not found.";
                  }
                }
            });
          }
        } else {
          document.getElementById("target-name").innerHTML = "Target not found.";
          document.title = "Target not found.";
        }
      } else {
        document.getElementById("target-name").innerHTML = "Target not found.";
        document.title = "Target not found.";
      }
    } else {
      document.getElementById("target-name").innerHTML = "Target not found.";
      document.title = "Target not found.";
    }
  </script>

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-3', 'auto');
    ga('send', 'pageview');

  </script>