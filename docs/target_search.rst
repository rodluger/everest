.. raw:: html
  
  <h2><span id="target"></span></h2><br>

  <script>
    
  function isInt(value) {
    return !isNaN(value) && (function(x) { return (x | 0) === x; })(parseFloat(value))
  }

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
  
  var id = getParameterByName("id");
  var mission = getParameterByName("mission");
  if (id != null)
    if (mission == "k2")
      if ((id.length == 9) && isInt(id))
        document.getElementById("target").innerHTML = "EPIC " + id;
      else
        document.getElementById("target").innerHTML = "Target not found."
    else
     document.getElementById("target").innerHTML = "Target not found."
  else
    document.getElementById("target").innerHTML = "Target not found." 
  
  </script>
   
   