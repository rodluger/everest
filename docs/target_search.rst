.. rubric: Search

.. raw:: html
   
   <h1><span id="EPIC"></span></h1><br>
   
   <script>
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
   var epic = getParameterByName("epic");
   if (epic != null)
     document.getElementById("EPIC").innerHTML = "EPIC " + epic;
   </script>
   
   