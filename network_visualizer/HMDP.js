var labelType, useGradients, nativeTextSupport, animate;

(function() {
  var ua = navigator.userAgent,
      iStuff = ua.match(/iPhone/i) || ua.match(/iPad/i),
      typeOfCanvas = typeof HTMLCanvasElement,
      nativeCanvasSupport = (typeOfCanvas == 'object' || typeOfCanvas == 'function'),
      textSupport = nativeCanvasSupport 
        && (typeof document.createElement('canvas').getContext('2d').fillText == 'function');
  //I'm setting this based on the fact that ExCanvas provides text support for IE
  //and that as of today iPhone/iPad current text support is lame
  labelType = (!nativeCanvasSupport || (textSupport && !iStuff))? 'Native' : 'HTML';
  nativeTextSupport = labelType == 'Native';
  useGradients = nativeCanvasSupport;
  animate = !(iStuff || !nativeCanvasSupport);
})();

var Log = {
  elem: false,
  write: function(text){
    if (!this.elem) 
      this.elem = document.getElementById('log');
    this.elem.innerHTML = text;
    this.elem.style.left = (500 - this.elem.offsetWidth / 2) + 'px';
  }
};


function init(){
  // init data
  // end
  // init ForceDirected
  var fd = new $jit.ForceDirected({
    //id of the visualization container
    injectInto: 'infovis',
    //Enable zooming and panning
    //by scrolling and DnD
    Navigation: {
      enable: true,
      //Enable panning events only if we're dragging the empty
      //canvas (and not a node).
      panning: 'avoid nodes',
      zooming: 10 //zoom speed. higher is more sensible
    },
    // Change node and edge styles such as
    // color and width.
    // These properties are also set per node
    // with dollar prefixed data-properties in the
    // JSON structure.
    Node: {
      overridable: true
    },
    Edge: {
      overridable: true,
      color: '#23A4FF',
      lineWidth: 0.4
    },
    //Native canvas text styling
    Label: {
      type: labelType, //Native or HTML
      size: 10,
      style: 'bold'
    },
    //Add Tips
    Tips: {
      enable: true,
      onShow: function(tip, node) {
        //count connections
        var count = 0;
        node.eachAdjacency(function() { count++; });
        //display node info in tooltip
        tip.innerHTML = "<div class=\"tip-title\">" + node.name + "</div>"
          + "<div class=\"tip-text\"><b>Associations:</b> " + count + "</div>";
      }
    },
    // Add node events
    Events: {
      enable: true,
      type: 'Native',
      //Change cursor style when hovering a node
      onMouseEnter: function() {
        fd.canvas.getElement().style.cursor = 'move';
      },
      onMouseLeave: function() {
        fd.canvas.getElement().style.cursor = '';
      },
      //Update node positions when dragged
      onDragMove: function(node, eventInfo, e) {
          var pos = eventInfo.getPos();
          node.pos.setc(pos.x, pos.y);
          fd.plot();
      },
      //Implement the same handler for touchscreens
      onTouchMove: function(node, eventInfo, e) {
        $jit.util.event.stop(e); //stop default touchmove event
        this.onDragMove(node, eventInfo, e);
      },
      //Add also a click handler to nodes
      onClick: function(node) {
        if(!node) return;
        // make the node the selected node for searching
	document.getElementById("input_name").value = node.name;
        document.getElementById("input_id").value = node.id;
        // Build the right column relations list.
        // This is done by traversing the clicked node connections.
        var html =
		"<form id=\"selectedelement\" method=\"post\" action=\"\"><h4><a href=\"#\" id=\"elementidname\" onclick=\"traitSelectedChange(\'elementidname\');\">" + node.name +
		"</a></h4><input type=\"hidden\" name=\"elementid\" id=\"elementid\" value=\"" + node.id + "\"></form><b> Associations:</b><ul><li>",
            list = [];
        node.eachAdjacency(function(adj){
          list.push(adj.nodeTo.name);
        });
        //append connections information
        $jit.id('infoclick').innerHTML = html + list.join("</li><li>") + "</li></ul>";
        document.getElementById("elementpval").focus();
      }
    },
    //Number of iterations for the FD algorithm
    iterations: 100,
    //Edge length
    levelDistance: 130,
    // Add text to the labels. This method is only triggered
    // on label creation and only for DOM labels (not native canvas ones).
    onCreateLabel: function(domElement, node){
      domElement.innerHTML = node.name;
      var style = domElement.style;
      style.fontSize = "0.8em";
      style.color = "#ddd";
    },
    // Change node styles when DOM labels are placed
    // or moved.
    onPlaceLabel: function(domElement, node){
      var style = domElement.style;
      var left = parseInt(style.left);
      var top = parseInt(style.top);
      var w = domElement.offsetWidth;
      style.left = (left - w / 2) + 'px';
      style.top = (top + 10) + 'px';
      style.display = '';
    }
  });
  // load JSON data.
  fd.loadJSON(json);
  // compute positions incrementally and animate.
  fd.computeIncremental({
    iter: 40,
    property: 'end',
    onStep: function(perc){
      Log.write(perc + '% loaded...');
    },
    onComplete: function(){
      Log.write('');
      fd.animate({
        modes: ['linear'],
        transition: $jit.Trans.Elastic.easeOut,
        duration: 50
      });
    }
  });
  // end
}

function getClickedElement(which_elementid)
{
  var elementpval = document.getElementById("elementpval").value;
  var elementid = document.getElementById(which_elementid).value;

  var newForm = document.createElement("form");
  newForm.action = "HMDP.php";
  newForm.method = "post";

  var newItem1 = document.createElement("input");
  newItem1.type = "hidden";
  newItem1.name = "elementpval";
  newItem1.value = elementpval;

  var newItem2 = document.createElement("input");
  newItem2.type = "hidden";
  newItem2.name = "elementid";
  newItem2.value = elementid;

  newForm.appendChild(newItem1);
  newForm.appendChild(newItem2);
  document.getElementById("hidden_form_container").appendChild(newForm);
  newForm.submit();
}

function traitListChange()
{
  document.getElementById("input_name").value = document.getElementById("traitlistid").options[document.getElementById("traitlistid").selectedIndex].innerHTML;
  document.getElementById("input_id").value = document.getElementById("traitlistid").options[document.getElementById("traitlistid").selectedIndex].value;
}

function traitSelectedChange(name)
{
  id = "";
  if(name == "selectedidname"){
    id = "selectedid";
  }
  else if(name == "elementidname"){
    id = "elementid";
  }
  document.getElementById("input_name").value = document.getElementById(name).innerHTML;
  document.getElementById("input_id").value = document.getElementById(id).value;
  document.getElementById("elementpval").focus();
}
